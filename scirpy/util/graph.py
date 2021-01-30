from scanpy import logging
import igraph as ig
import numpy as np
from .._compat import Literal
from scipy.sparse import spmatrix, csr_matrix


def igraph_from_sparse_matrix(
    matrix: spmatrix,
    *,
    matrix_type: Literal["connectivity", "distance"] = "distance",
    max_value: float = None,
) -> ig.Graph:
    """
    Get an igraph object from an adjacency or distance matrix.

    Parameters:
    -----------
    matrix
        A sparse matrix that represents the connectivity or distance matrix for the graph.
        Zero-entries mean "no edge between the two nodes".
    matrix_type
        Whether the `sparse_matrix` represents connectivities (higher value = smaller distance)
        or distances (higher value = higher distance). Distance matrices will be
        converted into connectivities. A connectivity matrix is also known as
        weighted adjacency matrix.
    max_value
        When converting distances to connectivities, this will be considered the
        maximum distance. This defaults to `numpy.max(sparse_matrix)`.

    Returns
    -------
    igraph object
    """
    matrix = matrix.tocsr()

    if matrix_type == "distance":
        matrix = _distance_to_connectivity(matrix, max_value=max_value)

    return _get_igraph_from_adjacency(matrix)


def _distance_to_connectivity(distances: csr_matrix, *, max_value=None) -> csr_matrix:
    """Get a weighted adjacency matrix from a distance matrix.

    A distance of 1 (in the sparse matrix) corresponds to an actual distance of 0.
    An actual distance of 0 corresponds to a connectivity of 1.

    A distance of 0 (in the sparse matrix) corresponds to an actual distance of
    infinity. An actual distance of infinity corresponds to a connectivity of 0.

    Parameters
    ----------
    distances
        sparse distance matrix
    max_value
        The max_value is used to normalize the distances, i.e. distances
        are divided by this value. If not specified it will
        be the max. of the input matrix.
    """
    if not isinstance(distances, csr_matrix):
        raise ValueError("Distance matrix must be in CSR format.")

    if max_value is None:
        max_value = np.max(distances)

    connectivities = distances.copy()
    d = connectivities.data - 1

    # structure of the matrix stays the same, we can safely change the data only
    connectivities.data = (max_value - d) / max_value
    connectivities.eliminate_zeros()

    return connectivities


def _get_igraph_from_adjacency(adj: csr_matrix):
    """Get an undirected igraph graph from adjacency matrix.
    Better than Graph.Adjacency for sparse matrices.

    Parameters
    ----------
    adj
        sparse, weighted, symmetrical adjacency matrix.
    """
    sources, targets = adj.nonzero()
    weights = adj[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    if isinstance(weights, csr_matrix):
        # this is the case when len(sources) == len(targets) == 0, see #236
        weights = weights.toarray()

    g = ig.Graph(directed=False)
    g.add_vertices(adj.shape[0])  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))

    g.es["weight"] = weights

    if g.vcount() != adj.shape[0]:
        logging.warning(
            f"The constructed graph has only {g.vcount()} nodes. "
            "Your adjacency matrix contained redundant nodes."
        )  # type: ignore

    # since we start from a symmetrical matrix, and the graph is undirected,
    # it is fine to take either of the two edges when simplifying.
    # g.simplify(combine_edges="first")

    return g


def _get_sparse_from_igraph(graph, weight_attr=None):
    # TODO remove unless I end up using this for testing.
    edges = graph.get_edgelist()
    if weight_attr is None:
        weights = [1] * len(edges)
    else:
        weights = graph.es[weight_attr]
    shape = graph.vcount()
    shape = (shape, shape)
    if len(edges) > 0:
        return csr_matrix((weights, zip(*edges)), shape=shape)
    else:
        return csr_matrix(shape)


def layout_components(
    graph: ig.Graph,
    component_layout: str = "fr",
    arrange_boxes: Literal["size", "rpack", "squarify"] = "squarify",
    pad_x: float = 1.0,
    pad_y: float = 1.0,
) -> np.ndarray:
    """
    Compute a graph layout by layouting all connected components individually.

    Adapted from https://stackoverflow.com/questions/53120739/lots-of-edges-on-a-graph-plot-in-python

    Parameters
    ----------
    graph
        The igraph object to plot.
    component_layout
        Layout function used to layout individual components.
        Can be anything that can be passed to `igraph.Graph.layout`
    arrange_boxes
        How to arrange the individual components. Can be "size"
        to arange them by the component size, or "rpack" to pack them as densly
        as possible, or "squarify" to arrange them using a treemap algorithm.
    pad_x
        Padding between subgraphs in the x dimension.
    pad_y
        Padding between subgraphs in the y dimension.

    Returns
    -------
    pos
        n_nodes x dim array containing the layout coordinates

    """
    # assign the original vertex id, it will otherwise get lost by decomposition
    for i, v in enumerate(graph.vs):
        v["id"] = i
    components = np.array(graph.decompose(mode="weak"))
    component_sizes = np.array([component.vcount() for component in components])
    order = np.argsort(component_sizes)
    components = components[order]
    component_sizes = component_sizes[order]
    vertex_ids = [v["id"] for comp in components for v in comp.vs]
    vertex_sorter = np.argsort(vertex_ids)

    bbox_fun = {"rpack": _bbox_rpack, "size": _bbox_sorted, "squarify": _bbox_squarify}[
        arrange_boxes
    ]
    bboxes = bbox_fun(component_sizes, pad_x, pad_y)

    component_layouts = [
        _layout_component(component, bbox, component_layout)
        for component, bbox in zip(components, bboxes)
    ]
    # get vertexes back into their original order
    coords = np.vstack(component_layouts)[vertex_sorter, :]
    return coords


def _bbox_rpack(component_sizes, pad_x=1.0, pad_y=1.0):
    """Compute bounding boxes for individual components
    by arranging them as densly as possible.

    Depends on `rectangle-packer`.
    """
    try:
        import rpack
    except ImportError:
        raise ImportError(
            "Using the 'components layout' requires the installation of "
            "the `rectangle-packer`. You can install it with "
            "`pip install rectangle-packer`."
        )

    dimensions = [_get_bbox_dimensions(n, power=0.8) for n in component_sizes]
    # rpack only works on integers; sizes should be in descending order
    dimensions = [
        (int(width + pad_x), int(height + pad_y))
        for (width, height) in dimensions[::-1]
    ]
    origins = rpack.pack(dimensions)
    outer_dimensions = rpack.enclosing_size(dimensions, origins)
    aspect_ratio = outer_dimensions[0] / outer_dimensions[1]
    if aspect_ratio > 1:
        scale_width, scale_height = 1, aspect_ratio
    else:
        scale_width, scale_height = aspect_ratio, 1
    bboxes = [
        (
            x,
            y,
            width * scale_width - pad_x,
            height * scale_height - pad_y,
        )
        for (x, y), (width, height) in zip(origins, dimensions)
    ]
    return bboxes[::-1]


def _bbox_squarify(component_sizes, pad_x=10, pad_y=10):
    """Arrange bounding boxes using the `squarify` implementation for treemaps"""
    try:
        import squarify
    except ImportError:
        raise ImportError(
            "Using the 'components layout' requires the installation"
            "of the `squarify` package. You can install it with "
            "`pip install squarify`"
        )
    order = np.argsort(-component_sizes)
    undo_order = np.argsort(order)
    component_sizes = component_sizes[order]
    component_sizes = squarify.normalize_sizes(component_sizes, 100, 100)
    rects = squarify.padded_squarify(component_sizes, 0, 0, 100, 100)

    bboxes = []
    for r in rects:
        width = r["dx"]
        height = r["dy"]
        offset_x = r["x"]
        offset_y = r["y"]
        delta = abs(width - height)
        if width > height:
            width = height
            offset_x += delta / 2
        else:
            height = width
            offset_y += delta / 2
        bboxes.append((offset_x, offset_y, width - pad_x, height - pad_y))

    return [bboxes[i] for i in undo_order]


def _bbox_sorted(component_sizes, pad_x=1.0, pad_y=1.0):
    """Compute bounding boxes for individual components
    by arranging them by component size"""
    bboxes = []
    x, y = (0, 0)
    current_n = 1
    for n in component_sizes:
        width, height = _get_bbox_dimensions(n, power=0.8)

        if not n == current_n:  # create a "new line"
            x = 0  # reset x
            y += height + pad_y  # shift y up
            current_n = n

        bbox = x, y, width, height
        bboxes.append(bbox)
        x += width + pad_x  # shift x down the line
    return bboxes


def _get_bbox_dimensions(n, power=0.5):
    # return (np.sqrt(n), np.sqrt(n))
    return (n ** power, n ** power)


def _layout_component(component, bbox, component_layout_func):
    """Compute layout for an individual component"""
    layout = component.layout(component_layout_func)
    rescaled_pos = _rescale_layout(np.array(layout.coords), bbox)
    return rescaled_pos


def _rescale_layout(coords, bbox):
    """Transpose the layout of a component into its bounding box"""
    min_x, min_y = np.min(coords, axis=0)
    max_x, max_y = np.max(coords, axis=0)

    if not min_x == max_x:
        delta_x = max_x - min_x
    else:  # graph probably only has a single node
        delta_x = 1.0

    if not min_y == max_y:
        delta_y = max_y - min_y
    else:  # graph probably only has a single node
        delta_y = 1.0

    new_min_x, new_min_y, new_delta_x, new_delta_y = bbox

    new_coords_x = (coords[:, 0] - min_x) / delta_x * new_delta_x + new_min_x
    new_coords_y = (coords[:, 1] - min_y) / delta_y * new_delta_y + new_min_y

    return np.vstack([new_coords_x, new_coords_y]).T