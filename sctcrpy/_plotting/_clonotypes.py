import scanpy as sc
from anndata import AnnData
import numpy as np
from typing import Union, Collection, Sequence
import matplotlib.pyplot as plt
from contextlib import contextmanager
import collections.abc as cabc
import warnings
from scanpy import settings
from matplotlib import cycler
from cycler import Cycler
import matplotlib


@contextmanager
def _patch_plot_edges(neighbors_key):
    """Monkey-patch scanpy's plot_edges to take our adjacency matrices"""
    scanpy_plot_edges = sc.plotting._utils.plot_edges

    def plot_edges(*args, **kwargs):
        return _plot_edges(*args, neighbors_key=neighbors_key, **kwargs)

    sc.plotting._utils.plot_edges = plot_edges
    try:
        yield
    finally:
        sc.plotting._utils.plot_edges = scanpy_plot_edges


@contextmanager
def _no_matplotlib_warnings():
    """Temporarily suppress matplotlib warnings"""
    import logging

    mpl_logger = logging.getLogger("matplotlib")
    log_level = mpl_logger.level
    mpl_logger.setLevel(logging.ERROR)
    try:
        yield
    finally:
        mpl_logger.setLevel(log_level)


def _plot_edges(axs, adata, basis, edges_width, edges_color, neighbors_key=None):
    import networkx as nx

    if not isinstance(axs, cabc.Sequence):
        axs = [axs]

    if neighbors_key is None:
        neighbors_key = "neighbors"
    if neighbors_key not in adata.uns:
        raise ValueError("`edges=True` requires `pp.neighbors` to be run before.")
    neighbors = adata.uns[neighbors_key]
    idx = np.where(~np.any(np.isnan(adata.obsm["X_" + basis]), axis=1))[0]
    g = nx.Graph(neighbors["connectivities"]).subgraph(idx)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for ax in axs:
            edge_collection = nx.draw_networkx_edges(
                g,
                adata.obsm["X_" + basis],
                ax=ax,
                width=edges_width,
                edge_color=edges_color,
            )
            edge_collection.set_zorder(-2)
            edge_collection.set_rasterized(settings._vector_friendly)


def clonotype_network(
    adata: AnnData,
    *,
    color: Union[str, Collection[str], None] = None,
    ax: Union[plt.Axes, Sequence[plt.Axes], None] = None,
    legend_loc: str = "right margin",
    palette: Union[str, Sequence[str], Cycler, None] = None,
    neighbors_key="tcr_neighbors",
    basis="clonotype_network",
    **kwargs
):
    """Plot the clonotype network"""
    # larger default size for figures when only one color is selected
    if (isinstance(color, str) or color is None) and ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    # use a cycler palette for many categories
    if isinstance(color, str) and adata.obs[color].unique().size > 40:
        if palette is None:
            palette = cycler(color=matplotlib.cm.Set3(range(12)))

    # for clonotype, use "on data" as default
    if isinstance(color, str) and color == "clonotype":
        if legend_loc == "right margin":
            legend_loc = "on data"

    with _patch_plot_edges(neighbors_key):
        with _no_matplotlib_warnings():
            return sc.pl.embedding(
                adata,
                basis="clonotype_network",
                color=color,
                ax=ax,
                edges=True,
                legend_loc=legend_loc,
                palette=palette,
                **kwargs,
            )
