---
jupyter:
  jupytext:
    formats: md,ipynb
    notebook_metadata_filter: -kernelspec
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.4.2
---

# Analysis of 3k T cells from cancer

<!-- #raw raw_mimetype="text/restructuredtext" -->
In this tutorial, we re-analyze single-cell TCR/RNA-seq data from Wu et al. (:cite:`Wu2020`)
generated on the 10x Genomics platform. The original dataset consists of >140k T cells
from 14 treatment-naive patients across four different types of cancer.
For this tutorial, to speed up computations, we use a downsampled version of 3k cells.
<!-- #endraw -->

```python
%load_ext autoreload
%autoreload 2
import sys
import warnings

import numpy as np
import pandas as pd

import scanpy as sc
import scirpy as ir
from matplotlib import pyplot as plt

sys.path.insert(0, "../..")

warnings.filterwarnings("ignore", category=FutureWarning)
```

```python
sc.logging.print_versions()
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
The dataset ships with the `scirpy` package. We can conveniently load it from the :mod:`~scirpy.datasets` module:
<!-- #endraw -->

```python
adata = ir.datasets.wu2020_3k()
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
`adata` is a regular :class:`~anndata.AnnData` object with additional, TCR-specific columns in `obs`. 
For more information, check the page about Scirpy's :ref:`data structure <data-structure>`. 

.. note:: For more information about our T-cell receptor model, see :ref:`tcr-model`. 
<!-- #endraw -->

```python
adata.shape
```

```python
adata.obs
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
.. note:: **Importing data**

    `scirpy` natively supports reading TCR data from `Cellranger <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger>`_ (10x), `TraCeR <https://github.com/Teichlab/tracer>`_ (Smart-seq2) 
    or the `AIRR rearrangement schema <https://docs.airr-community.org/en/latest/datarep/rearrangements.html>`__ and provides helper functions to import other data types. We provide a :ref:`dedicated tutorial on data loading <importing-data>` with more details. 

    This particular dataset has been imported using :func:`scirpy.io.read_10x_vdj` and merged
    with transcriptomics data using :func:`scirpy.pp.merge_with_tcr`. The exact procedure
    is described in :func:`scirpy.datasets.wu2020`.

<!-- #endraw -->

## Preprocess Transcriptomics data

Transcriptomics data needs to be filtered and preprocessed as with any other single-cell dataset.
We recommend following the [scanpy tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html)
and the best practice paper by [Luecken et al.](https://www.embopress.org/doi/10.15252/msb.20188746)

```python
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.filter_cells(adata, min_genes=100)
```

```python
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1000)
sc.pp.log1p(adata)
```

For the _Wu2020_ dataset, the authors already provide clusters and UMAP coordinates.
Instead of performing clustering and cluster annotation ourselves, we will just use
the provided data. The clustering and annotation methodology is 
described in [their paper](https://doi.org/10.1038/s41586-020-2056-8). 

```python
adata.obsm["X_umap"] = adata.obsm["X_umap_orig"]
```

```python
mapping = {
    "3.1-MT": "other",
    "4.1-Trm": "CD4_Trm",
    "4.2-RPL32": "CD4_RPL32",
    "4.3-TCF7": "CD4_TCF7",
    "4.4-FOS": "CD4_FOSS",
    "4.5-IL6ST": "CD4_IL6ST",
    "4.6a-Treg": "CD4_Treg",
    "4.6b-Treg": "CD4_Treg",
    "8.1-Teff": "CD8_Teff",
    "8.2-Tem": "CD8_Tem",
    "8.3a-Trm": "CD8_Trm",
    "8.3b-Trm": "CD8_Trm",
    "8.3c-Trm": "CD8_Trm",
    "8.4-Chrom": "other",
    "8.5-Mitosis": "other",
    "8.6-KLRB1": "other",
}
adata.obs["cluster"] = [mapping[x] for x in adata.obs["cluster_orig"]]
```

Let's inspect the UMAP plots. The first three panels show the UMAP plot colored by sample, patient and cluster.
We don't observe any clustering of samples or patients that could hint at batch effects.

The last three panels show the UMAP colored by the T cell markers _CD8_, _CD4_, and _FOXP3_.
We can confirm that the markers correspond to their respective cluster labels.

```python
sc.pl.umap(
    adata,
    color=["sample", "patient", "cluster", "CD8A", "CD4", "FOXP3"],
    ncols=2,
    wspace=0.5,
)
```

## TCR Quality Control

<!-- #raw raw_mimetype="text/restructuredtext" -->
While most of T cell receptors have exactly one pair of α and β chains, up to one third of 
T cells can have *dual TCRs*, i.e. two pairs of receptors originating from different alleles (:cite:`Schuldt2019`).

Using the :func:`scirpy.tl.chain_pairing` function, we can add a summary
about the T cell receptor compositions to `adata.obs`. We can visualize it using :func:`scirpy.pl.group_abundance`.

.. note:: **chain pairing**

    - *Orphan chain* refers to cells that have either a single alpha or beta receptor chain.
    - *Extra chain* refers to cells that have a full alpha/beta receptor pair, and an additional chain.
    - *Multichain* refers to cells with more than two receptor pairs detected. These cells are likely doublets.
<!-- #endraw -->

```python
ir.tl.chain_pairing(adata)
```

```python
ir.pl.group_abundance(
    adata, groupby="chain_pairing", target_col="source",
)
```

Indeed, in this dataset, ~6% of cells have more than 
one pair of productive T-cell receptors:

```python
print(
    "Fraction of cells with more than one pair of TCRs: {:.2f}".format(
        np.sum(
            adata.obs["chain_pairing"].isin(
                ["Extra beta", "Extra alpha", "Two full chains"]
            )
        )
        / adata.n_obs
    )
)
```

Next, we visualize the _Multichain_ cells on the UMAP plot and exclude them from downstream analysis:

```python
sc.pl.umap(adata, color="multi_chain")
```

```python
adata = adata[adata.obs["multi_chain"] != "True", :].copy()
```

## Define clonotypes and clonotype clusters

<!-- #raw raw_mimetype="text/restructuredtext" -->

In this section, we will define and visualize :term:`clonotypes <Clonotype>` and :term:`clonotype clusters <Clonotype cluster>`.

*Scirpy* implements a network-based approach for clonotype definition. The steps to create and visualize the clonotype-network are analogous to the construction of a neighborhood graph from transcriptomics data with *scanpy*.

.. list-table:: Analysis steps on transcriptomics data
    :widths: 40 60
    :header-rows: 1

    * - scanpy function
      - objective
    * - :func:`scanpy.pp.neighbors`
      - Compute a nearest-neighbor graph based on gene expression.
    * - :func:`scanpy.tl.leiden`
      - Cluster cells by the similarity of their transcriptional profiles.
    * - :func:`scanpy.tl.umap`
      - Compute positions of cells in UMAP embedding.
    * - :func:`scanpy.pl.umap`
      - Plot UMAP colored by different parameters.

.. list-table:: Analysis steps on TCR data
    :widths: 40 60
    :header-rows: 1

    - - scirpy function
      - objective
    - - :func:`scirpy.pp.tcr_neighbors`
      - Compute a neighborhood graph of CDR3-sequences.
    - - :func:`scirpy.tl.define_clonotypes`
      - Define :term:`clonotypes <Clonotype>` by nucleotide 
        sequence identity.
    - - :func:`scirpy.tl.define_clonotype_clusters`
      - Cluster cells by the similarity of their CDR3-sequences
    - - :func:`scirpy.tl.clonotype_network`
      - Compute positions of cells in clonotype network.
    - - :func:`scirpy.pl.clonotype_network`
      - Plot clonotype network colored by different parameters.

<!-- #endraw -->

### Compute CDR3 neighborhood graph and define clonotypes

<!-- #raw raw_mimetype="text/restructuredtext" -->
:func:`scirpy.pp.tcr_neighbors` computes a neighborhood graph based on :term:`CDR3 <CDR>` nucleotide (`nt`) or amino acid (`aa`) sequences, either based on sequence identity or similarity. 

Here, we define :term:`clonotypes <Clonotype>` based on nt-sequence identity. 
In a later step, we will define :term:`clonotype clusters <Clonotype cluster>` based on 
amino-acid similarity. 

The function :func:`scirpy.tl.define_clonotypes` will detect connected modules
in the graph and annotate them as clonotypes. This will add a `clonotype` and
`clonotype_size` column to `adata.obs`.
<!-- #endraw -->

```python
# using default parameters, `tcr_neighbors` will compute nucleotide sequence identity
ir.pp.tcr_neighbors(adata, receptor_arms="all", dual_tcr="primary_only")
ir.tl.define_clonotypes(adata)
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
To visualize the network we first call :func:`scirpy.tl.clonotype_network` to compute the layout.
We can then visualize it using :func:`scirpy.pl.clonotype_network`. We recommend setting the
`min_size` parameter to `>=2`, to prevent the singleton clonotypes from cluttering the network.
<!-- #endraw -->

```python
ir.tl.clonotype_network(adata, min_size=2)
ir.pl.clonotype_network(adata, color="clonotype", legend_loc="none")
```

### Re-compute CDR3 neighborhood graph and define clonotype clusters

<!-- #raw raw_mimetype="text/restructuredtext" -->
We can now re-compute the neighborhood graph based on amino-acid sequence similarity 
and define :term:`clonotype clusters <Clonotype cluster>`.

To this end, we need to change set `metric="alignment"` and specify a `cutoff` parameter. 
The distance is based on the `BLOSUM62 <https://en.wikipedia.org/wiki/BLOSUM>`__ matrix. 
For instance, a distance of `10` is equivalent to 2 `R`s mutating into `N`. 
This appoach was initially proposed as *TCRdist* by Dash et al. (:cite:`TCRdist`).

All cells with a distance between their CDR3 sequences lower than `cutoff` will be connected in the network. 
<!-- #endraw -->

```python
sc.settings.verbosity = 4
```

```python
ir.pp.tcr_neighbors(
    adata,
    metric="alignment",
    sequence="aa",
    cutoff=15,
    receptor_arms="all",
    dual_tcr="all",
)
ir.tl.define_clonotype_clusters(
    adata, partitions="connected", sequence="aa", metric="alignment"
)
```

```python
ir.tl.clonotype_network(adata, min_size=4, sequence="aa", metric="alignment")
```

Compared to the previous plot, we observe slightly larger clusters that are not necessarily fully connected any more. 

```python
ir.pl.clonotype_network(
    adata,
    color="ct_cluster_aa_alignment",
    legend_fontoutline=3,
    size=80,
    panel_size=(6, 6),
    legend_loc="on data",
)
```

Now we show the same graph, colored by patient.
We observe that for instance clonotypes 104 and 160 (center-left) are _private_, i.e. they contain cells from
a single sample only. On the other hand, for instance clonotype 233 (top-left) is _public_, i.e.
it is shared across patients _Lung5_ and _Lung1_ and _Lung3_. 

```python
ir.pl.clonotype_network(adata, color="patient", size=80, panel_size=(6, 6))
```

We can now extract information (e.g. CDR3-sequences) from a specific clonotype cluster by subsetting `AnnData`. 
For instance, we can find out that clonotype `233` does not have a detected alpha chain. 

```python
adata.obs.loc[
    adata.obs["ct_cluster_aa_alignment"] == "233",
    ["TRA_1_cdr3", "TRA_2_cdr3", "TRB_1_cdr3", "TRB_2_cdr3"],
]
```

### Including the V-gene in clonotype definition

<!-- #raw raw_mimetype="text/restructuredtext" -->
Using the paramter `use_v_gene` in :func:`~scirpy.tl.define_clonotypes`, we can enforce
clonotypes (or clonotype clusters) to have the same :term:`V-gene <V(D)J>`, and, therefore, the same :term:`CDR1 and 2 <CDR>`
regions. Let's look for clonotype clusters with different V genes:
<!-- #endraw -->

```python
ir.tl.define_clonotype_clusters(
    adata,
    sequence="aa",
    metric="alignment",
    same_v_gene="primary_only",
    key_added="ct_cluster_aa_alignment_same_v",
)
```

```python
# find clonotypes with more than one `clonotype_same_v`
ct_different_v = adata.obs.groupby("ct_cluster_aa_alignment").apply(
    lambda x: x["ct_cluster_aa_alignment_same_v"].unique().size > 1
)
ct_different_v = ct_different_v[ct_different_v].index.values
ct_different_v
```

```python
# Display the first 2 clonotypes with different v genes
adata.obs.loc[
    adata.obs["ct_cluster_aa_alignment"].isin(ct_different_v[:2]),
    [
        "ct_cluster_aa_alignment",
        "ct_cluster_aa_alignment_same_v",
        "TRA_1_v_gene",
        "TRB_1_v_gene",
    ],
].sort_values("ct_cluster_aa_alignment").drop_duplicates().reset_index(drop=True)
```

## Clonotype analysis

### Clonal expansion

<!-- #raw raw_mimetype="text/restructuredtext" -->
Let's visualize the number of expanded clonotypes (i.e. clonotypes consisting
of more than one cell) by cell-type. The first option is to add a column with the :func:`scirpy.tl.clonal_expansion` 
to `adata.obs` and overlay it on the UMAP plot. 
<!-- #endraw -->

```python
ir.tl.clonal_expansion(adata)
```

`clonal_expansion` refers to expansion categories, i.e singleton clonotypes, clonotypes with 2 cells and more than 2 cells. 
The `clonotype_size` refers to the absolute number of cells in a clonotype. 

```python
sc.pl.umap(adata, color=["clonal_expansion", "clonotype_size"])
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
The second option is to show the number of cells belonging to an expanded clonotype per category
in a stacked bar plot, using the :func:`scirpy.pl.clonal_expansion` plotting function. 
<!-- #endraw -->

```python
ir.pl.clonal_expansion(adata, groupby="cluster", clip_at=4, normalize=False)
```

The same plot, normalized to cluster size. Clonal expansion is a sign of positive selection
for a certain, reactive T-cell clone. It, therefore, makes sense that CD8+ effector T-cells 
have the largest fraction of expanded clonotypes. 

```python
ir.pl.clonal_expansion(adata, "cluster")
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
Expectedly, the CD8+ effector T cells have the largest fraction of expanded clonotypes. 

Consistent with this observation, they have the lowest :func:`scirpy.pl.alpha_diversity` of clonotypes.  
<!-- #endraw -->

```python
ax = ir.pl.alpha_diversity(adata, groupby="cluster")
```

### Clonotype abundance

<!-- #raw raw_mimetype="text/restructuredtext" -->
The function :func:`scirpy.pl.group_abundance` allows us to create bar charts for
arbitrary categorial from `obs`. Here, we use it to show the distribution of the 
ten largest clonotypes across the cell-type clusters.
<!-- #endraw -->

```python
ir.pl.group_abundance(adata, groupby="clonotype", target_col="cluster", max_cols=10)
```

It might be beneficial to normalize the counts
to the number of cells per sample to mitigate biases due to different sample sizes: 

```python
ir.pl.group_abundance(
    adata, groupby="clonotype", target_col="cluster", max_cols=10, normalize="sample"
)
```

Coloring the bars by patient gives us information about public and private clonotypes: 
Some clonotypes are *private*, i.e. specific to a certain tissue, 
others are *public*, i.e. they are shared across different tissues. 

```python
ax = ir.pl.group_abundance(
    adata, groupby="clonotype", target_col="source", max_cols=15, fig_kws={"dpi": 100}
)
```

However, clonotypes that are shared between *patients* are rare: 

```python
ax = ir.pl.group_abundance(
    adata, groupby="clonotype", target_col="patient", max_cols=15, fig_kws={"dpi": 100}
)
```

## Gene usage

<!-- #raw raw_mimetype="text/restructuredtext" -->
:func:`scirpy.tl.group_abundance` can also give us some information on VDJ usage. 
We can choose any of the `{TRA,TRB}_{1,2}_{v,d,j,c}_gene` columns to make a stacked bar plot. 
We use `max_col` to limit the plot to the 10 most abundant V-genes. 
<!-- #endraw -->

```python
ir.pl.group_abundance(
    adata, groupby="TRB_1_v_gene", target_col="cluster", normalize=True, max_cols=10
)
```

We can pre-select groups by filtering `adata`:

```python
ir.pl.group_abundance(
    adata[
        adata.obs["TRB_1_v_gene"].isin(
            ["TRBV20-1", "TRBV7-2", "TRBV28", "TRBV5-1", "TRBV7-9"]
        ),
        :,
    ],
    groupby="cluster",
    target_col="TRB_1_v_gene",
    normalize=True,
)
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
The exact combinations of VDJ genes can be visualized as a Sankey-plot using :func:`scirpy.pl.vdj_usage`. 
<!-- #endraw -->

```python
ir.pl.vdj_usage(adata, full_combination=False, top_n=30)
```

We can also use this plot to investigate the exact VDJ composition of one (or several) clonotypes: 

```python
ir.pl.vdj_usage(
    adata[adata.obs["clonotype"].isin(["274", "277", "211", "106"]), :], top_n=None
)
```

### Spectratype plots

<!-- #raw raw_mimetype="text/restructuredtext" -->
:func:`~scirpy.pl.spectratype` plots give us information about the length distribution of CDR3 regions. 
<!-- #endraw -->

```python
ir.pl.spectratype(adata, color="cluster", viztype="bar", fig_kws={"dpi": 120})
```

The same chart visualized as "ridge"-plot: 

```python
ir.pl.spectratype(
    adata,
    color="cluster",
    viztype="curve",
    curve_layout="shifted",
    fig_kws={"dpi": 120},
    kde_kws={"kde_norm": False},
)
```

A spectratype-plot by gene usage. To pre-select specific genes, we can simply filter the `adata` object before plotting. 

```python
ir.pl.spectratype(
    adata[
        adata.obs["TRB_1_v_gene"].isin(
            ["TRBV20-1", "TRBV7-2", "TRBV28", "TRBV5-1", "TRBV7-9"]
        ),
        :,
    ],
    cdr3_col="TRB_1_cdr3",
    color="TRB_1_v_gene",
    normalize="sample",
    fig_kws={"dpi": 120},
)
```

## Comparing repertoires

### Repertoire simlarity and overlaps

<!-- #raw raw_mimetype="text/restructuredtext" -->
Overlaps in the adaptive immune receptor repertoire of samples or sample groups enables to pinpoint important clonotype groups, as well as to provide a measure of similarity between samples.  
Running Scirpy's :func:`~scirpy.tl.repertoire_overlap` tool creates a matrix featuring the abundance of clonotypes in each sample. Additionally, it also computes a (Jaccard) distance matrix of samples as well as the linkage of hierarchical clustering. 
<!-- #endraw -->

```python
df, dst, lk = ir.tl.repertoire_overlap(adata, "sample", inplace=False)
df.head()
```

The distance matrix can be shown as a heatmap, while samples are reordered based on hierarchical clustering.

```python
ir.pl.repertoire_overlap(adata, "sample", heatmap_cats=["patient", "source"])
```

A specific pair of samples can be compared on a scatterplot, where dot size corresponds to the number of clonotypes at a given coordinate.

```python
ir.pl.repertoire_overlap(
    adata, "sample", pair_to_plot=["LN2", "LT2"], fig_kws={"dpi": 120}
)
```

### Clonotypes preferentially occuring in a group

<!-- #raw raw_mimetype="text/restructuredtext" -->
Clonotypes associated with an experimental group (a given cell type, samle or diagnosis) might be important candidates as biomarkers or disease drivers. Scirpy offers :func:`~scirpy.tl.clonotype_imbalance` to rank clonotypes based on Fisher's exact test comparing the fractional presence of a given clonotype in two groups.
<!-- #endraw -->

A possible grouping criterion could be Tumor vs. Control, separately for distinct tumor types. The site of the tumor can be extracted from patient metadata.

```python
adata.obs["site"] = adata.obs["patient"].str.slice(stop=-1)
```

```python
ir.pl.clonotype_imbalance(
    adata,
    replicate_col="sample",
    groupby="source",
    case_label="Tumor",
    additional_hue="site",
    plot_type="strip",
)
```

To get an idea how the above, top-ranked clonotypes compare to the bulk of all clonotypes, a Volcano plot is genereated, showing the `-log10 p-value` of the Fisher's test as a function of `log2(fold-change)` of the normalized proportion of a given clonotype in the test group compared to the control group. To avoid zero division, `0.01*(global minimum proportion)` was added to every normalized clonotype proportions.

```python
ir.pl.clonotype_imbalance(
    adata,
    replicate_col="sample",
    groupby="source",
    case_label="Tumor",
    additional_hue="diagnosis",
    plot_type="volcano",
    fig_kws={"dpi": 120},
)
```

## Integrating gene expression
### Clonotype imbalance among cell clusters

Leveraging the opportunity offered by close integeration with scanpy, transcriptomics-based data can be utilized directly. Using cell type annotation inferred from gene expression clusters, for example, clonotypes belonging to CD8+ effector T-cells and CD8+ tissue-resident memory T cells, can be compared.

```python
freq, stat = ir.tl.clonotype_imbalance(
    adata,
    replicate_col="sample",
    groupby="cluster",
    case_label="CD8_Teff",
    control_label="CD8_Trm",
    inplace=False,
)
top_differential_clonotypes = stat["clonotype"].tolist()[:5]
```

Showing top clonotypes on a UMAP clearly shows that clonotype 163 is featured by CD8+ tissue-resident memory T cells, while clonotype 277 by CD8+ effector T-cells.

```python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw={"wspace": 0.6})
sc.pl.umap(adata, color="cluster", ax=ax1, show=False)
sc.pl.umap(
    adata,
    color="clonotype",
    groups=top_differential_clonotypes,
    ax=ax2,
    # increase size of highlighted dots
    size=[
        80 if c in top_differential_clonotypes else 30 for c in adata.obs["clonotype"]
    ],
)
```

### Repertoire overlap of cell types

Just like comparing repertoire overlap among samples, Scirpy also offers comparison between gene expression clusters or cell subpopulations. As an example, repertoire overlap of the two cell types compared above is shown.

```python
ir.tl.repertoire_overlap(adata, "cluster")
ir.pl.repertoire_overlap(
    adata, "cluster", pair_to_plot=["CD8_Teff", "CD8_Trm"], fig_kws={"dpi": 120}
)
```

### Marker genes in top clonotypes

Gene expression of cells belonging to individual clonotypes can also be compared using Scanpy. As an example, differential gene expression of two clonotypes, found to be specific to cell type clusters can also be analysed. 

```python
sc.tl.rank_genes_groups(
    adata, "clonotype", groups=["163"], reference="277", method="wilcoxon"
)
sc.pl.rank_genes_groups_violin(adata, groups="163", n_genes=15)
```
