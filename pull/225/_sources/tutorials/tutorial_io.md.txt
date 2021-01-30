---
jupyter:
  jupytext:
    notebook_metadata_filter: -kernelspec
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.5.0.rc1
---

```python
%load_ext autoreload
%autoreload 2
import sys

sys.path.insert(0, "../..")
import scirpy as ir
import scanpy as sc
from glob import glob
import pandas as pd
import tarfile
import anndata
import warnings
from numba import NumbaPerformanceWarning

# ignore numba performance warnings
warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)

# suppress "storing XXX as categorical" warnings.
anndata.logging.anndata_logger.setLevel("ERROR")
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
.. _importing-data:

Loading adaptive Immune Receptor (:term:`IR`)-sequencing data with Scirpy
=========================================================================

In this notebook, we demonstrate how single-cell :term:`IR`-data can be imported into
an :class:`~anndata.AnnData` object for the use with Scirpy. To learn more about
AnnData and how Scirpy makes use of it, check out the :ref:`data-structure` section.

The example data used in this notebook are available from the
`Scirpy repository <https://github.com/icbi-lab/scirpy/tree/master/docs/tutorials/example_data>`__.


.. important:: **The Scirpy data model**

    Currently, the Scirpy data model has the following constraints:

     * BCR and TCR chains are supported. Chain loci must be valid :term:`Chain locus`,
       i.e. one of `TRA`, `TRG`, `IGK`, or `IGL` (chains with a :term:`VJ<V(D)J>` junction) or
       `TRB`, `TRD`, or `IGH` (chains with a :term:`VDJ<V(D)J>` junction). Other chains are discarded.
     * Non-productive chains are removed. *CellRanger*, *TraCeR*, and the *AIRR rearrangment format*
       flag these cells appropriately. When reading :ref:`custom formats <importing-custom-formats>`,
       you need to pass the flag explicitly or filter the chains beforehand.
     * Each chain can contain up to two `VJ` and two `VDJ` chains (:term:`Dual IR`).
       Excess chains are removed (those with lowest read count/:term:`UMI` count)
       and cells flagged as :term:`Multichain-cell`.

    For more information, see :ref:`receptor-model`.


.. note:: **:term:`IR` quality control**

     * After importing the data, we recommend running the :func:`scirpy.tl.chain_qc` function.
       It will

           1. identify the :term:`Receptor type` and :term:`Receptor subtype` and flag cells
              as `ambiguous` that cannot unambigously be assigned to a certain receptor (sub)type, and
           2. flag cells with :term:`orphan chains <Orphan chain>` (i.e. cells with only a single detected cell)
              and :term:`multichain-cells <Multichain-cell>` (i.e. cells with more than two full pairs of VJ- and VDJ-chains).
     * We recommend excluding multichain- and ambiguous cells as these likely represent doublets
     * Based on the *orphan chain* flags, the corresponding cells can be excluded. Alternatively,
       these cells can be matched to clonotypes on a single chain only, by using the `receptor_arms="any"`
       parameter when running :func:`scirpy.tl.define_clonotypes`.


Loading data from *10x Genomics CellRanger*, *TraCeR*, *BraCer* or AIRR-compliant tools
---------------------------------------------------------------------------------------

We provide convenience functions to load data from *CellRanger*, *TraCeR*, or *BraCeR* with a single function call,
supporting both data generated on the *10x* and *Smart-seq2* sequencing platforms, respectively.
Moreover, we support importing data in the community-standard
`AIRR rearrangement schema <https://docs.airr-community.org/en/latest/datarep/rearrangements.html>`__.

.. module:: scirpy.io
   :noindex:

.. autosummary::
   :toctree: ../generated

   read_10x_vdj
   read_tracer
   read_bracer
   read_airr

Read 10x data
^^^^^^^^^^^^^

With :func:`~scirpy.io.read_10x_vdj` we can load `filtered_contig_annotations.csv` or `contig_annotations.json` files as they are produced by *CellRanger*.
Here, we demonstrate how to load paired single cell transcriptomics and TCR sequencing data from COVID19 patients
from `GSE145926 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926>`_ (:cite:`Liao2020`).
<!-- #endraw -->

```python
# Load the TCR data
adata_tcr = ir.io.read_10x_vdj(
    "example_data/liao-2019-covid19/GSM4385993_C144_filtered_contig_annotations.csv.gz"
)

# Load the associated transcriptomics data
adata = sc.read_10x_h5(
    "example_data/liao-2019-covid19/GSM4339772_C144_filtered_feature_bc_matrix.h5"
)
```

This particular sample only has a detected TCR for a small fraction of the cells:

```python
adata_tcr.shape
```

```python
adata.shape
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
Next, we integrate both the TCR and the transcriptomics data into a single :class:`anndata.AnnData` object
using :func:`scirpy.pp.merge_with_ir`:
<!-- #endraw -->

```python
ir.pp.merge_with_ir(adata, adata_tcr)
```

Now, we can use TCR-related variables together with the gene expression data.
Here, we visualize the cells with a detected TCR on the UMAP plot.
It is reassuring that the TCRs coincide with the T-cell marker gene CD3.

```python
sc.pp.log1p(adata)
sc.pp.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["has_ir", "CD3E"])
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
Read Smart-seq2 data processed with TraCeR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`TraCeR <https://github.com/Teichlab/tracer>`__ (:cite:`Stubbington2016-kh`) is a method commonly used
to extract TCR sequences from data generated with Smart-seq2 or other full-length single-cell sequencing protocols.
`Nf-core <https://nf-co.re/>`_ provides a full `pipeline for processing Smart-seq2 sequencing data <https://github.com/nf-core/smartseq2/>`__.

The :func:`scirpy.io.read_tracer` function obtains its TCR information from the `.pkl` file
in the `filtered_TCR_seqs` folder TraCeR generates for each cell.

For this example, we load the ~500 cells from triple-negative breast cancer patients from
`GSE75688 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688>`_ (:cite:`Chung2017`).
The raw data has been processed using the aforementioned `Smart-seq2 pipeline <https://github.com/nf-core/smartseq2/>`__ from nf-core.
<!-- #endraw -->

```python
# extract data
with tarfile.open("example_data/chung-park-2017.tar.bz2", "r:bz2") as tar:
    tar.extractall("example_data/chung-park-2017")
```

```python
# Load transcriptomics data from count matrix
expr_chung = pd.read_csv("example_data/chung-park-2017/counts.tsv", sep="\t")
# anndata needs genes in columns and samples in rows
expr_chung = expr_chung.set_index("Geneid").T
adata = sc.AnnData(expr_chung)
adata.shape
```

```python
# Load TCR data and merge it with transcriptomics data
adata_tcr = ir.io.read_tracer("example_data/chung-park-2017/tracer/")
ir.pp.merge_with_ir(adata, adata_tcr)
```

```python
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=3000)
sc.pp.log1p(adata)
sc.pp.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["has_ir", "CD3E"])
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
Read an AIRR-compliant rearrangement table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We generated example data using `immuneSIM <https://immunesim.readthedocs.io/en/latest/>`__ (:cite:`Weber2020`).
The data consists of 100 cells and does not include transcriptomics data.

The rearrangement tables are often organized into separate tables per chain. Therefore, :func:`scirpy.io.read_airr` supports
specifiying multiple `tsv` files at once. This would have the same effect as concatenating them before
the import.
<!-- #endraw -->

```python
adata = ir.io.read_airr(
    [
        "example_data/immunesim_airr/immunesim_tra.tsv",
        "example_data/immunesim_airr/immunesim_trb.tsv",
    ]
)
ir.tl.chain_qc(adata)
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
The dataset does not come with transcriptomics data. We can, therefore, not
show the UMAP plot highlighting cells with TCRs, but we can still use scirpy
to analyse it. Below, we visualize the clonotype network
connecting cells with similar :term:`CDR3` sequences.

**Note:** The cutoff of 25 was chosen for demonstration purposes on this small sample dataset. Usually a smaller cutoff
is more approriate.
<!-- #endraw -->

```python
ir.pp.ir_neighbors(
    adata,
    metric="alignment",
    sequence="aa",
    cutoff=25,
    receptor_arms="any",
    dual_ir="primary_only",
)
```

```python
ir.tl.define_clonotype_clusters(adata, metric="alignment", sequence="aa")
ir.tl.clonotype_network(adata, layout="fr", metric="alignment", sequence="aa")
ir.pl.clonotype_network(adata, color="ct_cluster_aa_alignment", panel_size=(4, 4))
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
.. _importing-custom-formats:

Creating AnnData objects from other formats
-------------------------------------------

Often, immune receptor (IR) data are just provided as a simple table listing the :term:`CDR3` sequences for each cell.
We provide a generic data structure for cells with IRs, which can then be converted into
an :class:`~anndata.AnnData` object.

.. module:: scirpy.io
   :noindex:

.. autosummary::
   :toctree: ../generated

   IrCell
   IrChain
   from_ir_objs

If you believe you are working with a commonly used format, consider sending a `feature request <https://github.com/icbi-lab/scirpy/issues>`_
for a `read_XXX` function.

For this example, we again load the triple-negative breast cancer data from :cite:`Chung2017`. However, this
time, we retrieve the TCR data from a separate summary table containing the TCR information
(we generated this table for the sake of the example, but it could as well
be a supplementary file from the paper).

Such a table typically contains information about

 * CDR3 sequences (amino acid and/or nucleotide)
 * The :term:`Chain locus`, e.g. `TRA`, `TRB`, or `IGH`.
 * expression of the receptor chain (e.g. count, :term:`UMI`, transcripts per million (TPM))
 * the :term:`V(D)J` genes for each chain
 * information if the chain is :term:`productive <Productive chain>`.
<!-- #endraw -->

```python
tcr_table = pd.read_csv(
    "example_data/chung-park-2017/tcr_table.tsv",
    sep="\t",
    index_col=0,
    na_values=["None"],
    true_values=["True"],
)
tcr_table
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
Our task is now to dissect the table into :class:`~scirpy.io.IrCell` and :class:`~scirpy.io.IrChain` objects.
Each :class:`~scirpy.io.IrCell` can have an arbitrary number of chains.
When converting the :class:`~scirpy.io.IrCell` objects into an :class:`~anndata.AnnData` object,
scirpy will only retain at most two alpha and two beta chains per cell and flag cells which exceed
this number as :term:`multichain cells <Multichain-cell>`. For more information, check the page about our :ref:`receptor-model`.
<!-- #endraw -->

```python
tcr_cells = []
for idx, row in tcr_table.iterrows():
    cell = ir.io.IrCell(cell_id=row["cell_id"])
    alpha_chain = ir.io.IrChain(
        locus="TRA",
        cdr3=row["cdr3_alpha"],
        cdr3_nt=row["cdr3_nt_alpha"],
        expr=row["count_alpha"],
        v_gene=row["v_alpha"],
        j_gene=row["j_alpha"],
        is_productive=row["productive_alpha"],
    )
    beta_chain = ir.io.IrChain(
        locus="TRB",
        cdr3=row["cdr3_beta"],
        cdr3_nt=row["cdr3_nt_beta"],
        expr=row["count_beta"],
        v_gene=row["v_beta"],
        d_gene=row["d_beta"],
        j_gene=row["j_beta"],
        is_productive=row["productive_beta"],
    )
    cell.add_chain(alpha_chain)
    cell.add_chain(beta_chain)
    tcr_cells.append(cell)
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
Now, we can convert the list of :class:`~scirpy.io.IrCell` objects using :func:`scirpy.io.from_ir_objs`.
<!-- #endraw -->

```python
adata_tcr = ir.io.from_ir_objs(tcr_cells)
```

```python
# We can re-use the transcriptomics data from above...
adata = sc.AnnData(expr_chung)
# ... and merge it with the TCR data
ir.pp.merge_with_ir(adata, adata_tcr)
```

```python
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=3000)
sc.pp.log1p(adata)
sc.pp.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["has_ir", "CD3E"])
```

<!-- #raw raw_mimetype="text/restructuredtext" -->
Combining multiple samples
--------------------------

It is quite common that the sequncing data is split up in multiple samples.
To combine them into a single object, we load each sample independently using one of the approaches described
in this document. Then, we combine them using :meth:`anndata.AnnData.concatenate`.

Here is a full example loading and combining three samples from the COVID19 study by :cite:`Liao2020`.
<!-- #endraw -->

```python
# define sample metadata. Usually read from a file.
samples = {
    "C144": {"group": "mild"},
    "C146": {"group": "severe"},
    "C149": {"group": "healthy control"},
}
```

```python
# Create a list of AnnData objects (one for each sample)
adatas = []
for sample, sample_meta in samples.items():
    gex_file = glob(f"example_data/liao-2019-covid19/*{sample}*.h5")[0]
    tcr_file = glob(f"example_data/liao-2019-covid19/*{sample}*.csv.gz")[0]
    adata = sc.read_10x_h5(gex_file)
    adata_tcr = ir.io.read_10x_vdj(tcr_file)
    ir.pp.merge_with_ir(adata, adata_tcr)
    adata.obs["sample"] = sample
    adata.obs["group"] = sample_meta["group"]
    # concatenation only works with unique gene names
    adata.var_names_make_unique()
    adatas.append(adata)
```

```python
# Merge anndata objects
adata = adatas[0].concatenate(adatas[1:])
```

The data is now integrated in a single object.
Again, the detected TCRs coincide with `CD3E` gene expression.
We clearly observe batch effects between the samples -- for a meaningful downstream analysis further
processing steps such as highly-variable gene filtering and batch correction are necessary.

```python
sc.pp.log1p(adata)
sc.pp.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["has_ir", "CD3E", "sample"])
```
