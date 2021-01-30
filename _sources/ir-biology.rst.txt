.. _receptor-model:

Immune receptor (IR) model
==========================

What is a clonotype?
^^^^^^^^^^^^^^^^^^^^

A clonotype designates a collection of T or B cells that descend from a
common, antecedent cell, and therefore, bear the same adaptive
immune receptors (:term:`IR`) and recognize the same :term:`epitopes<Epitope>`.

In single-cell RNA-sequencing (scRNA-seq) data, T or B cells sharing identical
complementarity-determining regions 3 (:term:`CDR3`) nucleotide sequences of both
:term:`VJ<V(D)J>` and :term:`VDJ<V(D)J>` chains (e.g. both α and β :term:`TCR` chains)
make up a clonotype. Scirpy provides an option to additionally
require clonotypes to have the same :term:`V-gene <V(D)J>`, enforcing the CDR 1
and 2 regions to be the same.


Dual TCRs and allelically included B cells
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Contrary to what would be expected based on the previously described
mechanism of allelic exclusion (:cite:`Brady2010-gh`), scRNA-seq datasets can
feature a considerable number of cells with more than one pair of
:term:`VJ<V(D)J>` and :term:`VDJ<V(D)J>` chains. Since cells with more than
one productive CDR3 sequence for each chain did not fit into our
understanding of adaptive immune cell biology, most TCR analysis tools ignore
these cells (:cite:`Fischer2019`, :cite:`Zhang2018-ip`), or select the CDR3
sequence with the highest expression level (:cite:`Afik2017-sg`). While in
some cases these cells might represent artifacts (e.g. doublets of a CD8+ and
a CD4+ T cell engaged in an immunological synapse), there is an increasing
amount of evidence in support of a bone fide :term:`Dual IR` population
(:cite:`Schuldt2019`, :cite:`Ji2010-bn`, :cite:`Vettermann2010`).

Scirpy allows investigating the composition and phenotypes of both single-
and dual-IR cells by leveraging a immune cell model similar to the one
proposed in :cite:`Stubbington2016-kh`, where immune cells are allowed to
have a primary and a secondary pair of :term:`VJ<V(D)J>` and
:term:`VDJ<V(D)J>` chains. For each cell, the primary pair consists of the VJ-
and VDJ-chain with the highest read count. Likewise, the secondary pair is the
pair of VJ/VDJ-chains with the second highest expression level. Based on the
assumption that each cell has only two copies of the underlying chromosome
set, if more than two variants of a chain are recovered for the same cell,
the excess IR chains are ignored by Scirpy and the corresponding cells
flagged as :term:`Multichain-cell`. This filtering strategy leaves the choice
of discarding or including multichain cells in downstream analyses.


Clonotype definition
^^^^^^^^^^^^^^^^^^^^

Scirpy implements a network-based clonotype definition that enables
clustering cells into :term:`clonotypes<Clonotype>` based on *identical*
:term:`CDR3` *nucleotide sequences* and into :term:`clonotype clusters<Clonotype cluster>`
based on

- *identical CDR3 amino acid sequences*, or
- *similar CDR3 amino acid sequences* based on pairwise sequence alignment.

The latter approach is inspired by studies showing that similar TCR sequences
also share epitope targets (:cite:`Fischer2019`, :cite:`Glanville2017-ay`,
:cite:`TCRdist`). Based on these approaches, Scirpy constructs a global,
epitope-focused cell similarity network. While convergence of the
nucleotide-based clonotype definition to the amino acid based one hints at
selection pressure, sequence alignment based networks offer the opportunity
to identify cells that might recognize the same epitopes.
