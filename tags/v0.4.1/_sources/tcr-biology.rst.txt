.. _tcr-model:

T-cell receptor model
=========================

What is a clonotype?
^^^^^^^^^^^^^^^^^^^^

A clonotype designates a collection of T or B cells that descend from a 
common, antecedent cell, and therefore, bear the same adaptive
immune receptors and recognize the same :term:`epitopes<Epitope>`.

In single-cell RNA-sequencing (scRNA-seq) data, T cells sharing identical
complementarity-determining regions 3 (:term:`CDR3`) nucleotide sequences of both α and β 
:term:`TCR` chains make up a clonotype. Scirpy provides an option to additionally
require clonotypes to have the same :term:`V-gene <V(D)J>`, enforcing the CDR 1 
and 2 regions to be the same. 


Dual TCRs
^^^^^^^^^

Contrary to what would be expected based on the previously described mechanism of
allelic exclusion (:cite:`Brady2010-gh`), scRNA-seq datasets can feature a considerable
number of cells with more than one TCR α and β pair. Since cells with more than one 
productive CDR3 sequence for each chain did not fit into our understanding of 
T cell biology, most TCR analysis tools ignore these cells (:cite:`Fischer2019`, 
:cite:`Zhang2018-ip`), or select the CDR3 sequence with the highest expression level
(:cite:`Afik2017-sg`). While in some cases these double-TCR cells might represent 
artifacts (e.g. doublets of a CD8+ and a CD4+ T cell engaged in an immunological 
synapse), there is an increasing amount of evidence in support of a bone fide
dual-TCR population (:cite:`Schuldt2019`, :cite:`Ji2010-bn`).

Scirpy allows investigating the composition and phenotypes of both single- and dual-TCR 
T cells by leveraging a T cell model similar to the one proposed in 
:cite:`Stubbington2016-kh`, where T cells are allowed to have a primary and a secondary 
pair of α- and β chains. For each cell, the primary pair consists of the α- and β-chain 
with the highest read count. Likewise, the secondary pair is the pair of α/β-chains with
the second highest expression level. Based on the assumption that each cell has only two
copies of the underlying chromosome set, if more than two variants of a chain are 
recovered for the same cell, the excess TCR chains are ignored by Scirpy and the 
corresponding cells flagged as :term:`Multichain-cell`. This filtering strategy leaves the choice 
of discarding or including multichain cells in downstream analyses.


Clonotype definition
^^^^^^^^^^^^^^^^^^^^

Scirpy implements a network-based clonotype definition that enables clustering cells
into clonotypes based on the following options:

 - identical :term:`CDR3` nucleotide sequences;
 - identical CDR3 amino acid sequences;
 - similar CDR3 amino acid sequences based on pairwise sequence alignment.

The latter approach is inspired by studies showing that similar TCR sequences also 
share epitope targets (:cite:`Fischer2019`, :cite:`Glanville2017-ay`, :cite:`TCRdist`).
Based on these approaches, Scirpy constructs a global, epitope-focused cell similarity
network. While convergence of the nucleotide-based clonotype definition to the amino 
acid based one hints at selection pressure, sequence alignment based networks offer
the opportunity to identify cells that might recognize the same epitopes.



