
Glossary
========

.. glossary::
    TCR
        T-cell receptor. A TCR consists of one α and one β chain
        (or, alternatively, one γ and one δ chain). Each chain consists of a 
        constant and a variable region. The variable region is responsible for antigen
        recognition, mediated by :term:`CDR` regions. 

        *Scirpy* currently only supports α/β-TCRs. For more information, see the
        page about our :ref:`TCR model<tcr-model>`. 

        .. figure:: img/tcr.jpg
           :width: 310px

           Image from `Wikimedia commons <https://commons.wikimedia.org/wiki/File:2215_Alpha-Beta_T_Cell_Receptor.jpg>`_
           under the `CC BY-3.0 <https://creativecommons.org/licenses/by/3.0/deed.en>`__ license. 

    Clonotype
        A clonotype designates a collection of T or B cells that descend from a 
        common, antecedent cell, and therefore, bear the same adaptive
        immune receptors and recognize the same :term:`epitopes<Epitope>`.

        In single-cell RNA-sequencing (scRNA-seq) data, T cells sharing identical
        complementarity-determining regions 3 (:term:`CDR3`) nucleotide sequences of both α and β 
        :term:`TCR` chains make up a clonotype.

        *Scirpy* provides a flexible approach to clonotype definition based on 
        :term:`CDR3<CDR>` sequence identity or similarity. Additionally, it is possible
        to require clonotypes to have the same :term:`V-gene <V(D)J>`, enforcing the CDR 1 
        and 2 regions to be the same. 

        For more details, see the page about our :ref:`TCR model<tcr-model>` and
        the API documentation of :func:`scirpy.tl.define_clonotypes`. 

    Clonotype cluster
        A higher-order aggregation of :term:`clonotypes <Clonotype>` that have different 
        :term:`CDR3<CDR>` nucleotide sequences, but might recognize the same antigen 
        because they have the same or similar CDR3 amino acid sequence. 

        See also: :func:`scirpy.tl.define_clonotype_clusters`. 

    Private clonotype
        A clonotype that is specific for a certain patient.

    Public clonotype
        A clonotype that is shared across multiple 
        patients, e.g. a clonotype recognizing common viral epitope. 

        .. figure:: img/public-private.jpg
           :width: 420px

           Image from :cite:`Setliff2018` under the `CC BY-4.0 <https://creativecommons.org/licenses/by/4.0/>`__ license. 

    Tissue-specific clonotype
        A clonotype that only occurs in a certain tissue of a certain patient. 

    Multi-tissue clonotype
        A clonotype that occurs in multiple tissues of the same patient.

    Epitope
        The part of an antigen that is recognized by the :term:`TCR`
        (or B-cell receptor, or antibody). 

    CDR3
        Complementary-determining region 3. See :term:`CDR`. 

    CDR
        Complementary-determining region. The diversity and, therefore, antigen-specificity
        of :term:`TCRs<TCR>` is predominanly determined by three hypervariable loops 
        (CDR1, CDR2, and CDR3) on each of the α- and β receptor arms. 

        CDR1 and CDR2 are fully encoded in germline V genes. In contrast, 
        the CDR3 loops are assembled from V and J segments (TCR-α) and V, D and 
        J segments (TCR-β) and comprise random additions and deletions at the junction
        sites (see also :term:`V(D)J`). Thus, CDR3 regions make up a large part of the 
        TCR variability and are therefore thought to be particularly important for antigen 
        specificity (reviewed in :cite:`Attaf2015`). 

        .. figure:: img/tcr_cdr3.png
           :width: 310px
           
           Image from :cite:`Attaf2015` under the `CC BY-NC-SA-3.0 <https://creativecommons.org/licenses/by-nc-sa/3.0/>`__ license.

    V(D)J
        The variability of :term:`TCR` chain sequences originates from the genetic recombination
        of **V**\ ariable, **D**\ iversity and **J**\ oining gene segments. The TCR-α
        chain gets assembled from V and J loci only, the TCR-β chain from all three
        V, D and J loci. 

        As an example, the figure below shows how a TCR-α chain is assembed from 
        the *tra* locus. V to J recombination joins one of many TRAV segments to one of 
        many TRAJ segments. Next, introns are spliced out, resulting in a TCR-α chain 
        transcript with V, J and C segments directly next to each other (reviewed in :cite:`Attaf2015`).  

        .. figure:: img/vdj.png
           :width: 600px
           
           Image from :cite:`Attaf2015` under the `CC BY-NC-SA-3.0 <https://creativecommons.org/licenses/by-nc-sa/3.0/>`__ license.
        
    Dual TCR 
        :term:`TCRs<TCR>` with more than one pair of α- and β chains. While this was
        previously thought to be impossible due to the mechanism of allelic exclusion
        (:cite:`Brady2010-gh`), there is an increasing amound of evidence for a *bona fide*
        dual-TCR population (:cite:`Schuldt2019`, :cite:`Ji2010-bn`). 

        For more information on how *Scirpy* handles dual TCRs, see the
        page about our :ref:`TCR model<tcr-model>`.  
     

    Multichain-cell
        Cells with more than two α- and β chains that do not fit into the
        :term:`Dual TCR` model. These are usually rare and could be explained 
        by doublets/multiplets, i.e. two ore more cells that were captured 
        in the same droplet.

        .. figure:: img/multichain.png
           :width: 450px 

           (a) UMAP plot of 96,000 cells from :cite:`Wu2020` with at least one detected
           :term:`CDR3` sequence with multichain-cells (n=474) highlighted in green. 
           (b)  Comparison of detected reads per cell in multichain-cells and other cells.
           Multichain cells comprised significantly more reads per cell
           (p = 9.45 × 10−251, Wilcoxon-Mann-Whitney-test), supporting the hypothesis that
           (most of) multichain cells are technical artifacts arising from 
           cell-multiplets (:cite:`Ilicic2016`).

    Orphan chain
        A :term:`TCR` chain is called *orphan*, if its corresponding counterpart 
        has not been detected. For instance, if a cell has only a TCR-α chain,
        but no TCR-β chain, the cell will be flagged as "Orphan alpha". 

        Orphan chains are most likley the effect of stochastic dropouts due to 
        sequencing inefficiencies. 

        See also :func:`scirpy.tl.chain_pairing`. 

    UMI
        Unique molecular identifier. Some single-cell RNA-seq protocols
        label each RNA with a unique barcode prior to PCR-amplification to mitigate
        PCR bias. With these protocols, UMI-counts replace the read-counts 
        generally used with RNA-seq.

    productive chain
        Productive chains are TCR chains with a :term:`CDR3` sequence that produces
        a functional peptide. Scirpy relies on the preprocessing tools (e.g. 
        CellRanger or TraCeR) for flagging non-productive chains. 
        Typically chains are flagged as non-productive if they contain
        a stop codon or are not within the reading frame. 
