
Glossary
========

.. glossary::
    :sorted:

    TCR
        T-cell receptor. A TCR consists of one α and one β chain
        (or, alternatively, one γ and one δ chain). Each chain consists of a
        constant and a variable region. The variable region is responsible for antigen
        recognition, mediated by :term:`CDR` regions.

        For more information on how Scirpy represents TCRs, see the page about our
        :ref:`receptor model<receptor-model>`.

        .. figure:: img/tcr.jpg
           :width: 310px

           Image from `Wikimedia commons <https://commons.wikimedia.org/wiki/File:2215_Alpha-Beta_T_Cell_Receptor.jpg>`_
           under the `CC BY-3.0 <https://creativecommons.org/licenses/by/3.0/deed.en>`__ license.

    Clonotype
        A clonotype designates a collection of T or B cells that descend from a
        common, antecedent cell, and therefore, bear the same adaptive
        immune receptors and recognize the same :term:`epitopes<Epitope>`.

        In single-cell RNA-sequencing (scRNA-seq) data, T or B cells sharing identical
        complementarity-determining regions 3 (:term:`CDR3`) nucleotide sequences of both
        :term:`VJ<V(D)J>` and :term:`VDJ<V(D)J>` chains (e.g. both α and β :term:`TCR` chains)
        make up a clonotype..

        *Scirpy* provides a flexible approach to clonotype definition based on
        :term:`CDR3<CDR>` sequence identity or similarity. Additionally, it is possible
        to require clonotypes to have the same :term:`V-gene <V(D)J>`, enforcing the CDR 1
        and 2 regions to be the same.

        For more details, see the page about our :ref:`IR model<receptor-model>` and
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

    Convergent evolution of clonotypes
        It has been proposed that :term:`IRs <IR>` are subject to convergent evolution,
        i.e. a selection pressure that leads to IRs recognizing the same antigen
        (:cite:`Venturi2006`).

        Evidence of convergent evolution could be clonotypes with the same :term:`CDR3`
        amino acid sequence, but different CDR3 nucleotide sequences (due to synonymous
        codons) or clonotypes with highly similar CDR3 amino acid sequences that
        recognize the same antigen.


    Epitope
        The part of an antigen that is recognized by the :term:`TCR`, :term:`BCR`,
        or antibody.

    CDR3
        Complementary-determining region 3. See :term:`CDR`.

    CDR
        Complementary-determining region. The diversity and, therefore, antigen-specificity
        of :term:`IRs<IR>` is predominanly determined by three hypervariable loops
        (CDR1, CDR2, and CDR3) on each of the α- and β receptor arms.

        CDR1 and CDR2 are fully encoded in germline V genes. In contrast,
        the CDR3 loops are assembled from :term:`V, (D), and J segments<V(D)J>` and
        comprise random additions and deletions at the junction
        sites. Thus, CDR3 regions make up a large part of the
        adpative immune receptor variability and are therefore thought to
        be particularly important for antigen specificity
        (reviewed in :cite:`Attaf2015`).

        .. figure:: img/tcr_cdr3.png
           :width: 310px

           Image from :cite:`Attaf2015` under the `CC BY-NC-SA-3.0 <https://creativecommons.org/licenses/by-nc-sa/3.0/>`__ license.

    V(D)J
        The variability of :term:`IR` chain sequences originates from the genetic recombination
        of **V**\ ariable, **D**\ iversity and **J**\ oining gene segments. The :term:`TCR`-α,
        TCR-ɣ, :term:`IG`-κ, and IG-λ chains get assembled from V and J loci only. We refer
        to these chains as `VJ` chains in Scirpy. The TCR-β, TCR-δ, and IG-heavy chains
        get assembled from all three segments. We refer to these chains as `VDJ`-chains
        in Scirpy.

        As an example, the figure below shows how a TCR-α chain is assembed from
        the *tra* locus. V to J recombination joins one of many `TRAV` segments to one of
        many `TRAJ` segments. Next, introns are spliced out, resulting in a TCR-α chain
        transcript with V, J and C segments directly next to each other (reviewed in :cite:`Attaf2015`).

        .. figure:: img/vdj.png
           :width: 600px

           Image from :cite:`Attaf2015` under the `CC BY-NC-SA-3.0 <https://creativecommons.org/licenses/by-nc-sa/3.0/>`__ license.

    Dual TCR
        :term:`TCRs<TCR>` with more than one pair of α- and β (or γ- and δ) chains.
        See :term:`Dual IR`.

    Multichain-cell
        Cells with more than two pairs of :term:`VJ<V(D)J>` and
        :term:`VDJ<V(D)J>` sequences that do not fit into the :term:`Dual IR`
        model. These are usually rare and could be explained by
        doublets/multiplets, i.e. two ore more cells that were captured in
        the same droplet.

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
        A :term:`IR` chain is called *orphan*, if its corresponding counterpart
        has not been detected. For instance, if a cell has only a :term:`VJ<V(D)J>` chain,
        (e.g. TCR-alpha), but no :term:`VDJ<V(D)J>` chain (e.g. TCR-beta),
        the cell will be flagged as "Orphan VJ".

        Orphan chains are most likely the effect of stochastic dropouts due to
        sequencing inefficiencies.

        See also :func:`scirpy.tl.chain_qc`.

    UMI
        Unique molecular identifier. Some single-cell RNA-seq protocols
        label each RNA with a unique barcode prior to PCR-amplification to mitigate
        PCR bias. With these protocols, UMI-counts replace the read-counts
        generally used with RNA-seq.

    Productive chain
        Productive chains are :term:`IR` chains with a :term:`CDR3` sequence that produces
        a functional peptide. Scirpy relies on the preprocessing tools (e.g.
        CellRanger or TraCeR) for flagging non-productive chains.
        Typically chains are flagged as non-productive if they contain
        a stop codon or are not within the reading frame.

    Receptor type
        Classification of immune receptors into :term:`BCR` and :term:`TCR`.

        See also :func:`scirpy.tl.chain_qc`.

    Receptor subtype
        More fine-grained classification of the :term:`receptor type<Receptor type>`
        into

        * α/β T cells
        * γ/δ T cells
        * IG-heavy/IG-κ B cells
        * IG-heavy/IG-λ B cells

        See also :func:`scirpy.tl.chain_qc`.


    IR
        Immune receptor.

    BCR
        B-cell receptor. A BCR consiste of two Immunoglobulin (IG) heavy chains and
        two IG light chains. The two light chains contain a variable region, which is
        responsible for antigen recognition.

        .. figure:: img/bcr.jpg
           :width: 310px

           Image By CNX `OpenStax <http://cnx.org/contents/GFy_h8cu@10.53:rZudN6XP@2/Introduction>`__
           under the `CC BY-4.0 <https://creativecommons.org/licenses/by/4.0/deed.en>`__ license,
           obtained from `wikimedia commons <https://commons.wikimedia.org/w/index.php?curid=49935883>`__

    Dual IR
        :term:`IRs<IR>` with more than one pair of :term:`VJ<V(D)J>` and
        :term:`VDJ<V(D)J>` sequences. While this was
        previously thought to be impossible due to the mechanism of allelic exclusion
        (:cite:`Brady2010-gh`), there is an increasing amound of evidence for a *bona fide*
        dual-IR population (:cite:`Schuldt2019`, :cite:`Ji2010-bn`, :cite:`Vettermann2010`).

        For more information on how *Scirpy* handles dual IRs, see the
        page about our :ref:`IR model<receptor-model>`.

    AIRR
        Adaptive Immune Receptor Repertoire.
        See also the `AIRR community <https://www.antibodysociety.org/the-airr-community/>`_.

        Within the Scirpy documentation, we simply speak of :term:`immune receptors (IR)<IR>`.

    Chain locus
        Scirpy supports all valid `IGMT locus names <http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGTnomenclature.html>`_:

        Loci with a :term:`VJ<V(D)J>` junction:
            * `TRA` (T-cell receptor alpha)
            * `TRG` (T-cell receptor gamma)
            * `IGL` (Immunoglobulin lambda)
            * `IGK` (Immunoglobulin kappa)

        Loci with a :term:`VDJ<V(D)J>` junction:
            * `TRB` (T-cell receptor beta)
            * `TRD` (T cell receptor delta)
            * `IGH` (Immunoglobulin heavy chain)

    IG
        Immunoglobulin

    Alellically included B-cells
        A B cell with two pairs of :term:`IG` chains. See :term:`Dual IR`.
