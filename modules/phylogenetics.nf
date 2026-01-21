/*
 * Phylogenetic Analysis Subworkflow
 * Workflow: Core gene extraction → MAFFT → IQ-TREE 3
 *
 * Takes panta_dir from PANGENOMICS workflow as input
 * For batch analysis of bacterial genomes (~50-500 genomes)
 */

include { EXTRACT_CORE_GENES } from './phylogenetics/extract_core_genes.nf'
include { ALIGN_CORE_GENES } from './phylogenetics/align_core_genes.nf'
include { CONCATENATE_ALIGNMENT } from './phylogenetics/concatenate_alignment.nf'
include { IQTREE3 } from './phylogenetics/iqtree3.nf'
include { TREE_VISUALIZATION } from './phylogenetics/tree_visualization.nf'

workflow PHYLOGENETICS {
    take:
    panta_channel       // tuple(cohort_id, panta_dir) - from PANGENOMICS workflow
    sample_count_ch     // tuple(cohort_id, sample_count) - number of samples in pangenome
    internal_genus_map  // path to TSV file with internal sample_id -> genus mapping

    main:

    // =========================================================================
    // Step 1: Extract core genes from pangenome (from PANGENOMICS)
    // =========================================================================

    panta_channel
        .join(sample_count_ch)
        .set { extract_input }

    EXTRACT_CORE_GENES(extract_input)
    
    // =========================================================================
    // Step 2: Align core genes with MAFFT
    // =========================================================================

    EXTRACT_CORE_GENES.out.core_gene_seqs
        .join(EXTRACT_CORE_GENES.out.sample_count)
        .set { align_input }

    ALIGN_CORE_GENES(align_input)

    // =========================================================================
    // Step 3: Concatenate alignments into supermatrix
    // =========================================================================

    ALIGN_CORE_GENES.out.aligned_genes
        .join(ALIGN_CORE_GENES.out.sample_count)
        .set { concat_input }

    CONCATENATE_ALIGNMENT(concat_input)

    // =========================================================================
    // Step 4: Build phylogenetic tree with IQ-TREE 3
    // =========================================================================

    CONCATENATE_ALIGNMENT.out.core_alignment
        .join(CONCATENATE_ALIGNMENT.out.partition_file)
        .join(CONCATENATE_ALIGNMENT.out.sample_count)
        .set { tree_input }

    IQTREE3(tree_input)

    // =========================================================================
    // Step 5: Generate tree visualizations (PNG, SVG, PDF, HTML)
    // =========================================================================

    // Combine tree with genus map for visualization
    IQTREE3.out.consensus_tree
        .combine(internal_genus_map)
        .set { ch_viz_input }

    TREE_VISUALIZATION(ch_viz_input)

    emit:
    // Core gene extraction outputs
    core_genes = EXTRACT_CORE_GENES.out.core_gene_list

    // Alignment outputs
    core_alignment = CONCATENATE_ALIGNMENT.out.core_alignment
    alignment_stats = CONCATENATE_ALIGNMENT.out.alignment_stats
    partition_file = CONCATENATE_ALIGNMENT.out.partition_file

    // Tree outputs (IQ-TREE 3)
    tree = IQTREE3.out.tree
    consensus_tree = IQTREE3.out.consensus_tree
    iqtree_report = IQTREE3.out.iqtree_report
    iqtree_stats = IQTREE3.out.iqtree_stats

    // Tree visualization outputs
    tree_png = TREE_VISUALIZATION.out.png
    tree_svg = TREE_VISUALIZATION.out.svg
    tree_pdf = TREE_VISUALIZATION.out.pdf
    taxonomy_map = TREE_VISUALIZATION.out.taxonomy_map
}
