/*
 * Phylogenetics Block
 *
 * Core gene extraction, alignment, IQ-TREE phylogeny, tree visualization,
 * and optional genus/species subtree pruning + visualization.
 */

include { PHYLOGENETICS } from '../modules/phylogenetics/phylogenetics.nf'
include { BUILD_SUBTREE_GROUPS } from '../modules/phylogenetics/build_subtree_groups.nf'
include { SUBTREE_VISUALIZATION } from '../modules/phylogenetics/subtree_visualization.nf'

workflow PHYLOGENETICS_BLOCK {

    take:
    ch_panta_dir          // tuple(cohort_id, panta_dir)
    ch_sample_count       // tuple(cohort_id, sample_count)
    ch_internal_taxonomy  // tuple(sample_id, family, genus, species)
    ch_taxonomy_csv       // path (value channel)
    ch_pipeline_summary   // path (value channel)

    main:
    // Build genus map for tree visualization
    ch_genus_map_file = ch_internal_taxonomy
        .map { sample_id, _family, genus, _species -> [sample_id, genus] }
        .collectFile(name: 'internal_genus_map.tsv', storeDir: "${params.outdir}/phylogenetics") { sample_id, genus ->
            "${sample_id}\t${genus}\n"
        }

    // --- Core phylogenetics pipeline ---
    PHYLOGENETICS(
        ch_panta_dir,
        ch_sample_count,
        ch_genus_map_file
    )

    // --- Subtree visualization (prune main tree by genus/species) ---
    if (params.subset_analysis_enable) {
        // Use ML tree (.treefile) instead of consensus tree (.contree)
        // to avoid inflated branch lengths from bootstrap averaging
        BUILD_SUBTREE_GROUPS(
            PHYLOGENETICS.out.tree.map { _id, tree -> tree }.first(),
            ch_taxonomy_csv,
            ch_pipeline_summary
        )

        ch_subtree_groups = BUILD_SUBTREE_GROUPS.out.tip_files
            .flatten()
            .splitCsv(sep: '\t', strip: true)
            .map { row -> tuple(row[0], row[1]) }

        SUBTREE_VISUALIZATION(
            ch_subtree_groups,
            PHYLOGENETICS.out.tree.map { _id, tree -> tree }.first(),
            ch_taxonomy_csv,
            ch_pipeline_summary
        )
    }

    emit:
    tree           = PHYLOGENETICS.out.tree
    consensus_tree = PHYLOGENETICS.out.consensus_tree
    tree_png       = PHYLOGENETICS.out.tree_png
}
