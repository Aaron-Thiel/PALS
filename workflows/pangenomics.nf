/*
 * Pangenomics Block
 *
 * EggNOG functional annotation, family-level pangenome plots,
 * and genus/species-level subset analysis.
 *
 * Auto-detects cached EggNOG results from params.outdir.
 */

include { EGGNOG } from '../modules/pangenomics/eggnog.nf'
include { PANGENOME_PLOTS } from '../modules/pangenomics/plotting.nf'
include { SUBSET_RTAB } from '../modules/pangenomics/subset_rtab.nf'
include { PANGENOME_PLOTS as SUBSET_PANGENOME_PLOTS } from '../modules/pangenomics/plotting.nf'

// HQ-only plot aliases
include { FILTER_RTAB_HQ } from '../modules/pangenomics/filter_rtab_hq.nf'
include { PANGENOME_PLOTS as HQ_PANGENOME_PLOTS } from '../modules/pangenomics/plotting.nf'
include { SUBSET_RTAB as HQ_SUBSET_RTAB } from '../modules/pangenomics/subset_rtab.nf'
include { PANGENOME_PLOTS as HQ_SUBSET_PANGENOME_PLOTS } from '../modules/pangenomics/plotting.nf'

workflow PANGENOMICS_BLOCK {

    take:
    ch_panta_dir         // tuple(cohort_id, panta_dir)
    ch_rtab              // tuple(cohort_id, rtab)
    ch_pipeline_summary  // path (value channel)
    ch_taxonomy_csv      // path (value channel)

    main:
    // --- EggNOG functional annotation (auto-detect cache) ---
    def cached_eggnog = file("${params.outdir}/pangenomics/eggnog/pangenomics_cohort.emapper.annotations")
    if (cached_eggnog.exists()) {
        log.info "  [pangenomics] Loading cached EggNOG annotations from ${params.outdir}/pangenomics/eggnog/"
        ch_eggnog_annotations = channel.of(tuple("pangenomics_cohort", cached_eggnog))
    } else {
        EGGNOG(ch_panta_dir)
        ch_eggnog_annotations = EGGNOG.out.annotations
    }

    // --- Family-level pangenome plots ---
    ch_rtab
        .join(ch_eggnog_annotations)
        .map { cohort_id, rtab_file, eggnog_annot ->
            tuple(cohort_id, rtab_file, eggnog_annot)
        }
        .set { ch_plotting_input }

    PANGENOME_PLOTS(ch_plotting_input)

    // --- Genus/species-level subset analysis ---
    if (params.subset_analysis_enable) {
        SUBSET_RTAB(
            ch_rtab,
            ch_pipeline_summary,
            ch_taxonomy_csv
        )

        SUBSET_RTAB.out.rtab_files
            .flatten()
            .map { rtab -> tuple(rtab.baseName, rtab) }
            .combine(ch_eggnog_annotations.map { _id, ann -> ann })
            .map { group_id, rtab, ann -> tuple(group_id, rtab, ann) }
            .set { ch_subset_plotting }

        SUBSET_PANGENOME_PLOTS(ch_subset_plotting)
    }

    // --- HQ-only plots (optional) ---
    if (params.hq_plots_enable) {
        FILTER_RTAB_HQ(
            ch_rtab,
            ch_pipeline_summary,
            params.scaffold_min_completeness
        )

        // Family-level HQ plots
        FILTER_RTAB_HQ.out.rtab
            .join(ch_eggnog_annotations)
            .map { cohort_id, rtab_file, eggnog_annot ->
                tuple(cohort_id, rtab_file, eggnog_annot)
            }
            .set { ch_hq_plotting_input }

        HQ_PANGENOME_PLOTS(ch_hq_plotting_input)

        // HQ genus/species subsets
        if (params.subset_analysis_enable) {
            HQ_SUBSET_RTAB(
                FILTER_RTAB_HQ.out.rtab,
                ch_pipeline_summary,
                ch_taxonomy_csv
            )

            HQ_SUBSET_RTAB.out.rtab_files
                .flatten()
                .map { rtab -> tuple(rtab.baseName, rtab) }
                .combine(ch_eggnog_annotations.map { _id, ann -> ann })
                .map { group_id, rtab, ann -> tuple(group_id, rtab, ann) }
                .set { ch_hq_subset_plotting }

            HQ_SUBSET_PANGENOME_PLOTS(ch_hq_subset_plotting)
        }
    }

    emit:
    eggnog_annotations = ch_eggnog_annotations
    plots              = PANGENOME_PLOTS.out.plots_dir
}
