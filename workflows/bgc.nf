/*
 * BGC Analysis Block
 *
 * antiSMASH BGC detection on internal samples.
 * When subset_analysis_enable is true, also runs taxonomy-based grouping,
 * BiG-SCAPE similarity network analysis, and network visualization.
 */

include { ANTISMASH } from '../modules/bgc/antismash.nf'
include { COLLECT_SUBSET_BGCS } from '../modules/bgc/collect_subset_bgcs.nf'
include { BIGSCAPE } from '../modules/bgc/bigscape.nf'
include { BIGSCAPE_NETWORK_PLOT } from '../modules/bgc/bigscape_network_plot.nf'

workflow BGC_BLOCK {

    take:
    ch_internal_fna      // tuple(sample_id, fna)
    ch_internal_gff      // tuple(sample_id, gff)
    ch_pipeline_summary  // path (value channel)

    main:
    // --- antiSMASH BGC detection (always runs) ---
    ch_internal_fna
        .join(ch_internal_gff)
        .map { sample_id, fna, gff -> tuple(sample_id, fna, gff) }
        .set { ch_antismash_input }

    ANTISMASH(ch_antismash_input)

    // --- BiG-SCAPE taxonomy grouping + networks (requires subset_analysis_enable) ---
    if (params.subset_analysis_enable) {
        ANTISMASH.out.results_dir
            .map { _sample_id, dir -> dir }
            .collect()
            .set { ch_all_antismash_dirs }

        COLLECT_SUBSET_BGCS(
            ch_all_antismash_dirs,
            ch_pipeline_summary
        )

        def pfam_file = file(params.pfam_db)

        COLLECT_SUBSET_BGCS.out.group_dirs
            .flatten()
            .filter { dir -> dir.isDirectory() }
            .map { dir -> tuple(dir.name, dir) }
            .set { ch_bigscape_input }

        BIGSCAPE(ch_bigscape_input, pfam_file)

        BIGSCAPE_NETWORK_PLOT(
            BIGSCAPE.out.results,
            ch_pipeline_summary
        )
    }

    emit:
    antismash_results = ANTISMASH.out.results_dir
}
