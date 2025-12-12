/*
 * Pangenomics Analysis Subworkflow
 *
 * Workflow:
 *   1. BAKTA_PANGENOMICS: Re-annotate PASA scaffolds (per-sample)
 *   2. PANTA_COHORT: Build pangenome from all .gff3 files (cohort-level)
 *   3. EGGNOG: Annotate Panta representative proteins (cohort-level, after Panta) - optional
 *   3. ANTISMASH: Detect BGCs using .fna + .gff3 (per-sample, parallel with Panta)
 *   4. PANGENOME_PLOTS: Generate visualizations (after Panta + EggNOG) - requires EggNOG
 *
 * Execution order:
 *   BAKTA_PANGENOMICS (per-sample)
 *           ↓
 *   ┌───────┴───────┐
 *   ↓               ↓
 *   PANTA_COHORT    ANTISMASH (per-sample)
 *   ↓
 *   EGGNOG (optional)
 *   ↓
 *   PANGENOME_PLOTS (optional, requires EggNOG)
 */

include { BAKTA_PANGENOMICS } from './pangenomics/bakta_pangenomics.nf'
include { PANTA_COHORT } from './pangenomics/panta_cohort.nf'
include { EGGNOG } from './pangenomics/eggnog.nf'
include { ANTISMASH } from './pangenomics/antismash.nf'
include { PANGENOME_PLOTS } from './pangenomics/plotting.nf'

workflow PANGENOMICS {
    take:
    scaffold_channel    // tuple(sample_id, fasta) - PASA selected scaffolds

    main:

    // =========================================================================
    // Step 1: Re-annotate PASA scaffolds with Bakta (per-sample)
    // Generates .gff3, .fna, .faa for downstream tools
    // =========================================================================

    BAKTA_PANGENOMICS(scaffold_channel)

    // =========================================================================
    // Step 2a: Collect all GFF3 files for Panta pangenome analysis (cohort)
    // =========================================================================

    BAKTA_PANGENOMICS.out.gff
        .map { sample_id, gff -> gff }
        .collect()
        .map { gff_list ->
            def count = gff_list.size()
            log.info "Collected ${count} genomes for pangenomics analysis"
            tuple("pangenomics_${count}_samples", gff_list, count)
        }
        .set { all_gffs }

    PANTA_COHORT(all_gffs)

    // =========================================================================
    // Step 2b: Run antiSMASH on each sample (per-sample, parallel with Panta)
    // Uses .fna + .gff3 from Bakta
    // =========================================================================

    BAKTA_PANGENOMICS.out.fna
        .join(BAKTA_PANGENOMICS.out.gff)
        .map { sample_id, fna, gff ->
            tuple(sample_id, fna, gff)
        }
        .set { ch_antismash_input }

    ANTISMASH(ch_antismash_input)

    // =========================================================================
    // Step 3: Run EggNOG on Panta representative proteins (cohort, after Panta)
    // Only runs if eggnog_enable is true (requires downloaded database)
    // =========================================================================

    // Initialize empty channels for conditional outputs
    ch_eggnog_annotations = Channel.empty()
    ch_eggnog_hits = Channel.empty()
    ch_eggnog_orthologs = Channel.empty()
    ch_plots_dir = Channel.empty()
    ch_pangenome_summary = Channel.empty()
    ch_pangenome_statistics = Channel.empty()

    if (params.eggnog_enable) {
        EGGNOG(PANTA_COHORT.out.panta_dir)

        ch_eggnog_annotations = EGGNOG.out.annotations
        ch_eggnog_hits = EGGNOG.out.hits
        ch_eggnog_orthologs = EGGNOG.out.seed_orthologs

        // =========================================================================
        // Step 4: Run plotting after Panta + EggNOG complete
        // Uses Rtab file from Panta and COG annotations from EggNOG
        // =========================================================================

        // Join Rtab file with EggNOG annotations for plotting
        PANTA_COHORT.out.rtab
            .join(EGGNOG.out.annotations)
            .map { cohort_id, rtab_file, eggnog_annot ->
                tuple(cohort_id, rtab_file, eggnog_annot)
            }
            .set { ch_plotting_input }

        PANGENOME_PLOTS(ch_plotting_input)

        ch_plots_dir = PANGENOME_PLOTS.out.plots_dir
        ch_pangenome_summary = PANGENOME_PLOTS.out.summary
        ch_pangenome_statistics = PANGENOME_PLOTS.out.statistics
    } else {
        log.info "EggNOG annotation disabled (params.eggnog_enable = false)"
        log.info "Skipping PANGENOME_PLOTS (requires EggNOG annotations)"
    }

    emit:
    // Bakta outputs (per-sample)
    bakta_gff = BAKTA_PANGENOMICS.out.gff
    bakta_fna = BAKTA_PANGENOMICS.out.fna
    bakta_faa = BAKTA_PANGENOMICS.out.faa

    // Panta outputs (cohort)
    panta_dir = PANTA_COHORT.out.panta_dir
    panta_sample_count = PANTA_COHORT.out.sample_count
    gene_annotation = PANTA_COHORT.out.gene_annotation
    gene_position = PANTA_COHORT.out.gene_position
    clusters = PANTA_COHORT.out.clusters
    samples = PANTA_COHORT.out.samples

    // EggNOG outputs (cohort) - empty if disabled
    eggnog_annotations = ch_eggnog_annotations
    eggnog_hits = ch_eggnog_hits
    eggnog_orthologs = ch_eggnog_orthologs

    // antiSMASH outputs (per-sample)
    antismash_results = ANTISMASH.out.results_dir
    antismash_genbank = ANTISMASH.out.genbank
    antismash_summary = ANTISMASH.out.summary

    // Plotting outputs (cohort) - empty if EggNOG disabled
    plots = ch_plots_dir
    pangenome_summary = ch_pangenome_summary
    pangenome_statistics = ch_pangenome_statistics
}
