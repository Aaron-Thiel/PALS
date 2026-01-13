#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Custom Pangenomics + Phylogenetics Pipeline
 *
 * Takes pre-existing Bakta annotation results as input and runs:
 *   1. PANTA_COHORT - Pangenome analysis
 *   2. ANTISMASH - BGC detection (parallel with PANTA)
 *   3. PHYLOGENETICS - Core gene tree (after PANTA)
 *   4. EGGNOG - Functional annotation (after PHYLOGENETICS - runs last as it's slow)
 *   5. PANGENOME_PLOTS - Visualizations (after EGGNOG)
 */

// Import pangenomics modules (skip BAKTA_PANGENOMICS)
include { PANTA_COHORT } from './modules/pangenomics/panta_cohort.nf'
include { EGGNOG } from './modules/pangenomics/eggnog.nf'
include { ANTISMASH } from './modules/pangenomics/antismash.nf'
include { PANGENOME_PLOTS } from './modules/pangenomics/plotting.nf'

// Import phylogenetics workflow
include { PHYLOGENETICS } from './modules/phylogenetics.nf'

// =========================================================================
// Parameters
// =========================================================================

params.outdir = 'results_panphylo'

// Input paths for bakta results
params.samn_bakta = 'validation_results/SAMN*/pangenomics/bakta'
params.gca_annotations = '/home/azureuser/BGC-link/master_thesis/external_data/lactoannotated/annotations/GCA_*'

// Enable/disable modules
params.antismash_enable = true
params.eggnog_enable = true
params.phylogenetics_enable = true

// =========================================================================
// Main Workflow
// =========================================================================

workflow {

    log.info """
    ╔═══════════════════════════════════════════════════════════════════╗
    ║     Custom Pangenomics + Phylogenetics Pipeline                   ║
    ╠═══════════════════════════════════════════════════════════════════╣
    ║  Input sources:                                                   ║
    ║    - SAMN bakta: ${params.samn_bakta}
    ║    - GCA annotations: ${params.gca_annotations}
    ║                                                                   ║
    ║  Modules enabled:                                                 ║
    ║    - PANTA_COHORT: always                                         ║
    ║    - ANTISMASH: ${params.antismash_enable}
    ║    - EGGNOG: ${params.eggnog_enable} (uses DIAMOND)
    ║    - PANGENOME_PLOTS: ${params.eggnog_enable} (requires EGGNOG)
    ║    - PHYLOGENETICS: ${params.phylogenetics_enable}
    ╚═══════════════════════════════════════════════════════════════════╝
    """.stripIndent()

    // =========================================================================
    // Step 1: Collect all GFF3 files from both sources
    // =========================================================================

    // SAMN samples: validation_results/SAMN*/pangenomics/bakta/*.gff3
    ch_samn_gff = Channel
        .fromPath("${params.samn_bakta}/*.gff3")
        .map { gff ->
            def sample_id = gff.baseName  // e.g., SAMN11618738
            tuple(sample_id, gff)
        }

    // GCA samples: /path/annotations/GCA_*/GCA_*.gff3
    ch_gca_gff = Channel
        .fromPath("${params.gca_annotations}/*.gff3")
        .map { gff ->
            def sample_id = gff.baseName  // e.g., GCA_000008065.1_ASM806v1
            tuple(sample_id, gff)
        }

    // Merge all GFF3 files
    ch_all_gff = ch_samn_gff.mix(ch_gca_gff)

    // =========================================================================
    // Step 2: Collect all FNA files (for antiSMASH)
    // =========================================================================

    // SAMN samples
    ch_samn_fna = Channel
        .fromPath("${params.samn_bakta}/*.fna")
        .map { fna ->
            def sample_id = fna.baseName
            tuple(sample_id, fna)
        }

    // GCA samples
    ch_gca_fna = Channel
        .fromPath("${params.gca_annotations}/*.fna")
        .map { fna ->
            def sample_id = fna.baseName
            tuple(sample_id, fna)
        }

    // Merge all FNA files
    ch_all_fna = ch_samn_fna.mix(ch_gca_fna)

    // =========================================================================
    // Step 3: Prepare PANTA input (collect all GFF3s into cohort)
    // =========================================================================

    ch_all_gff
        .map { sample_id, gff -> gff }
        .collect()
        .map { gff_list ->
            def count = gff_list.size()
            log.info "Collected ${count} genomes for pangenomics analysis"
            tuple("pangenomics_cohort", gff_list, count)
        }
        .set { ch_panta_input }

    // =========================================================================
    // Step 4: Run PANTA_COHORT (pangenome analysis)
    // =========================================================================

    PANTA_COHORT(ch_panta_input)

    // =========================================================================
    // Step 5: Run ANTISMASH (parallel with PANTA)
    // =========================================================================

    if (params.antismash_enable) {
        // Join FNA and GFF by sample_id
        ch_all_fna
            .join(ch_all_gff)
            .map { sample_id, fna, gff ->
                tuple(sample_id, fna, gff)
            }
            .set { ch_antismash_input }

        ANTISMASH(ch_antismash_input)
    }

    // =========================================================================
    // Step 6: Run PHYLOGENETICS (after PANTA)
    // =========================================================================

    ch_phylo_done = Channel.empty()

    if (params.phylogenetics_enable) {
        PHYLOGENETICS(
            PANTA_COHORT.out.panta_dir,
            PANTA_COHORT.out.sample_count
        )
        // Use tree output as signal that phylogenetics is done
        ch_phylo_done = PHYLOGENETICS.out.tree
    }

    // =========================================================================
    // Step 7: Run EGGNOG (after PHYLOGENETICS - runs last as it's slow)
    // =========================================================================

    ch_eggnog_annotations = Channel.empty()

    if (params.eggnog_enable) {
        if (params.phylogenetics_enable) {
            // Wait for phylogenetics to complete before starting EggNOG
            // Combine panta_dir with phylo completion signal
            PANTA_COHORT.out.panta_dir
                .combine(ch_phylo_done.map { it -> true }.first())
                .map { cohort_id, panta_dir, done -> tuple(cohort_id, panta_dir) }
                .set { ch_eggnog_input }

            EGGNOG(ch_eggnog_input)
        } else {
            // If phylogenetics disabled, run EggNOG directly after PANTA
            EGGNOG(PANTA_COHORT.out.panta_dir)
        }
        ch_eggnog_annotations = EGGNOG.out.annotations

        // =====================================================================
        // Step 8: Run PANGENOME_PLOTS (after EGGNOG)
        // =====================================================================

        PANTA_COHORT.out.rtab
            .join(EGGNOG.out.annotations)
            .map { cohort_id, rtab_file, eggnog_annot ->
                tuple(cohort_id, rtab_file, eggnog_annot)
            }
            .set { ch_plotting_input }

        PANGENOME_PLOTS(ch_plotting_input)
    }
}

// =========================================================================
// Workflow completion handler
// =========================================================================

workflow.onComplete {
    log.info """
    ╔═══════════════════════════════════════════════════════════════════╗
    ║                    Pipeline Complete!                             ║
    ╠═══════════════════════════════════════════════════════════════════╣
    ║  Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    ║  Duration: ${workflow.duration}
    ║  Output directory: ${params.outdir}
    ╚═══════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}
