#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    BGC-link Analysis Pipeline
========================================================================================
    Downstream analysis pipeline for pangenomics and phylogenetics.
    Takes PASA scaffold FASTAs from assembly pipeline, annotates with Bakta,
    then performs cohort-level pangenomics and phylogenetics.

    Pipeline steps:
      1. BAKTA           - Annotate scaffolds (or load existing with --spades)
      2. CODON_QC        - Codon bias analysis (detects assembly anomalies)
      3. MERGE_EXTERNAL  - Merge with external annotations (optional)
      4. PREPARE_PANTA   - Collect GFF3s for cohort
      5. PANTA_COHORT    - Cohort-level pangenome analysis
      6. ANTISMASH       - BGC detection (parallel with PANTA)
      7. PHYLOGENETICS   - Core gene tree (after PANTA)
      8. EGGNOG          - Functional annotation (after PHYLOGENETICS)
      9. PANGENOME_PLOTS - Visualizations (after EGGNOG)

    Run with: nextflow run analysis.nf -c analysis.config --internal results/assembly
    SPAdes mode: nextflow run analysis.nf -c analysis.config --internal results/assembly --spades

    Normal mode:  reads PASA scaffolds from <internal>/SAMPLE/pasa_filter/ and runs Bakta
    SPAdes mode:  reads existing Bakta from <internal>/SAMPLE/bakta/ (skips Bakta annotation)

    Author: Created for BGC-link project
========================================================================================
*/

// Import annotation module
include { BAKTA } from './modules/bakta'

// Import codon QC workflow
include { CODON_BIAS_QC } from './modules/codon.nf'

// Import pangenomics modules
include { PANTA_COHORT } from './modules/pangenomics/panta_cohort.nf'
include { EGGNOG } from './modules/pangenomics/eggnog.nf'
include { ANTISMASH } from './modules/pangenomics/antismash.nf'
include { PANGENOME_PLOTS } from './modules/pangenomics/plotting.nf'

// Import phylogenetics workflow
include { PHYLOGENETICS } from './modules/phylogenetics.nf'

// =========================================================================
// Parameters (defaults - override in analysis.config)
// =========================================================================

params.outdir = 'results/analysis'

// Input: Directory containing assembly pipeline results
// Will search for */pasa_filter/selected_scaffold.fasta within this directory
// Also reads pipeline_summary.tsv from this directory for genus information
params.internal = 'results/assembly_validation'

// Optional: Pre-annotated external references directory (GFF3/FNA files, skip Bakta)
params.external = null

// Optional: Use pre-existing Bakta results from assembly pipeline (skip Bakta annotation)
// When true, reads bakta results from params.internal/*/bakta/ instead of running Bakta
params.spades = false

// Enable/disable modules
params.codon_qc_enable = true
params.antismash_enable = true
params.eggnog_enable = true
params.phylogenetics_enable = true

// =========================================================================
// Main Workflow
// =========================================================================

workflow {

    log.info """
    ===================================
    BGC-link Analysis Pipeline
    ===================================
    Input sources:
      - Internal: ${params.internal}
      - External: ${params.external ?: 'None'}
      - SPAdes mode: ${params.spades} ${params.spades ? '(using pre-existing Bakta from ' + params.internal + '/*/bakta/)' : ''}

    Modules enabled:
      - BAKTA:           ${params.spades ? 'SKIPPED (using existing)' : 'yes (annotate scaffolds)'}
      - CODON_QC:        ${params.codon_qc_enable}
      - PANTA_COHORT:    always
      - ANTISMASH:       ${params.antismash_enable} (internal samples only)
      - EGGNOG:          ${params.eggnog_enable}
      - PANGENOME_PLOTS: ${params.eggnog_enable} (requires EGGNOG)
      - PHYLOGENETICS:   ${params.phylogenetics_enable}

    Output directory: ${params.outdir}
    ===================================
    """.stripIndent()

    // =========================================================================
    // Step 1: Collect inputs and run/load Bakta annotation
    // =========================================================================

    // Read pipeline_summary.tsv as a lookup table for taxonomy info
    ch_internal_taxonomy = channel
        .fromPath("${params.internal}/pipeline_summary.tsv", checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.sample_id,
                row.family ?: 'Unknown',
                row.genus ?: 'Unknown',
                row.species ?: 'Unknown'
            )
        }

    if (params.spades) {
        // =====================================================================
        // SPAdes mode: Use pre-existing Bakta results from assembly pipeline
        // =====================================================================
        log.info "SPAdes mode enabled: Loading pre-existing Bakta annotations from ${params.internal}/*/bakta/"

        // Read existing GFF3 files from assembly bakta output
        ch_internal_gff = channel
            .fromPath("${params.internal}/*/bakta/*.gff3", checkIfExists: true)
            .map { gff ->
                // Extract sample_id from filename (remove .gff3 extension)
                def sample_id = gff.baseName
                tuple(sample_id, gff)
            }

        // Read existing FNA files from assembly bakta output
        ch_internal_fna = channel
            .fromPath("${params.internal}/*/bakta/*.fna", checkIfExists: true)
            .map { fna ->
                def sample_id = fna.baseName
                tuple(sample_id, fna)
            }

        // Read existing FFN files from assembly bakta output (for codon QC)
        ch_internal_ffn = channel
            .fromPath("${params.internal}/*/bakta/*.ffn", checkIfExists: true)
            .map { ffn ->
                def sample_id = ffn.baseName
                tuple(sample_id, ffn)
            }

    } else {
        // =====================================================================
        // Normal mode: Run Bakta annotation on PASA scaffolds
        // =====================================================================

        // Read PASA scaffolds from internal directory
        // Searches for */pasa_filter/selected_scaffold.fasta within the internal directory
        ch_scaffolds = channel
            .fromPath("${params.internal}/*/pasa_filter/selected_scaffold.fasta", checkIfExists: true)
            .map { fasta ->
                // Extract sample_id from path: .../{sample_id}/pasa_filter/...
                def sample_id = fasta.parent.parent.name
                tuple(sample_id, fasta)
            }

        // Run Bakta annotation on scaffolds
        BAKTA(ch_scaffolds)

        // Get GFF3 and FNA outputs from Bakta
        ch_internal_gff = BAKTA.out.gff
        ch_internal_fna = BAKTA.out.fna
        ch_internal_ffn = BAKTA.out.nucleotides
    }

    // =========================================================================
    // Step 2: Codon Bias QC (optional - detects assembly anomalies)
    // =========================================================================

    if (params.codon_qc_enable) {
        // Prepare samples with taxonomy info for codon QC
        // Join BAKTA FFN outputs with taxonomy from pipeline_summary.tsv
        ch_codon_samples = ch_internal_ffn
            .join(ch_internal_taxonomy)
            .map { sample_id, ffn, family, genus, species ->
                def f = (family && family != 'NA' && family != '') ? family : 'Unknown'
                def g = (genus && genus != 'NA' && genus != '') ? genus : 'Unknown'
                def s = (species && species != 'NA' && species != '') ? species : 'Unknown'
                [[id: sample_id, family: f, genus: g, species: s], ffn]
            }

        // Run codon bias QC workflow
        // Compares samples against species-level (preferred) or genus-level references
        CODON_BIAS_QC(ch_codon_samples)
    }

    // =========================================================================
    // Step 3: Merge with external annotations (if provided)
    // =========================================================================

    // External annotations (optional) - these are pre-annotated, skip Bakta
    // Supports both flat structure (*.gff3) and nested structure (*/*.gff3)
    if (params.external) {
        ch_external_gff = channel
            .fromPath(["${params.external}/*.gff3", "${params.external}/**/*.gff3"], checkIfExists: false)
            .map { gff ->
                def sample_id = gff.baseName
                tuple(sample_id, gff)
            }
        ch_external_fna = channel
            .fromPath(["${params.external}/*.fna", "${params.external}/**/*.fna"], checkIfExists: false)
            .map { fna ->
                def sample_id = fna.baseName
                tuple(sample_id, fna)
            }
        ch_all_gff = ch_internal_gff.mix(ch_external_gff)
        ch_all_fna = ch_internal_fna.mix(ch_external_fna)
    } else {
        ch_all_gff = ch_internal_gff
        ch_all_fna = ch_internal_fna
    }

    // =========================================================================
    // Step 4: Prepare PANTA input (collect all GFF3s into cohort)
    // =========================================================================

    ch_all_gff
        .map { _sample_id, gff -> gff }
        .collect()
        .map { gff_list ->
            def count = gff_list.size()
            tuple("pangenomics_cohort", gff_list, count)
        }
        .set { ch_panta_input }

    // =========================================================================
    // Step 5: Run PANTA_COHORT (pangenome analysis)
    // =========================================================================

    PANTA_COHORT(ch_panta_input)

    // =========================================================================
    // Step 6: Run ANTISMASH (parallel with PANTA) - internal samples only
    // =========================================================================

    if (params.antismash_enable) {
        // Join FNA and GFF by sample_id - only internal samples (not external)
        ch_internal_fna
            .join(ch_internal_gff)
            .map { sample_id, fna, gff ->
                tuple(sample_id, fna, gff)
            }
            .set { ch_antismash_input }

        ANTISMASH(ch_antismash_input)
    }

    // =========================================================================
    // Step 7: Run PHYLOGENETICS (after PANTA)
    // =========================================================================

    ch_phylo_done = channel.empty()

    if (params.phylogenetics_enable) {
        // Collect internal genus mapping to a single TSV file for tree visualization
        ch_internal_taxonomy
            .map { sample_id, _family, genus, _species -> [sample_id, genus] }
            .collectFile(name: 'internal_genus_map.tsv', storeDir: "${params.outdir}/phylogenetics") { sample_id, genus ->
                "${sample_id}\t${genus}\n"
            }
            .set { ch_genus_map_file }

        PHYLOGENETICS(
            PANTA_COHORT.out.panta_dir,
            PANTA_COHORT.out.sample_count,
            ch_genus_map_file
        )
        // Use tree output as signal that phylogenetics is done
        ch_phylo_done = PHYLOGENETICS.out.tree
    }

    // =========================================================================
    // Step 8: Run EGGNOG (after PHYLOGENETICS - runs last as it's slow)
    // =========================================================================

    ch_eggnog_annotations = channel.empty()

    if (params.eggnog_enable) {
        if (params.phylogenetics_enable) {
            // Wait for phylogenetics to complete before starting EggNOG
            // Combine panta_dir with phylo completion signal
            PANTA_COHORT.out.panta_dir
                .combine(ch_phylo_done.map { _it -> true }.first())
                .map { cohort_id, panta_dir, _done -> tuple(cohort_id, panta_dir) }
                .set { ch_eggnog_input }

            EGGNOG(ch_eggnog_input)
        } else {
            // If phylogenetics disabled, run EggNOG directly after PANTA
            EGGNOG(PANTA_COHORT.out.panta_dir)
        }
        ch_eggnog_annotations = EGGNOG.out.annotations

        // =====================================================================
        // Step 9: Run PANGENOME_PLOTS (after EGGNOG)
        // =====================================================================

        PANTA_COHORT.out.rtab
            .join(EGGNOG.out.annotations)
            .map { cohort_id, rtab_file, eggnog_annot ->
                tuple(cohort_id, rtab_file, eggnog_annot)
            }
            .set { ch_plotting_input }

        PANGENOME_PLOTS(ch_plotting_input)
    }

    // Log completion
    log.info "Analysis pipeline submitted - check Nextflow output for progress"
}
