#!/usr/bin/env nextflow

/*
========================================================================================
    SPAdes-Bakta-Panta-PASA-RagTag Pipeline
========================================================================================
    A Nextflow pipeline for bacterial genome assembly, annotation, pangenome analysis,
    and reference-guided scaffolding.

    Pipeline steps:
    Optional preprocessing (--run_preprocessing):
      1. FastQC:     Initial read quality assessment
      2. fastp:      Read quality control and trimming
      3. Kraken2:    Taxonomic classification of reads
      4. KrakenTools: Filter for Lactobacillaceae and unclassified reads
      5. FastQC:     Final read quality assessment

    Core pipeline:
    1. SPAdes:    De novo genome assembly from raw or filtered reads
    2. QC_SPADES: Quality control with CheckM2 (completeness/contamination)
    3. HQ_FILTER: Filter samples based on CheckM2 quality metrics
                  Only high-quality samples (>90% complete, <5% contamination) proceed
    4. CLASSIFICATION: Taxonomic classification with sourmash (only HQ samples)
                       - Identifies genus/family
                       - Filters non-Lactobacillaceae samples
                       - Downloads reference genomes for the genus
                       - Runs CheckM2/Bakta on reference genomes
                       - Selects best references via skani ANI
    5. Bakta:     Bacterial genome annotation (HQ + Lactobacillaceae samples only)
    6. Panta:     Pangenome analysis comparing sample against reference genomes
    7. PASA:      Genome scaffolding using pangenome synteny information
    8. QC + PASA Filter: Intelligent scaffold selection (optional)

    Optional subworkflows:
    - Pangenomics: Cohort-level pangenome analysis (Panta + EggNOG + antiSMASH)
    - Phylogenetics: Core gene alignment + IQ-TREE 3 phylogeny (requires Pangenomics)

    Author: Created for BGC-link project
    Date: November 2025
========================================================================================
*/

nextflow.enable.dsl=2

/*
========================================================================================
    PARAMETER VALIDATION
========================================================================================
*/

// Parameters are defined in nextflow.config

/*
========================================================================================
    INCLUDE MODULES
========================================================================================
*/

include { PREPROCESSING } from './modules/preprocessing.nf'
include { SPADES } from './modules/spades'
include { BAKTA } from './modules/bakta'
include { PANTA } from './modules/panta'
include { PASA } from './modules/pasa'
include { QC } from './modules/qc.nf'
include { QC as QC_SPADES } from './modules/qc.nf'
include { QC as QC_STANDARD } from './modules/qc.nf'
include { QC as QC_SENSITIVE } from './modules/qc.nf'
include { CLASSIFICATION } from './modules/classification.nf'
include { HQ_FILTER } from './modules/hq_filter.nf'
include { PASA_FILTER } from './modules/pasa_filter'
// REFERENCE_SELECTOR now integrated into CLASSIFICATION workflow

// Import subworkflows
include { PHYLOGENETICS } from './modules/phylogenetics.nf'
include { PANGENOMICS } from './modules/pangenomics.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Check mandatory parameters
    if (!params.input) {
        error "Please provide an input samplesheet with --input"
    }
    
    
    // Print pipeline parameters
    log.info """\
        ===================================
        BGC-link Pipeline - Assembly, Annotation, and QC
        ===================================
        Input samplesheet : ${params.input}
        Output directory  : ${params.outdir}
        Sourmash classify : ${params.sourmash_classify}
        HQ Filter         : Enabled (CheckM2-based)
        HQ min complete   : >${params.filter_min_completeness}%
        HQ max contam     : <${params.filter_max_contamination}%
        Family filter     : ${params.filter_family} (${params.filter_target_family})
        Panta identity    : ${params.panta_identity}
        Panta core        : ${params.panta_core}
        Panta soft        : ${params.panta_soft}
        Phylogeny enabled : ${params.phylo_enable}
        Pangenome identity: ${params.pangenome_identity}
        Phylogeny model   : ${params.phylo_model}
        Bootstrap reps    : ${params.phylo_bootstrap_replicates}
        Pangenomics       : ${params.pangenomics_enable}
        PASA mode         : ${params.pasa_mode}
        PASA filter       : ${params.pasa_filter}
        Auto references   : ${params.sourmash_classify ? 'Enabled (via classification)' : 'Disabled'}
        Database source   : ${params.database_source}
        Run preprocessing : ${params.run_preprocessing}
        QC after PASA     : ${params.qc_after_pasa}
        ===================================
        """.stripIndent()
    
    // Read samplesheet
    // Expected format: sample_id,read1,read2,reference_dir
    // Note: reference_dir can be empty when sourmash_classify=true (references downloaded automatically)
    channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample_id
            def read1 = file(row.read1, checkIfExists: true)
            def read2 = file(row.read2, checkIfExists: true)

            // Handle reference directory - can be empty when using sourmash classification
            def ref_dir_value = row.reference_dir?.trim()
            def reference_dir
            if (params.sourmash_classify) {
                // When classification is enabled, reference_dir is optional (will be auto-downloaded)
                reference_dir = ref_dir_value ? file(ref_dir_value) : file('NO_REF_DIR')
            } else {
                // When classification is disabled, reference_dir is required
                if (!ref_dir_value) {
                    error "Sample ${sample_id}: reference_dir is required when sourmash_classify=false"
                }
                reference_dir = file(ref_dir_value, checkIfExists: true)
            }

            tuple(sample_id, [read1, read2], reference_dir)
        }
        .set { ch_samples }
    
    // Split channel for parallel processing
    ch_samples
        .multiMap { sample_id, reads, reference_dir ->
            reads: tuple(sample_id, reads)
            references: tuple(sample_id, reference_dir)
        }
        .set { ch_split }
    
    // Optional preprocessing steps
    if (params.run_preprocessing) {
        // STEP 1: Complete preprocessing workflow with quality assessment
        PREPROCESSING(ch_split.reads)
        
        // Use filtered and quality-controlled reads for assembly
        ch_assembly_reads = PREPROCESSING.out.filtered_reads
    } else {
        // Use original reads for assembly
        ch_assembly_reads = ch_split.reads
    }
    
    // STEP 4: Run SPAdes assembly
    SPADES(ch_assembly_reads)

    //=========================================================================
    // STEP 5: Quality Control with CheckM2 (always run for HQ filtering)
    //=========================================================================
    QC_SPADES(SPADES.out.contigs, "spades")

    //=========================================================================
    // STEP 6: High Quality Filter based on CheckM2 results
    // Only samples passing completeness/contamination thresholds proceed
    //=========================================================================
    // Prepare input for HQ filter: join contigs with CheckM2 JSON results
    ch_hq_filter_input = SPADES.out.contigs
        .join(QC_SPADES.out.checkm2_json)
        .map { sample_id, contigs, checkm2_json ->
            tuple(sample_id, contigs, checkm2_json)
        }

    HQ_FILTER(ch_hq_filter_input)

    // Only HQ samples proceed to classification
    ch_hq_contigs = HQ_FILTER.out.hq_genome

    //=========================================================================
    // STEP 7: Taxonomic Classification (only for HQ samples)
    // Includes Lactobacillaceae family filter via sourmash
    //=========================================================================
    if (params.sourmash_classify) {
        // Only classify high-quality samples
        CLASSIFICATION(ch_hq_contigs)
        ch_classification_results = CLASSIFICATION.out.classification
        ch_taxonomy_results = CLASSIFICATION.out.taxonomy
        ch_selected_references = CLASSIFICATION.out.selected_references
    } else {
        ch_classification_results = channel.empty()
        ch_taxonomy_results = channel.empty()
        ch_selected_references = channel.empty()
    }

    // Samples that pass both HQ filter AND Lactobacillaceae filter proceed to Bakta
    // When classification is enabled, use samples that made it through classification
    // When disabled, use all HQ samples
    if (params.sourmash_classify) {
        // Get sample IDs that passed classification (have selected references)
        ch_filtered_contigs = ch_hq_contigs
            .join(ch_selected_references)
            .map { sample_id, contigs, _references -> tuple(sample_id, contigs) }
    } else {
        ch_filtered_contigs = ch_hq_contigs
    }

    //=========================================================================
    // STEP 8: Run Bakta annotation on filtered contigs
    //=========================================================================
    BAKTA(ch_filtered_contigs)
    
    // Prepare reference GFF files for Panta
    if (params.sourmash_classify) {
        // Use automatically selected references from REFERENCE_SELECTOR
        ch_selected_references
            .map { sample_id, selected_csv ->
                // Read the CSV file to get GFF3 paths
                // CSV columns: accession,species,ani,af_query,completeness,contamination,quality_score,composite_score,quality_source,selection_reason,ref_file,gff_file
                def gff_files = []
                selected_csv.readLines().drop(1).each { line ->  // Skip header
                    def columns = line.split(',')
                    if (columns.size() >= 12) {  // Ensure we have the gff_file column
                        def gff_path = columns[11]  // gff_file is column 12 (0-indexed = 11)
                        if (gff_path && !gff_path.isEmpty()) {
                            def gff_file = file(gff_path)
                            if (gff_file.exists()) {
                                gff_files << gff_file
                            }
                        }
                    }
                }
                tuple(sample_id, gff_files)
            }
            .set { ch_reference_gffs }
    } else {
        // Use user-provided reference directories from samplesheet
        ch_split.references
            .map { sample_id, reference_dir ->
                // Find all GFF3 files in reference subdirectories
                def gff_files = []
                if (reference_dir.isDirectory()) {
                    reference_dir.listFiles().each { subdir ->
                        if (subdir.isDirectory()) {
                            subdir.listFiles().findAll { file -> file.name.endsWith('.gff3') }.each { gff ->
                                gff_files << gff
                            }
                        }
                    }
                }
                tuple(sample_id, gff_files)
            }
            .set { ch_reference_gffs }
    }
    
    // Combine sample GFF with reference GFFs
    BAKTA.out.gff
        .join(ch_reference_gffs)
        .map { sample_id, sample_gff, reference_gffs ->
            tuple(sample_id, sample_gff, reference_gffs)
        }
        .set { ch_panta_input }
    
    // STEP 6: Run Panta pangenome analysis
    PANTA(ch_panta_input)
    
    // STEP 7: Run PASA scaffolding using PANTA (sample vs references)
    // Combine PANTA outputs with SPADES outputs for PASA
    PANTA.out.gene_annotation
        .join(PANTA.out.gene_position)
        .join(PANTA.out.clusters)
        .join(PANTA.out.samples)
        .join(SPADES.out.contigs)
        .join(SPADES.out.graph)
        .join(SPADES.out.paths)
        .set { ch_pasa_input }
    
    
    PASA(ch_pasa_input)
    
    // STEP 8: Select best PASA scaffold using BUSCO
    // Combine PASA scaffold outputs
    PASA.out.standard_scaffold
        .join(PASA.out.sensitive_scaffold, remainder: true)
        .map { sample_id, standard, sensitive ->
            // Handle cases where one or both scaffolds may be missing
            def validStandard = standard != null && standard.exists() && standard.size() > 0
            def validSensitive = sensitive != null && sensitive.exists() && sensitive.size() > 0

            if (!validStandard && !validSensitive) {
                log.warn "Sample ${sample_id}: No valid PASA scaffolds generated - skipping QC analysis"
                return null
            }

            if (!validStandard) {
                log.warn "Sample ${sample_id}: Standard PASA scaffold empty - will only use sensitive scaffold"
            }

            if (!validSensitive) {
                log.warn "Sample ${sample_id}: Sensitive PASA scaffold empty - will only use standard scaffold"
            }

            tuple(sample_id, standard ?: file('NO_FILE'), sensitive ?: file('NO_FILE'))
        }
        .filter { tuple -> tuple != null }
        .set { ch_pasa_scaffolds }
    
    // Optional: QC analysis and intelligent PASA selection
    if (params.pasa_filter || params.qc_after_pasa) {
        // Run QC on standard PASA scaffolds - only if file is valid
        QC_STANDARD(
            ch_pasa_scaffolds
                .filter { _sample_id, std, _sens -> std.exists() && std.size() > 0 }
                .map { sample_id, std, _sens -> tuple(sample_id, std) },
            "pasa_standard"
        )
        ch_qc_standard = QC_STANDARD.out.json_reports

        // Run QC on sensitive PASA scaffolds - only if file is valid
        QC_SENSITIVE(
            ch_pasa_scaffolds
                .filter { _sample_id, _std, sens -> sens.exists() && sens.size() > 0 }
                .map { sample_id, _std, sens -> tuple(sample_id, sens) },
            "pasa_sensitive"
        )
        ch_qc_sensitive = QC_SENSITIVE.out.json_reports
        
        // Use QC results to select the best scaffold (only if pasa_filter is enabled)
        if (params.pasa_filter) {
            // Use left outer join to handle cases where only one PASA scaffold exists
            // remainder: true keeps samples even if they don't have a matching QC result
            ch_pasa_scaffolds
                .join(ch_qc_standard, remainder: true)
                .join(ch_qc_sensitive, remainder: true)
                .map { sample_id, std, sens, qc_std, qc_sens ->
                    // Replace null QC results with placeholder empty files
                    def qc_std_file = qc_std ?: file('NO_QC_FILE')
                    def qc_sens_file = qc_sens ?: file('NO_QC_FILE')
                    tuple(sample_id, std, sens, qc_std_file, qc_sens_file)
                }
                .set { ch_pasa_filter_input }

            PASA_FILTER(ch_pasa_filter_input)
            ch_selected_scaffolds = PASA_FILTER.out.selected_scaffold
        } else {
            // Just run QC for reporting, but default to sensitive scaffolds
            ch_selected_scaffolds = ch_pasa_scaffolds.map { sample_id, _standard, sensitive ->
                tuple(sample_id, sensitive)
            }
        }
    } else {
        // Default to sensitive scaffolds if no QC analysis
        ch_selected_scaffolds = ch_pasa_scaffolds.map { sample_id, _standard, sensitive ->
            tuple(sample_id, sensitive)
        }
    }

    //=========================================================================
    // OPTIONAL: Pangenomics analysis - runs after PASA filter
    // Workflow: BAKTA → PANTA + ANTISMASH (parallel) → EGGNOG → Plotting
    //=========================================================================
    if (params.pangenomics_enable) {
        // Pangenomics subworkflow:
        // 1. Re-annotates PASA scaffolds with Bakta (per-sample)
        // 2. Runs Panta on all .gff3 files (cohort, collects all samples)
        // 3. Runs antiSMASH on .fna + .gff3 (per-sample, parallel with Panta)
        // 4. Runs EggNOG on Panta representative proteins (after Panta)
        // 5. Generates combined plots (after EggNOG + antiSMASH)

        PANGENOMICS(ch_selected_scaffolds)

        //=========================================================================
        // OPTIONAL: Phylogenetic analysis - uses pangenome from PANGENOMICS
        // Workflow: Extract core genes → MAFFT alignment → IQ-TREE 3
        //=========================================================================
        if (params.phylo_enable) {
            // PHYLOGENETICS takes panta_dir and sample_count from PANGENOMICS
            PHYLOGENETICS(
                PANGENOMICS.out.panta_dir,
                PANGENOMICS.out.panta_sample_count
            )
        }
    }

    // Final outputs are available through the process outputs
}

// Workflow completion handlers
def printCompletionMessage() {
    log.info """
    ===================================
    Pipeline completed!
    Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    Output directory: ${params.outdir}
    ===================================
    """.stripIndent()
}

def printErrorMessage() {
    log.error """
    ===================================
    Pipeline failed!
    Error: ${workflow.errorMessage}
    ===================================
    """.stripIndent()
}