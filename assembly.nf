#!/usr/bin/env nextflow

/*
========================================================================================
    BGC-link Assembly Pipeline
========================================================================================
    Genome assembly, annotation, and scaffolding pipeline.

    Run with: nextflow run assembly.nf -c assembly.config --input samplesheet.csv

    For downstream analysis (pangenomics, phylogenetics), use analysis.nf

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
    9. SUMMARY:   Generate sample summary TSV with key metrics

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
include { QC as QC_SPADES } from './modules/qc.nf'
include { QC as QC_PASA_STANDARD } from './modules/qc.nf'
include { QC as QC_PASA_SENSITIVE } from './modules/qc.nf'
include { CLASSIFICATION } from './modules/classification.nf'
include { HQ_FILTER } from './modules/hq_filter.nf'
include { PASA_FILTER } from './modules/pasa_filter'

// Import summary module
include { SUMMARY } from './modules/summary.nf'

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
        BGC-link Assembly Pipeline
        ===================================
        Input samplesheet : ${params.input}
        Output directory  : ${params.outdir}

        Preprocessing     : ${params.run_preprocessing}
        Sourmash classify : ${params.sourmash_classify}
        Database source   : ${params.database_source}

        HQ Filter:
          Min completeness: >${params.filter_min_completeness}%
          Max contamination: <${params.filter_max_contamination}%
          Family filter   : ${params.filter_family} (${params.filter_target_family})

        PASA scaffolding:
          Mode            : ${params.pasa_mode}
          QC filter       : ${params.pasa_filter}
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
                // CSV columns: sample_id,accession,species,ani,af_query,completeness,contamination,quality_score,composite_score,quality_source,selection_reason,ref_file,gff_file
                def gff_files = []
                selected_csv.readLines().drop(1).each { line ->  // Skip header
                    def columns = line.split(',')
                    if (columns.size() >= 13) {  // Ensure we have the gff_file column
                        def gff_path = columns[12]  // gff_file is column 13 (0-indexed = 12)
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
    
    //=========================================================================
    // STEP 8: Run QC on PASA scaffolds (always runs for summary comparison)
    //=========================================================================

    // Run QC on standard PASA scaffolds - only if file is valid
    QC_PASA_STANDARD(
        ch_pasa_scaffolds
            .filter { _sample_id, std, _sens -> std.exists() && std.size() > 0 }
            .map { sample_id, std, _sens -> tuple(sample_id, std) },
        "pasa_standard"
    )

    // Run QC on sensitive PASA scaffolds - only if file is valid
    QC_PASA_SENSITIVE(
        ch_pasa_scaffolds
            .filter { _sample_id, _std, sens -> sens.exists() && sens.size() > 0 }
            .map { sample_id, _std, sens -> tuple(sample_id, sens) },
        "pasa_sensitive"
    )

    //=========================================================================
    // Optional: Intelligent PASA scaffold selection based on QC results
    //=========================================================================
    if (params.pasa_filter) {
        // Use left outer join to handle cases where only one PASA scaffold exists
        ch_pasa_scaffolds
            .join(QC_PASA_STANDARD.out.json_reports, remainder: true)
            .join(QC_PASA_SENSITIVE.out.json_reports, remainder: true)
            .map { sample_id, std, sens, qc_std, qc_sens ->
                // Use unique placeholder names per sample to avoid Nextflow file staging collisions
                def qc_std_file = qc_std ?: file("NO_QC_FILE_${sample_id}_std")
                def qc_sens_file = qc_sens ?: file("NO_QC_FILE_${sample_id}_sens")
                tuple(sample_id, std, sens, qc_std_file, qc_sens_file)
            }
            .set { ch_pasa_filter_input }

        PASA_FILTER(ch_pasa_filter_input)
        ch_selected_scaffolds = PASA_FILTER.out.selected_scaffold
    } else {
        // Default to sensitive scaffolds without filtering
        ch_selected_scaffolds = ch_pasa_scaffolds.map { sample_id, _standard, sensitive ->
            tuple(sample_id, sensitive)
        }
    }

    //=========================================================================
    // STEP 9: Generate summary report
    // Collects key metrics from all pipeline stages into a single TSV
    // Compares SPAdes vs PASA QC metrics to show improvement/degradation
    // SUMMARY runs after PASA_FILTER to ensure all data is available
    //=========================================================================

    // Collect SPAdes QC JSON files (completeness/contamination before scaffolding)
    ch_spades_qc_files = QC_SPADES.out.checkm2_json
        .map { sample_id, json_file -> json_file }
        .collect()
        .ifEmpty(file('NO_SPADES_QC'))

    // Collect PASA QC JSON files (both standard and sensitive - summary will pick best per sample)
    // These channels are populated by QC_PASA_STANDARD and QC_PASA_SENSITIVE which run after PASA
    ch_pasa_qc_files = QC_PASA_STANDARD.out.checkm2_json
        .mix(QC_PASA_SENSITIVE.out.checkm2_json)
        .map { sample_id, json_file -> json_file }
        .collect()
        .ifEmpty(file('NO_PASA_QC'))

    // Collect Bakta JSON files for gene counts
    ch_bakta_files = BAKTA.out.json_report
        .map { sample_id, json_file -> json_file }
        .collect()
        .ifEmpty(file('NO_BAKTA'))

    // Collect selected reference CSV files
    ch_refs_files = ch_selected_references
        .map { sample_id, csv_file -> csv_file }
        .collect()
        .ifEmpty(file('NO_REFS'))

    // Collect taxonomy data as list of maps
    ch_taxonomy_data = ch_taxonomy_results
        .map { sample_id, genus, species, family ->
            [sample_id: sample_id, genus: genus, species: species, family: family]
        }
        .collect()
        .ifEmpty([])

    // Generate summary TSV comparing SPAdes vs PASA metrics
    // Dependencies: QC_SPADES, QC_PASA_STANDARD, QC_PASA_SENSITIVE, BAKTA, CLASSIFICATION
    SUMMARY(
        ch_spades_qc_files,
        ch_pasa_qc_files,
        ch_bakta_files,
        ch_refs_files,
        ch_taxonomy_data
    )
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