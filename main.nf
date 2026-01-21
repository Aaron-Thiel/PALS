#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    BGC-link Pipeline Orchestrator
========================================================================================
    Main entry point for running BGC-link pipelines.

    Available pipelines:
      - assembly:  Genome assembly, annotation, and scaffolding (assembly.nf)
      - analysis:  Pangenomics and phylogenetics (analysis.nf)

    Usage:
      nextflow run main.nf --pipeline assembly --input samplesheet.csv
      nextflow run main.nf --pipeline analysis
      nextflow run main.nf --pipeline all --input samplesheet.csv

    Or run pipelines directly:
      nextflow run assembly.nf -c assembly.config --input samplesheet.csv
      nextflow run analysis.nf -c analysis.config

    Author: Created for BGC-link project
========================================================================================
*/

// =========================================================================
// Parameters
// =========================================================================

params.pipeline = null  // 'assembly', 'analysis', or 'all'
params.help = false

// =========================================================================
// Help message
// =========================================================================

def helpMessage() {
    log.info """
    ===================================
    BGC-link Pipeline Orchestrator
    ===================================

    Usage:
      nextflow run main.nf --pipeline <pipeline_name> [options]

    Available pipelines:
      assembly    Run the assembly pipeline (reads → scaffolds)
                  Requires: --input samplesheet.csv
                  Config: assembly.config

      analysis    Run the analysis pipeline (annotations → pangenome/phylogeny)
                  Config: analysis.config

      all         Run assembly followed by analysis
                  Requires: --input samplesheet.csv

    Options:
      --help      Show this help message

    Examples:
      # Run assembly pipeline
      nextflow run main.nf --pipeline assembly --input samples.csv

      # Run analysis pipeline on existing results
      nextflow run main.nf --pipeline analysis

      # Run full workflow (assembly + analysis)
      nextflow run main.nf --pipeline all --input samples.csv

    Or run pipelines directly:
      nextflow run assembly.nf -c assembly.config --input samples.csv
      nextflow run analysis.nf -c analysis.config

    ===================================
    """.stripIndent()
}

// =========================================================================
// Main Workflow
// =========================================================================

workflow {
    if (params.help || !params.pipeline) {
        helpMessage()
        if (!params.help) {
            error "Please specify a pipeline with --pipeline (assembly, analysis, or all)"
        }
        return
    }

    log.info """
    ===================================
    BGC-link Pipeline Orchestrator
    ===================================
    Pipeline: ${params.pipeline}
    ===================================
    """.stripIndent()

    def pipeline = params.pipeline.toLowerCase()

    if (pipeline == 'assembly') {
        log.info "Starting assembly pipeline..."
        log.info "Run: nextflow run assembly.nf -c assembly.config --input ${params.input ?: '<samplesheet>'}"

        if (!params.input) {
            error "Assembly pipeline requires --input samplesheet.csv"
        }

    } else if (pipeline == 'analysis') {
        log.info "Starting analysis pipeline..."
        log.info "Run: nextflow run analysis.nf -c analysis.config"

    } else if (pipeline == 'all') {
        log.info "Starting full workflow (assembly → analysis)..."
        log.info "Step 1: nextflow run assembly.nf -c assembly.config --input ${params.input ?: '<samplesheet>'}"
        log.info "Step 2: nextflow run analysis.nf -c analysis.config"

        if (!params.input) {
            error "Full workflow requires --input samplesheet.csv"
        }

    } else {
        error "Unknown pipeline: ${params.pipeline}. Choose from: assembly, analysis, all"
    }

    log.info """

    Note: This orchestrator shows the commands to run.
    For actual execution, run the pipelines directly:

      nextflow run assembly.nf -c assembly.config [options]
      nextflow run analysis.nf -c analysis.config [options]

    ===================================
    """.stripIndent()
}
