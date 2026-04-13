/*
 * PREPROCESSING WORKFLOW
 *
 * Complete read preprocessing pipeline with quality assessment:
 * 1. FASTQC (raw): Initial quality assessment
 * 2. FASTP: Read trimming and quality control
 * 3. KRAKEN2: Taxonomic classification
 * 4. KRAKENTOOLS: Filter for target family reads (optional, controlled by params.filter_reads)
 * 5. FASTQC (processed): Final quality assessment
 */

include { FASTQC } from './preprocessing/fastqc.nf'
include { FASTQC as FASTQC_PROCESSED } from './preprocessing/fastqc.nf'
include { FASTP } from './preprocessing/fastp.nf'
include { KRAKEN2 } from './preprocessing/kraken2.nf'
include { KRAKENTOOLS } from './preprocessing/krakentools.nf'

workflow PREPROCESSING {
    take:
    reads_ch  // tuple val(sample_id), path(reads)

    main:

    log.info """
    ════════════════════════════════════════════════════════════════
    PREPROCESSING WORKFLOW
    ════════════════════════════════════════════════════════════════
    Kraken2 database: ${params.kraken2_db}
    Target family: ${params.filter_target_family}
    Filter reads: ${params.filter_reads}
    Output directory: ${params.outdir}
    ════════════════════════════════════════════════════════════════
    """

    //=========================================================================
    // STEP 1: Initial quality assessment with FastQC (raw reads)
    //=========================================================================
    FASTQC(reads_ch, "raw")

    //=========================================================================
    // STEP 2: Read quality control and trimming with fastp
    //=========================================================================
    FASTP(reads_ch)

    //=========================================================================
    // STEP 3: Taxonomic classification with Kraken2
    //=========================================================================
    KRAKEN2(FASTP.out.trimmed_reads)

    //=========================================================================
    // STEP 4: Filter reads for target family and unclassified reads (optional)
    //=========================================================================
    if (params.filter_reads) {
        // Prepare KrakenTools input
        kraken_filter_input = KRAKEN2.out.kraken_output
            .join(KRAKEN2.out.kraken_report)
            .join(FASTP.out.trimmed_reads)
            .map { row ->
                def sample_id = row[0]
                def kraken_output = row[1]
                def kraken_report = row[2]
                def trimmed_reads = row[3]
                tuple(sample_id, trimmed_reads, kraken_output, kraken_report)
            }

        KRAKENTOOLS(kraken_filter_input)

        ch_processed_reads = KRAKENTOOLS.out.filtered_reads
        ch_krakentools_summary = KRAKENTOOLS.out.summary
    } else {
        // Skip filtering - use fastp trimmed reads directly
        log.info "Skipping KrakenTools filtering (filter_reads=false)"
        ch_processed_reads = FASTP.out.trimmed_reads
        ch_krakentools_summary = channel.empty()
    }

    //=========================================================================
    // STEP 5: Final quality assessment with FastQC (processed reads)
    //=========================================================================
    FASTQC_PROCESSED(ch_processed_reads, "processed")

    //=========================================================================
    // OUTPUT CHANNELS
    //=========================================================================

    emit:
    // Primary output - filtered and quality-controlled reads
    filtered_reads = ch_processed_reads

    // Quality assessment reports
    raw_fastqc_html = FASTQC.out.html_reports.filter { _sample_id, reports ->
        reports.any { report -> report.name.contains('raw') }
    }
    processed_fastqc_html = FASTQC_PROCESSED.out.html_reports.filter { _sample_id, reports ->
        reports.any { report -> report.name.contains('processed') }
    }

    // Individual module outputs for reporting
    fastp_reports = FASTP.out.json_report
    kraken_reports = KRAKEN2.out.kraken_report
    krakentools_summary = ch_krakentools_summary

    // Summary outputs for MultiQC integration
    raw_fastqc_summary = FASTQC.out.summary.filter { _sample_id, summary ->
        def content = summary.text
        content.contains('Stage: raw')
    }
    processed_fastqc_summary = FASTQC_PROCESSED.out.summary.filter { _sample_id, summary ->
        def content = summary.text
        content.contains('Stage: processed')
    }
}