/*
 * FASTQC - Read Quality Assessment
 * 
 * Generates quality control reports for sequencing reads using FastQC
 */

process FASTQC {
    tag "${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}/preprocessing/fastqc", mode: 'copy'

    container 'staphb/fastqc:latest'

    input:
    tuple val(sample_id), path(reads)
    val(suffix) // e.g., "raw" or "trimmed"

    output:
    tuple val(sample_id), path("*.html"), emit: html_reports
    tuple val(sample_id), path("*.zip"), emit: zip_reports
    tuple val(sample_id), path("fastqc_summary.txt"), emit: summary

    script:
    def input_files = reads instanceof List ? reads.join(' ') : reads
    """
    echo "Running FastQC on ${suffix} reads for ${sample_id}..."
    
    # Run FastQC
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        --format fastq \\
        ${input_files}
    
    # Create summary report
    cat > fastqc_summary.txt << EOF
Sample ID: ${sample_id}
Stage: ${suffix}
Input Files: ${input_files}
Timestamp: \$(date)
EOF
    
    echo "FastQC analysis completed for ${sample_id} (${suffix} reads)"
    """
}