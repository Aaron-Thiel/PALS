/*
 * fastp read quality control and trimming
 * Performs adapter removal, quality trimming, and generates QC reports
 */

process FASTP {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/preprocessing/fastp", mode: 'copy'
    
    container 'staphb/fastp:latest'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.html"), emit: html_report
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json_report
    tuple val(sample_id), path("${sample_id}_fastp_summary.txt"), emit: summary
    
    script:
    def read1 = reads[0]
    def read2 = reads[1]
    
    """
    fastp \\
        --in1 ${read1} \\
        --in2 ${read2} \\
        --out1 ${sample_id}_R1_trimmed.fastq.gz \\
        --out2 ${sample_id}_R2_trimmed.fastq.gz \\
        --html ${sample_id}_fastp.html \\
        --json ${sample_id}_fastp.json \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --correction \\
        --cut_front \\
        --cut_tail \\
        --cut_window_size 4 \\
        --cut_mean_quality 20 \\
        --qualified_quality_phred 20 \\
        --unqualified_percent_limit 10 \\
        --length_required 50
    
    # Create summary report
    cat > ${sample_id}_fastp_summary.txt << EOF
fastp Quality Control Summary for ${sample_id}
==============================================

Input files:
  Read 1: ${read1}
  Read 2: ${read2}

Output files:
  Trimmed Read 1: ${sample_id}_R1_trimmed.fastq.gz
  Trimmed Read 2: ${sample_id}_R2_trimmed.fastq.gz

Settings:
  - Adapter detection: Enabled for paired-end
  - Quality trimming: Q20 sliding window (size 4)
  - Length filtering: Minimum 50bp
  - Error correction: Enabled
  
See ${sample_id}_fastp.html for detailed QC metrics and plots.
EOF
    
    echo "fastp processing completed for ${sample_id}"
    """
}