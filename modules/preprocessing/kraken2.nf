/*
 * Kraken2 taxonomic classification of reads
 * Classifies reads against database for subsequent filtering
 */

process KRAKEN2 {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/preprocessing/kraken2", mode: 'copy'
    
    container 'staphb/kraken2:latest'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_classified.fastq"), emit: classified_reads
    tuple val(sample_id), path("${sample_id}_unclassified.fastq"), emit: unclassified_reads
    tuple val(sample_id), path("${sample_id}_kraken2_output.txt"), emit: kraken_output
    tuple val(sample_id), path("${sample_id}_kraken2_report.txt"), emit: kraken_report
    tuple val(sample_id), path("${sample_id}_classification_summary.txt"), emit: summary
    
    script:
    def read1 = reads[0]
    def read2 = reads[1]
    
    """
    # Run Kraken2 classification on paired-end reads
    kraken2 \\
        --db /databases/kraken2 \\
        --threads ${task.cpus} \\
        --paired \\
        --output ${sample_id}_kraken2_output.txt \\
        --report ${sample_id}_kraken2_report.txt \\
        --classified-out ${sample_id}_classified#.fastq \\
        --unclassified-out ${sample_id}_unclassified#.fastq \\
        ${read1} ${read2}
    
    # Combine paired files into single files for easier processing
    cat ${sample_id}_classified_1.fastq ${sample_id}_classified_2.fastq > ${sample_id}_classified.fastq
    cat ${sample_id}_unclassified_1.fastq ${sample_id}_unclassified_2.fastq > ${sample_id}_unclassified.fastq
    
    # Clean up intermediate files
    rm -f ${sample_id}_classified_1.fastq ${sample_id}_classified_2.fastq
    rm -f ${sample_id}_unclassified_1.fastq ${sample_id}_unclassified_2.fastq
    
    # Extract classification statistics
    total_reads=\$(wc -l < ${sample_id}_kraken2_output.txt)
    classified_reads=\$(grep -v "^U" ${sample_id}_kraken2_output.txt | wc -l)
    unclassified_reads=\$(grep "^U" ${sample_id}_kraken2_output.txt | wc -l)
    
    # Get Lactobacillaceae statistics
    lactobacillaceae_reads=\$(grep "Lactobacillaceae" ${sample_id}_kraken2_report.txt | head -1 | awk '{print \$3}' || echo "0")
    lactobacillaceae_percent=\$(grep "Lactobacillaceae" ${sample_id}_kraken2_report.txt | head -1 | awk '{print \$1}' || echo "0.00")
    
    # Create summary report
    cat > ${sample_id}_classification_summary.txt << EOF
Kraken2 Classification Summary for ${sample_id}
===============================================

Total reads processed: \$total_reads
Classified reads: \$classified_reads (\$(awk 'BEGIN {printf "%.2f", '\$classified_reads' * 100 / '\$total_reads'}')%)
Unclassified reads: \$unclassified_reads (\$(awk 'BEGIN {printf "%.2f", '\$unclassified_reads' * 100 / '\$total_reads'}')%)

Lactobacillaceae reads: \$lactobacillaceae_reads (\$lactobacillaceae_percent%)

Next step: Use KrakenTools to extract Lactobacillaceae and unclassified reads
for assembly to remove contamination.
EOF
    
    echo "Kraken2 classification completed for ${sample_id}"
    """
}