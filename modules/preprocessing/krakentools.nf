/*
 * KrakenTools filtering to extract specific taxa
 * Extracts Lactobacillaceae and unclassified reads for clean assembly
 */

process KRAKENTOOLS {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/preprocessing/krakentools", mode: 'copy'
    
    container 'staphb/krakentools:latest'
    
    input:
    tuple val(sample_id), path(reads), path(kraken_output), path(kraken_report)
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered_R{1,2}.fastq.gz"), emit: filtered_reads
    tuple val(sample_id), path("${sample_id}_filtering_summary.txt"), emit: summary
    tuple val(sample_id), path("${sample_id}_extracted_taxa.txt"), emit: taxa_list
    
    script:
    def read1 = reads[0]
    def read2 = reads[1]
    
    """
    # Get Lactobacillaceae taxid from the report
    lactobacillaceae_taxid=\$(grep -E "\\sF\\s.*Lactobacillaceae" ${kraken_report} | head -1 | awk '{print \$5}' || echo "33958")
    
    # Extract read IDs for Lactobacillaceae and unclassified reads
    # Unclassified reads have taxid 0
    # Lactobacillaceae and its descendants
    python3 /KrakenTools/extract_kraken_reads.py \\
        -k ${kraken_output} \\
        -s1 ${read1} \\
        -s2 ${read2} \\
        -o ${sample_id}_filtered_R1.fastq \\
        -o2 ${sample_id}_filtered_R2.fastq \\
        -t \$lactobacillaceae_taxid 0 \\
        -r ${kraken_report} \\
        --include-children \\
        --fastq-output
    
    # Compress the filtered reads
    gzip ${sample_id}_filtered_R1.fastq ${sample_id}_filtered_R2.fastq
    
    # Count reads before and after filtering
    original_r1=\$(( \$(zcat ${read1} | wc -l) / 4 ))
    original_r2=\$(( \$(zcat ${read2} | wc -l) / 4 ))
    filtered_r1=\$(( \$(zcat ${sample_id}_filtered_R1.fastq.gz | wc -l) / 4 ))
    filtered_r2=\$(( \$(zcat ${sample_id}_filtered_R2.fastq.gz | wc -l) / 4 ))
    
    # Calculate retention percentages using awk (more portable than bc)
    retention_r1=\$(awk 'BEGIN {if('\$original_r1' > 0) printf "%.2f", '\$filtered_r1' * 100 / '\$original_r1'; else print "0.00"}')
    retention_r2=\$(awk 'BEGIN {if('\$original_r2' > 0) printf "%.2f", '\$filtered_r2' * 100 / '\$original_r2'; else print "0.00"}')
    
    # Create taxa list for reference
    cat > ${sample_id}_extracted_taxa.txt << EOF
Extracted Taxa for ${sample_id}
==============================

Target taxa:
- Unclassified reads (taxid: 0)
- Lactobacillaceae (taxid: \$lactobacillaceae_taxid) and descendants

This filtering strategy retains:
1. Unclassified reads (potentially novel Lactobacillaceae sequences)
2. All Lactobacillaceae family members and their descendants
3. Removes contaminating taxa from other bacterial families

EOF
    
    # Create summary report
    cat > ${sample_id}_filtering_summary.txt << EOF
KrakenTools Filtering Summary for ${sample_id}
==============================================

Extraction criteria:
- Unclassified reads (taxid: 0)
- Lactobacillaceae (taxid: \$lactobacillaceae_taxid) and descendants

Read counts:
                Original    Filtered    Retention
Read 1:         \$original_r1       \$filtered_r1      \$retention_r1%
Read 2:         \$original_r2       \$filtered_r2      \$retention_r2%

Output files:
- ${sample_id}_filtered_R1.fastq.gz
- ${sample_id}_filtered_R2.fastq.gz

These filtered reads will be used for assembly, removing contamination
from non-Lactobacillaceae organisms while retaining target sequences
and potentially novel unclassified sequences.
EOF
    
    echo "KrakenTools filtering completed for ${sample_id}"
    echo "Retained \$retention_r1% of R1 reads and \$retention_r2% of R2 reads"
    """
}