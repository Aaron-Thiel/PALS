/*
 * QUAST assembly quality assessment
 * Provides detailed assembly statistics and comparison
 */

process QUAST {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/qc_${stage_name}/quast", mode: 'copy'

    container 'staphb/quast:latest'

    input:
    tuple val(sample_id), path(assembly)
    path reference_genome, stageAs: 'reference.fasta'
    val stage_name
    
    output:
    tuple val(sample_id), path("quast_results"), emit: results
    tuple val(sample_id), path("quast_summary.txt"), emit: summary
    tuple val(sample_id), path("*.tsv"), emit: tsv_reports
    
    script:
    def reference_args = reference_genome.name != 'null' && reference_genome.name != 'empty_reference.fasta' && reference_genome.size() > 0 ? "--reference reference.fasta" : ""
    
    """
    echo "Running QUAST assembly assessment on ${assembly}..."
    
    # Run QUAST on single assembly
    quast.py \\
        ${assembly} \\
        --output-dir quast_results \\
        --threads ${task.cpus} \\
        --min-contig 500 \\
        ${reference_args} \\
        --no-plots \\
        --no-html \\
        --no-icarus
    
    # Extract key metrics and create summary
    if [ -f "quast_results/transposed_report.tsv" ]; then
        cp quast_results/transposed_report.tsv quast_transposed_report.tsv
        
        # Extract key metrics
        contigs=\$(grep "# contigs" quast_results/transposed_report.tsv | cut -f2)
        total_length=\$(grep "Total length" quast_results/transposed_report.tsv | cut -f2)
        n50=\$(grep "N50" quast_results/transposed_report.tsv | cut -f2)
        l50=\$(grep "L50" quast_results/transposed_report.tsv | cut -f2)
        largest_contig=\$(grep "Largest contig" quast_results/transposed_report.tsv | cut -f2)
        
        # Check for misassemblies if reference was provided
        misassemblies=""
        if [ -f "reference.fasta" ] && [ "${reference_genome.name}" != "null" ] && [ "${reference_genome.name}" != "empty_reference.fasta" ] && [ -s "reference.fasta" ]; then
            misassemblies=\$(grep "# misassemblies" quast_results/transposed_report.tsv | cut -f2 2>/dev/null || echo "N/A")
            misassembled_contigs=\$(grep "# misassembled contigs" quast_results/transposed_report.tsv | cut -f2 2>/dev/null || echo "N/A")
        fi
        
        # Create summary report
        cat > quast_summary.txt << EOF
QUAST Assembly Quality Assessment for ${sample_id}
=================================================

Assembly: ${assembly}
Number of contigs: \${contigs}
Total length: \${total_length}
N50: \${n50}
L50: \${l50}
Largest contig: \${largest_contig}
EOF

        # Add misassembly information if reference was used
        if [ -f "reference.fasta" ] && [ "${reference_genome.name}" != "null" ] && [ "${reference_genome.name}" != "empty_reference.fasta" ] && [ -s "reference.fasta" ]; then
            cat >> quast_summary.txt << EOF

Reference-based Quality Assessment:
Number of misassemblies: \${misassemblies}
Number of misassembled contigs: \${misassembled_contigs}
EOF
        fi
        
        echo "" >> quast_summary.txt
        echo "Quality Interpretation:" >> quast_summary.txt
        echo "- Higher N50 values indicate better contiguity" >> quast_summary.txt
        echo "- Lower L50 values indicate better contiguity" >> quast_summary.txt
        echo "- Fewer contigs generally indicate better assembly" >> quast_summary.txt
        if [ -f "reference.fasta" ] && [ "${reference_genome.name}" != "null" ] && [ "${reference_genome.name}" != "empty_reference.fasta" ] && [ -s "reference.fasta" ]; then
            echo "- Lower misassembly counts indicate higher accuracy" >> quast_summary.txt
        fi
        
    else
        echo "QUAST analysis failed - no results generated" > quast_summary.txt
    fi
    
    echo "QUAST analysis completed for ${sample_id}"
    """
}