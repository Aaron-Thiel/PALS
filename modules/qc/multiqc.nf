/*
 * MultiQC aggregation process
 * Creates comprehensive quality control report
 * Aggregates results from all QC tools
 */

process MULTIQC {
    tag "all_samples"
    
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    container 'staphb/multiqc:latest'
    
    input:
    path(qc_files)
    
    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"), emit: data
    path("pipeline_summary.txt"), emit: summary
    
    script:
    """
    # Create MultiQC config file
    cat > multiqc_config.yaml << EOF
title: "BGC-link Pipeline QC Report"
subtitle: "SPAdes-Bakta-Panta-PASA-RagTag Pipeline Quality Control"
intro_text: "This report aggregates quality control metrics from the bacterial genome assembly and annotation pipeline."

report_header_info:
  - Contact E-mail: 'your-email@example.com'
  - Pipeline: 'BGC-link'
  - Version: '1.0'

module_order:
  - busco
  - checkm2
  - quast
  - kraken
  - fastqc
  - custom_content

custom_data:
  assembly_comparison:
    file_format: 'json'
    section_name: 'Assembly Comparison'
    description: 'Comparison between PASA standard and sensitive scaffolds'
    plot_type: 'bargraph'
    pconfig:
      id: 'assembly_comparison'
      title: 'Assembly Quality Comparison'
      ylab: 'Quality Score'

sp:
  busco:
    fn: '*busco*summary*.json'
  checkm2:
    fn: '*checkm2*summary*.json'
EOF

    # Create pipeline summary
    cat > pipeline_summary.txt << EOF
BGC-link Pipeline Quality Control Summary
=========================================

Pipeline Steps Completed:
1. ✓ SPAdes - De novo genome assembly
2. ✓ Bakta - Bacterial genome annotation  
3. ✓ Panta - Pangenome analysis
4. ✓ PASA - Synteny-based scaffolding
5. ✓ MASH - Reference similarity analysis
6. ✓ RagTag - Reference-guided scaffolding (conditional)
7. ✓ BUSCO - Completeness assessment
8. ✓ CheckM2 - Quality assessment

Quality Control Tools:
- BUSCO: Genome completeness using lactobacillales_odb10 lineage
- CheckM2: Contamination and completeness assessment
- QUAST: Assembly statistics and comparison
- Kraken2: Contamination detection (if included)

Output Files:
- Assembly: SPAdes contigs and scaffolds
- Annotation: Bakta GFF3, GenBank, and protein files
- Pangenome: Panta cluster analysis
- Scaffolding: PASA improved scaffolds (standard and sensitive)
- Reference scaffolding: RagTag scaffolds (when reference similarity >95%)
- Quality reports: BUSCO and CheckM2 assessments

Recommendations:
Check the MultiQC report sections for:
1. Assembly quality metrics
2. Completeness scores
3. Contamination levels
4. Best scaffold selection (standard vs sensitive PASA)
EOF

    # Run MultiQC
    multiqc \\
        --config multiqc_config.yaml \\
        --title "BGC-link Pipeline QC Report" \\
        --filename multiqc_report.html \\
        --force \\
        --verbose \\
        .
    
    # Add timestamp to summary
    echo "" >> pipeline_summary.txt
    echo "Report generated: \$(date)" >> pipeline_summary.txt
    echo "MultiQC version: \$(multiqc --version | head -1)" >> pipeline_summary.txt
    
    echo "MultiQC report generated successfully"
    """
}