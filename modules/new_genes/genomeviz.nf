/*
 * GenomeViz Visualization
 *
 * Compares each internal genome (PASA scaffold and SPAdes contigs) against
 * its best reference genome using the aaronthiel/genomeviz Docker image.
 *
 * Supports different modes via the mode input parameter:
 *   - "all":        contig vs ref + scaffold vs ref + scaffold vs contig comparison
 *   - "comparison": scaffold vs contig only (produces gene_comparison_report.csv)
 *   - "scaffold":   scaffold vs reference only (produces alignment PNGs)
 *   - "contig":     contig vs reference only
 *
 * Uses stageAs to rename input files to match genomeviz's auto-detection
 * naming convention (reference.fna, scaffolds.fna, contigs.fna, etc.),
 * enabling the clean --input . invocation mode.
 */

process GENOME_VISUALIZATION {
    tag "$sample_id"
    label 'process_medium'

    container 'aaronthiel/genomeviz:latest'

    input:
    tuple val(sample_id),
          val(mode),
          path(scaffold_fna, stageAs: 'scaffolds.fna'),
          path(scaffold_gff, stageAs: 'scaffolds.gff3'),
          path(scaffold_faa, stageAs: 'scaffolds.faa'),
          path(ref_fna, stageAs: 'reference.fna'),
          path(ref_gff, stageAs: 'reference.gff3'),
          path(contig_fna, stageAs: 'contigs.fna'),
          path(contig_gff, stageAs: 'contigs.gff3'),
          path(contig_faa, stageAs: 'contigs.faa')

    output:
    tuple val(sample_id), path("${sample_id}"), emit: results

    script:
    """
    echo "=============================================="
    echo "GenomeViz: ${sample_id} (mode: ${mode})"
    echo "=============================================="
    echo "Scaffold: ${scaffold_fna}"
    echo "Scaffold GFF: ${scaffold_gff}"
    echo "Scaffold FAA: ${scaffold_faa}"
    echo "Reference: ${ref_fna}"
    echo "Reference GFF: ${ref_gff}"
    echo "Contigs: ${contig_fna}"
    echo "Contig GFF: ${contig_gff}"
    echo "Contig FAA: ${contig_faa}"
    echo ""

    python3 /app/genomeViz.py \\
        --input . \\
        --output ${sample_id} \\
        --mode ${mode} \\
        --no-interactive \\
        --no-gene-alignments \\
        --force

    echo "GenomeViz completed for ${sample_id} (mode: ${mode})"
    """
}
