/*
 * Lost Genes Analysis
 *
 * Identifies genes present in SPAdes contigs (pre-PASA) but absent from
 * PASA scaffolds, classifies their completeness (complete vs truncated
 * at contig boundaries), then maps them to COG categories via the SPAdes
 * PANTA pangenome and EggNOG annotations.
 *
 * Produces three subsets of outputs (all / high_completeness / low_completeness)
 * each containing TSV summaries, COG distribution plots, and per-species/genus breakdowns.
 *
 * Inputs are staged into separate subdirectories to avoid filename collisions
 * between contig and scaffold FAA files.
 */

process LOST_GENES {
    tag "lost_genes_analysis"
    label 'process_medium'

    publishDir "${params.outdir}/lost_genes", mode: 'copy'

    container 'aaronthiel/lost_genes:latest'

    input:
    path(contig_faa_files, stageAs: 'contig_faa/*')
    path(contig_gff_files, stageAs: 'contig_gff/*')
    path(scaffold_faa_files, stageAs: 'scaffold_faa/*')
    path(pipeline_summary)
    path(panta_csv)
    path(eggnog_annotations)
    path(rtab_file)

    output:
    path("results/all"),               emit: all_results
    path("results/high_completeness"), emit: high_results
    path("results/low_completeness"),  emit: low_results

    script:
    def min_scaffolding = params.lost_genes_min_scaffolding ?: '0.9'
    """
    lost_genes_analysis.py \
        --contig-faa-dir contig_faa \
        --contig-gff-dir contig_gff \
        --scaffold-faa-dir scaffold_faa \
        --pipeline-summary ${pipeline_summary} \
        --panta-csv ${panta_csv} \
        --eggnog ${eggnog_annotations} \
        --rtab ${rtab_file} \
        --outdir results \
        --min-scaffolding ${min_scaffolding}
    """
}
