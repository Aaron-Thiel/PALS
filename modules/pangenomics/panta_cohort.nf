/*
 * Panta pangenome analysis for cohort (all samples together)
 * Performs pangenome analysis on all annotated genomes as a single cohort
 */

process PANTA_COHORT {
    tag "pangenomics_${sample_count}_samples"

    publishDir "${params.outdir}/pangenomics/panta", mode: 'copy'

    container 'aaronthiel/panta-nextflow:latest'

    input:
    tuple val(cohort_id), path(gff_files), val(sample_count)

    output:
    tuple val(cohort_id), path("panta_results"), emit: panta_dir
    tuple val(cohort_id), path("panta_results/gene_presence_absence.Rtab"), emit: rtab
    tuple val(cohort_id), path("panta_results/gene_annotation.csv"), emit: gene_annotation
    tuple val(cohort_id), path("panta_results/gene_position.csv"), emit: gene_position
    tuple val(cohort_id), path("panta_results/annotated_clusters.json"), emit: clusters
    tuple val(cohort_id), path("panta_results/samples.json"), emit: samples
    tuple val(cohort_id), path("panta_results/representative_seqs.fasta"), emit: representative_seqs, optional: true
    tuple val(cohort_id), val(sample_count), emit: sample_count

    script:
    """
    # Create a directory for all GFF files
    mkdir -p input_gffs

    # Copy all GFF files to input directory
    for gff in ${gff_files}; do
        cp -L "\$gff" input_gffs/ 2>/dev/null || true
    done

    # Count input files
    GFF_COUNT=\$(ls -1 input_gffs/*.gff3 2>/dev/null | wc -l)
    echo "Running Panta pangenome analysis on \$GFF_COUNT genomes"

    # Run panta
    panta main \\
        --gff input_gffs/*.gff3 \\
        --outdir panta_results \\
        --threads ${task.cpus} \\
        --identity ${params.panta_identity} \\
        --LD ${params.panta_length_diff} \\
        --AL ${params.panta_align_long} \\
        --AS ${params.panta_align_short} \\
        --evalue ${params.panta_evalue} \\
        --blast ${params.panta_blast_method} \\
        \$([ "${params.panta_split_paralogs}" == "false" ] && echo "--dont-split" || echo "") \\
        --core ${params.panta_core} \\
        --soft ${params.panta_soft} \\
        --shell ${params.panta_shell}

    echo "Panta pangenome analysis completed for ${sample_count} samples"
    """
}
