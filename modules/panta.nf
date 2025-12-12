/*
 * Panta pangenome analysis process
 * Performs pangenome analysis on annotated genomes
 */

process PANTA {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/panta", mode: 'copy'

    container 'aaronthiel/panta-nextflow:latest'

    input:
    tuple val(sample_id), path(sample_gff), path(reference_gffs)

    output:
    tuple val(sample_id), path("gene_annotation.csv"), emit: gene_annotation
    tuple val(sample_id), path("gene_position.csv"), emit: gene_position
    tuple val(sample_id), path("annotated_clusters.json"), emit: clusters
    tuple val(sample_id), path("samples.json"), emit: samples
    
    script:
    """
    # Create a directory for all GFF files
    mkdir -p input_gffs
    
    # Copy all GFF files to input directory
    cp -L "${sample_gff}" input_gffs/ 2>/dev/null || true
    cp -L ${reference_gffs} input_gffs/ 2>/dev/null || true
    
    # Run panta - now works as regular command in Nextflow-compatible image
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
    
    # Move panta_results content to current directory
    if [ -d panta_results ]; then
        mv panta_results/* . 2>/dev/null || true
        rmdir panta_results 2>/dev/null || true
    fi
    """
}
