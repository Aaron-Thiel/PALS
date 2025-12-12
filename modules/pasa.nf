/*
 * PASA scaffolding process (Nextflow-friendly)
 * Scaffolds genome assemblies using pangenome synteny information.
 */

process PASA {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/pasa", mode: 'copy'

    container 'aaronthiel/pasa:latest'

    input:
    tuple val(sample_id),
          path(gene_annotation),      // from PANTA
          path(gene_position),        // from PANTA
          path(clusters),             // from PANTA
          path(samples_json),         // from PANTA
          path(contigs),              // from SPADES
          path(graph),                // from SPADES
          path(paths)                 // from SPADES

    output:
    tuple val(sample_id), path("scaffolds/pasa_contigs.fasta"), emit: standard_scaffold, optional: true
    tuple val(sample_id), path("scaffolds/pasa_sensitive_contigs.fasta"), emit: sensitive_scaffold, optional: true
    tuple val(sample_id), path("pasa_summary.json"), emit: summary, optional: true

    script:
    """
    # Create directory structure expected by PASA
    mkdir -p panta_data spades_data

    # Stage PANTA files
    cp ${gene_annotation} panta_data/
    cp ${gene_position} panta_data/
    cp ${clusters} panta_data/
    cp ${samples_json} panta_data/

    # Stage SPADES files
    cp ${contigs} spades_data/
    cp ${graph} spades_data/
    cp ${paths} spades_data/

    pasa \\
      --panta-dir panta_data \\
      --spades-dir spades_data \\
      --sample-name ${sample_id} \\
      --output-dir .
    """
}
