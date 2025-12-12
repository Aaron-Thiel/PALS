/*
 * EggNOG-mapper process for functional annotation (cohort-level)
 * Uses representative protein sequences from Panta pangenome clusters
 * Annotates with orthology assignments, GO terms, and KEGG pathways
 */

process EGGNOG {
    tag "$cohort_id"

    publishDir "${params.outdir}/pangenomics/eggnog", mode: 'copy'

    container 'nanozoo/eggnog-mapper:2.1.13--c16a7d2'

    input:
    tuple val(cohort_id), path(panta_dir)

    output:
    tuple val(cohort_id), path("${cohort_id}.emapper.annotations"), emit: annotations
    tuple val(cohort_id), path("${cohort_id}.emapper.hits"), emit: hits
    tuple val(cohort_id), path("${cohort_id}.emapper.seed_orthologs"), emit: seed_orthologs
    tuple val(cohort_id), path("eggnog_results"), emit: results_dir

    script:
    def cpu_opt = task.cpus ? "--cpu ${task.cpus}" : "--cpu 4"
    def taxon_opt = params.eggnog_taxon ? "-d ${params.eggnog_taxon}" : ""
    """
    # Create output directory
    mkdir -p eggnog_results

    # Find representative protein sequences from Panta output
    REP_PROTS=""
    if [ -f "${panta_dir}/representative_clusters_prot.fasta" ]; then
        REP_PROTS="${panta_dir}/representative_clusters_prot.fasta"
    elif [ -f "${panta_dir}/representative_seqs.fasta" ]; then
        REP_PROTS="${panta_dir}/representative_seqs.fasta"
    else
        # Search for protein fasta files
        REP_PROTS=\$(find ${panta_dir} -name "*prot*.fasta" -o -name "*protein*.fasta" | head -1)
    fi

    if [ -z "\$REP_PROTS" ] || [ ! -f "\$REP_PROTS" ]; then
        echo "ERROR: Could not find representative protein sequences in Panta output"
        echo "Contents of panta_dir:"
        ls -la ${panta_dir}/
        exit 1
    fi

    echo "Using protein sequences: \$REP_PROTS"
    echo "Number of sequences: \$(grep -c '^>' \$REP_PROTS)"
    echo "Search mode: ${params.eggnog_mode}"
    echo "Taxon database: ${params.eggnog_taxon ?: 'default'}"

    # Run eggnog-mapper on representative cluster proteins
    # Using HMMER mode with Lactobacillaceae-specific database for faster, focused annotation
    emapper.py \\
        -i "\$REP_PROTS" \\
        --output ${cohort_id} \\
        --output_dir eggnog_results \\
        ${cpu_opt} \\
        -m ${params.eggnog_mode} \\
        ${taxon_opt} \\
        --data_dir ${params.eggnog_db} \\
        --override

    # Move main output files to working directory for easier access
    mv eggnog_results/${cohort_id}.emapper.annotations .
    mv eggnog_results/${cohort_id}.emapper.hits .
    mv eggnog_results/${cohort_id}.emapper.seed_orthologs .

    echo "EggNOG-mapper annotation completed for cohort ${cohort_id}"
    """
}
