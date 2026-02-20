/*
 * BiG-SCAPE 2 — BGC similarity network clustering
 *
 * Runs BiG-SCAPE 2 in cluster mode on a directory of antiSMASH region .gbk files.
 * Compares against MiBIG database for known BGC family identification.
 * Produces similarity networks, GCF assignments, and an interactive HTML viewer.
 */

process BIGSCAPE {
    tag "$group_id"
    label 'process_medium'

    publishDir "${params.outdir}/pangenomics/bigscape/${group_id}", mode: 'copy'

    container 'bigscape:latest'

    input:
    tuple val(group_id), path(bgc_dir)
    path(pfam_db)

    output:
    tuple val(group_id), path("bigscape_output"), emit: results
    tuple val(group_id), path("bigscape_output/output_files"), emit: output_files, optional: true

    script:
    def cutoffs = params.bigscape_cutoffs ?: '0.3'
    """
    echo "=============================================="
    echo "BiG-SCAPE 2 — ${group_id}"
    echo "=============================================="
    echo "Input directory: ${bgc_dir}"
    echo "Pfam database: ${pfam_db}"
    echo "GCF cutoffs: ${cutoffs}"
    echo ""

    # Count input BGC files
    N_BGCS=\$(find -L ${bgc_dir} -name '*.gbk' | wc -l)
    echo "Total .gbk files found: \$N_BGCS"

    if [ "\$N_BGCS" -eq 0 ]; then
        echo "WARNING: No .gbk files found. Creating empty output."
        mkdir -p bigscape_output/output_files
        exit 0
    fi

    bigscape.py cluster \
        -i ${bgc_dir} \
        -o bigscape_output \
        -p ${pfam_db} \
        --mibig-version 3.1 \
        --gcf-cutoffs ${cutoffs} \
        --include-singletons \
        --mix \
        -c ${task.cpus} \
        -v

    echo ""
    echo "BiG-SCAPE completed for ${group_id}!"
    """
}
