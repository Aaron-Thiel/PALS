/*
 * SPAdes assembly process
 * Assembles Illumina paired-end reads into contigs
 */

process SPADES {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/spades", mode: 'copy'

    container 'staphb/spades:latest'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("contigs.fasta"), emit: contigs
    tuple val(sample_id), path("assembly_graph.fastg"), emit: graph
    tuple val(sample_id), path("contigs.paths"), emit: paths, optional: true
    tuple val(sample_id), path("spades.log"), emit: log
    
    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    spades.py \\
        --isolate \\
        --cov-cutoff auto \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o . \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()}
    """
}
