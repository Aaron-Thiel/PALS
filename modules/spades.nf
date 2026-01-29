/*
 * SPAdes assembly process
 * Assembles Illumina paired-end reads into contigs
 *
 * Configurable parameters (set in assembly.config):
 *   spades_mode       - Assembly mode: '--isolate', '--careful', '--meta', or ''
 *   spades_cov_cutoff - Coverage cutoff: 'auto', 'off', or a number
 *   spades_kmers      - K-mer sizes: '' for auto, or '21,33,55,77'
 *   spades_extra_args - Additional arguments
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
    def mode_arg = params.spades_mode ?: ''
    def cov_arg = params.spades_cov_cutoff ? "--cov-cutoff ${params.spades_cov_cutoff}" : ''
    def kmer_arg = params.spades_kmers ? "-k ${params.spades_kmers}" : ''
    def extra_args = params.spades_extra_args ?: ''
    """
    spades.py \\
        ${mode_arg} \\
        ${cov_arg} \\
        ${kmer_arg} \\
        ${extra_args} \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o . \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()}
    """
}
