/*
 * Concatenate Aligned Genes into Supermatrix
 * Creates concatenated alignment and partition file for phylogenetic analysis
 */

process CONCATENATE_ALIGNMENT {
    tag "$sample_id"
    label 'process_low'
    
    publishDir "${params.outdir}/phylogenetics/alignment", mode: 'copy'
    
    container 'python:3.9'
    
    input:
    tuple val(sample_id), path(aligned_genes_dir), val(sample_count)
    
    output:
    tuple val(sample_id), path("core_alignment.fasta"), emit: core_alignment
    tuple val(sample_id), path("alignment_stats.txt"), emit: alignment_stats
    tuple val(sample_id), path("partition.txt"), emit: partition_file
    tuple val(sample_id), val(sample_count), emit: sample_count
    
    script:
    def core_threshold = params.phylo_core_threshold ?: 0.95
    def mafft_method = params.mafft_method ?: 'auto'
    """
    #!/usr/bin/env python3
    
    import os
    from collections import defaultdict
    from pathlib import Path
    
    print("==============================================")
    print("Concatenating Alignments into Supermatrix")
    print("==============================================")
    print(f"Sample ID: ${sample_id}")
    print("")
    
    aligned_dir = Path("${aligned_genes_dir}")
    files = sorted([f for f in aligned_dir.iterdir() if f.suffix == '.aln'])
    
    print(f"Concatenating {len(files)} aligned genes")
    
    if len(files) == 0:
        print("WARNING: No aligned genes found")
        with open("core_alignment.fasta", 'w') as f:
            f.write(">no_core_genes_found\\nNNNNNNNNNNNNNNNNNNNN\\n")
        with open("alignment_stats.txt", 'w') as f:
            f.write("No core genes identified\\n")
        with open("partition.txt", 'w') as f:
            f.write("# Empty partition file\\n")
        exit(0)
    
    # Simple FASTA parser
    def parse_fasta(filepath):
        records = []
        current_id = None
        current_seq = []
        
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        records.append({'id': current_id, 'seq': ''.join(current_seq)})
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id:
                records.append({'id': current_id, 'seq': ''.join(current_seq)})
        return records
    
    # Collect sequences per sample
    sample_seqs = defaultdict(list)
    gene_info = []  # (name, length)
    
    for aln_file in files:
        gene = aln_file.stem
        records = parse_fasta(aln_file)
        
        if records:
            length = len(records[0]['seq'])
            gene_info.append((gene, length))
            
            for rec in records:
                # Extract sample name - handle various formats
                sample = rec['id']
                if '~~~' in sample:
                    sample = sample.split('~~~')[-1]
                sample_seqs[sample].append(rec['seq'])
    
    # Build concatenated alignment
    total_len = sum(l for _, l in gene_info)
    n_genes = len(gene_info)
    n_samples = len(sample_seqs)
    
    print(f"Total: {total_len} bp from {n_genes} genes across {n_samples} samples")
    
    # Write concatenated FASTA
    with open("core_alignment.fasta", 'w') as f:
        for sample, seqs in sorted(sample_seqs.items()):
            if len(seqs) == n_genes:
                concat = ''.join(seqs)
                f.write(f">{sample}\\n")
                for i in range(0, len(concat), 80):
                    f.write(concat[i:i+80] + "\\n")
            else:
                print(f"WARNING: {sample} has {len(seqs)} genes, expected {n_genes} - skipping")
    
    # Write partition file for IQ-TREE
    with open("partition.txt", 'w') as f:
        f.write("# Partition file for IQ-TREE\\n")
        f.write(f"# {n_genes} core genes, {total_len} total sites\\n")
        pos = 1
        for gene, length in gene_info:
            f.write(f"DNA, {gene} = {pos}-{pos + length - 1}\\n")
            pos += length
    
    # Write alignment statistics
    with open("alignment_stats.txt", 'w') as f:
        f.write("Core Genome Alignment Statistics\\n")
        f.write("=" * 50 + "\\n\\n")
        f.write(f"Core genes aligned: {n_genes}\\n")
        f.write(f"Samples in alignment: {n_samples}\\n")
        f.write(f"Total alignment length: {total_len} bp\\n")
        f.write(f"Average gene length: {total_len / n_genes if n_genes else 0:.1f} bp\\n")
        f.write(f"MAFFT method: ${mafft_method}\\n")
        f.write(f"Core threshold: ${core_threshold}\\n")
    
    print("Concatenation complete")
    print("==============================================")
    """
}