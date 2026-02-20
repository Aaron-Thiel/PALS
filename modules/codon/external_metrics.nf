/*
 * External Codon Metrics
 *
 * Batch process per genus: computes per-gene codon metrics (ENC, CAI, GC3s, GC_content)
 * from pre-extracted FFN (CDS nucleotide) files.
 *
 * Runs once per genus. Results are cached via storeDir.
 */

process EXTERNAL_CODON_METRICS {
    tag "$genus"
    label 'process_low'

    storeDir { "${params.outdir}/codon/external_metrics/${genus}" }

    conda 'conda-forge::pandas conda-forge::numpy'
    container null

    input:
    tuple val(genus), path(ffn_files), path(codon_reference)

    output:
    tuple val(genus), path("${genus}_external_codon_metrics.tsv"), emit: metrics

    script:
    """
    #!/usr/bin/env python3
    import json
    import numpy as np
    import pandas as pd
    from collections import defaultdict
    from pathlib import Path
    import glob

    # =========================================================================
    # Codon table and metric functions (shared with qc_analysis.nf)
    # =========================================================================

    CODON_TABLE = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    AA_TO_CODONS = defaultdict(list)
    for codon, aa in CODON_TABLE.items():
        AA_TO_CODONS[aa].append(codon)

    SYNONYMOUS_CODONS = [c for c in CODON_TABLE if CODON_TABLE[c] not in ['M', 'W', '*']]

    def calculate_enc(codon_counts):
        F_values = {}
        for aa, codons in AA_TO_CODONS.items():
            if aa in ['M', 'W', '*']:
                continue
            n = len(codons)
            total = sum(codon_counts.get(c, 0) for c in codons)
            if total <= 1:
                continue
            p_sum = sum((codon_counts.get(c, 0) / total) ** 2 for c in codons)
            F = (total * p_sum - 1) / (total - 1)
            if n not in F_values:
                F_values[n] = []
            F_values[n].append(F)
        F_avg = {n: np.mean(v) for n, v in F_values.items() if v}
        enc = 2
        enc += 9 / F_avg.get(2, 1) if F_avg.get(2, 0) > 0 else 9
        enc += 1 / F_avg.get(3, 1) if F_avg.get(3, 0) > 0 else 1
        enc += 5 / F_avg.get(4, 1) if F_avg.get(4, 0) > 0 else 5
        enc += 3 / F_avg.get(6, 1) if F_avg.get(6, 0) > 0 else 3
        return min(61, max(20, enc))

    def calculate_cai(codon_counts, ref_codon_counts):
        optimal_freqs = {}
        for aa, codons in AA_TO_CODONS.items():
            if aa in ['M', 'W', '*']:
                continue
            ref_counts = [ref_codon_counts.get(c, 0) for c in codons]
            max_count = max(ref_counts) if ref_counts else 0
            total = sum(ref_counts)
            if total > 0 and max_count > 0:
                max_freq = max_count / total
                for i, codon in enumerate(codons):
                    freq = ref_counts[i] / total
                    optimal_freqs[codon] = freq / max_freq if max_freq > 0 else 0
        log_sum = 0
        n_codons = 0
        for codon, count in codon_counts.items():
            if codon in optimal_freqs and count > 0:
                w = optimal_freqs[codon]
                if w > 0:
                    log_sum += count * np.log(w)
                    n_codons += count
        if n_codons > 0:
            return np.exp(log_sum / n_codons)
        return np.nan

    def calculate_gc3s(codon_counts):
        gc3 = sum(codon_counts.get(c, 0) for c in SYNONYMOUS_CODONS if c[2] in 'GC')
        total = sum(codon_counts.get(c, 0) for c in SYNONYMOUS_CODONS)
        return gc3 / total if total > 0 else 0

    def calculate_gc_content(seq):
        seq = seq.upper()
        gc = sum(1 for b in seq if b in 'GC')
        total = sum(1 for b in seq if b in 'ATGC')
        return gc / total if total > 0 else 0

    def process_ffn(ffn_path, genome_id, ref_codon_counts):
        gene_counts = {}
        gene_sequences = {}
        current_gene = None
        current_seq = []

        with open(ffn_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_gene and current_seq:
                        seq = ''.join(current_seq).upper().replace('U', 'T')
                        gene_sequences[current_gene] = seq
                        counts = defaultdict(int)
                        for i in range(0, len(seq) - 2, 3):
                            codon = seq[i:i+3]
                            if all(b in 'TCAG' for b in codon):
                                counts[codon] += 1
                        gene_counts[current_gene] = dict(counts)
                    current_gene = line[1:].split()[0]
                    current_seq = []
                elif not line.startswith('#'):
                    current_seq.append(line)
            if current_gene and current_seq:
                seq = ''.join(current_seq).upper().replace('U', 'T')
                gene_sequences[current_gene] = seq
                counts = defaultdict(int)
                for i in range(0, len(seq) - 2, 3):
                    codon = seq[i:i+3]
                    if all(b in 'TCAG' for b in codon):
                        counts[codon] += 1
                gene_counts[current_gene] = dict(counts)

        results = []
        for gene_id, counts in gene_counts.items():
            total_codons = sum(counts.values())
            if total_codons < 30:
                continue
            enc = calculate_enc(counts)
            cai = calculate_cai(counts, ref_codon_counts)
            gc3s = calculate_gc3s(counts)
            gc_content = calculate_gc_content(gene_sequences.get(gene_id, ''))
            contig = '_'.join(gene_id.split('_')[:-1]) if '_' in gene_id else gene_id
            results.append({
                'genome_id': genome_id,
                'gene_id': gene_id,
                'total_codons': total_codons,
                'enc': enc,
                'cai': cai,
                'gc3s': gc3s,
                'gc_content': gc_content,
                'contig': contig
            })
        return results

    # =========================================================================
    # Main: compute metrics from pre-extracted FFN files
    # =========================================================================

    # Load reference for CAI calculation
    with open("${codon_reference}") as f:
        ref_data = json.load(f)
    ref_codon_counts = ref_data.get('codon_counts', {})
    print(f"Loaded ${genus} reference: {ref_data.get('n_genomes', 0)} genomes")

    all_results = []
    processed = 0

    for ffn in sorted(glob.glob("*.ffn")):
        genome_id = Path(ffn).stem
        gene_results = process_ffn(ffn, genome_id, ref_codon_counts)
        all_results.extend(gene_results)
        processed += 1

        if processed % 50 == 0:
            print(f"  Processed {processed} genomes, {len(all_results)} genes so far...")

    # Write output
    if all_results:
        df = pd.DataFrame(all_results)
        df.to_csv("${genus}_external_codon_metrics.tsv", sep='\\t', index=False)
    else:
        pd.DataFrame(columns=['genome_id', 'gene_id', 'total_codons', 'enc', 'cai', 'gc3s', 'gc_content', 'contig']) \\
            .to_csv("${genus}_external_codon_metrics.tsv", sep='\\t', index=False)

    print(f"${genus}: processed {processed} genomes, {len(all_results)} total genes")
    """
}
