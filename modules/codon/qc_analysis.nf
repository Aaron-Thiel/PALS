/*
 * Codon Bias QC Analysis
 *
 * Compares sample codon usage against a reference (species or genus level)
 * to detect potential assembly anomalies or contamination.
 */

process CODON_QC {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/codon_qc/${meta.id}", mode: 'copy'

    container 'quay.io/biocontainers/pandas:2.2.1'
    conda 'conda-forge::pandas conda-forge::numpy'

    input:
    tuple val(meta), path(ffn), path(reference)

    output:
    tuple val(meta), path("${meta.id}.codon_qc_summary.tsv"),      emit: summary
    tuple val(meta), path("${meta.id}.codon_qc_report.json"),      emit: report
    tuple val(meta), path("${meta.id}.gene_codon_analysis.tsv"),   emit: gene_analysis
    tuple val(meta), path("${meta.id}.contig_codon_analysis.tsv"), emit: contig_analysis
    tuple val(meta), path("${meta.id}.codon_frequencies.tsv"),     emit: codon_frequencies

    script:
    def sample_id = meta.id
    def genus = meta.genus ?: 'unknown'
    def species = meta.species ?: 'unknown'
    def ref_level = meta.ref_level ?: 'unknown'
    def zscore_threshold = params.codon_zscore_threshold ?: 2.5
    def min_contig_genes = params.codon_min_contig_genes ?: 5  // Minimum genes for contig-level analysis
    def fail_outlier_pct = params.codon_fail_outlier_pct ?: 15.0
    def warn_contig_outlier_frac = params.codon_warn_contig_outlier_frac ?: 0.2
    def has_reference = reference.name != 'NO_FILE'
    """
    #!/usr/bin/env python3
    import json
    import pandas as pd
    import numpy as np
    from collections import defaultdict
    from pathlib import Path

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

    # Amino acid full names for readability
    AA_NAMES = {
        'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
        'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
        'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
        'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
        '*': 'Stop'
    }

    def calculate_codon_frequencies(codon_counts):
        \"\"\"Calculate frequency of each codon within its amino acid family.\"\"\"
        freqs = {}
        for aa, codons in AA_TO_CODONS.items():
            total = sum(codon_counts.get(c, 0) for c in codons)
            for codon in codons:
                count = codon_counts.get(codon, 0)
                freqs[codon] = count / total if total > 0 else 0.0
        return freqs

    def calculate_rscu(codon_counts):
        rscu = {}
        for aa, codons in AA_TO_CODONS.items():
            if aa in ['M', 'W', '*']:
                continue
            total = sum(codon_counts.get(c, 0) for c in codons)
            expected = total / len(codons) if len(codons) > 0 else 0
            for codon in codons:
                observed = codon_counts.get(codon, 0)
                rscu[codon] = observed / expected if expected > 0 else 0.0
        return rscu

    def calculate_enc(codon_counts):
        \"\"\"Calculate Effective Number of Codons (ENC).\"\"\"
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

    def calculate_gc3s(codon_counts):
        \"\"\"Calculate GC content at 3rd codon position.\"\"\"
        gc3 = sum(codon_counts.get(c, 0) for c in SYNONYMOUS_CODONS if c[2] in 'GC')
        total = sum(codon_counts.get(c, 0) for c in SYNONYMOUS_CODONS)
        return gc3 / total if total > 0 else 0

    def euclidean_distance(rscu1, rscu2):
        \"\"\"Calculate Euclidean distance between two RSCU profiles.\"\"\"
        codons = [c for c in SYNONYMOUS_CODONS if c in rscu1 and c in rscu2]
        if not codons:
            return np.nan
        vec1 = np.array([rscu1[c] for c in codons])
        vec2 = np.array([rscu2[c] for c in codons])
        return np.sqrt(np.sum((vec1 - vec2) ** 2))

    def calculate_cai(codon_counts, ref_codon_counts):
        \"\"\"Calculate Codon Adaptation Index (CAI).
        Uses the most frequent codon per amino acid in reference as 'optimal'.
        CAI ranges from 0 to 1, where 1 means perfect match to reference preferences.
        \"\"\"
        # Find optimal (most frequent) codon per amino acid in reference
        optimal_freqs = {}
        for aa, codons in AA_TO_CODONS.items():
            if aa in ['M', 'W', '*']:
                continue
            ref_counts = [ref_codon_counts.get(c, 0) for c in codons]
            max_count = max(ref_counts) if ref_counts else 0
            total = sum(ref_counts)
            if total > 0 and max_count > 0:
                # Relative adaptiveness: w_i = freq_i / max_freq
                max_freq = max_count / total
                for i, codon in enumerate(codons):
                    freq = ref_counts[i] / total
                    optimal_freqs[codon] = freq / max_freq if max_freq > 0 else 0

        # Calculate CAI as geometric mean of relative adaptiveness
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

    def chi_square_test(observed_counts, expected_freqs):
        \"\"\"Chi-square test comparing observed codon counts to expected frequencies.
        Returns chi-square statistic and p-value.
        Uses Wilson-Hilferty approximation for p-value (no scipy needed).
        \"\"\"
        total = sum(observed_counts.get(c, 0) for c in SYNONYMOUS_CODONS)
        if total == 0:
            return np.nan, np.nan

        observed = []
        expected = []
        for codon in SYNONYMOUS_CODONS:
            obs = observed_counts.get(codon, 0)
            exp = expected_freqs.get(codon, 0) * total
            if exp > 0:  # Only include codons with expected counts
                observed.append(obs)
                expected.append(exp)

        if len(observed) < 2:
            return np.nan, np.nan

        # Calculate chi-square statistic
        observed = np.array(observed)
        expected = np.array(expected)
        chi2 = float(np.sum((observed - expected) ** 2 / expected))
        df = len(observed) - 1

        # Calculate p-value using Wilson-Hilferty approximation
        # For chi-square X with df degrees of freedom:
        # z = ((X/df)^(1/3) - (1 - 2/(9*df))) / sqrt(2/(9*df))
        # Then p-value ≈ 1 - Φ(z)
        if df <= 0:
            return chi2, np.nan

        h = 2.0 / (9.0 * df)
        z = ((chi2 / df) ** (1.0/3.0) - (1.0 - h)) / np.sqrt(h)

        # Normal CDF approximation using error function approximation
        # Φ(z) ≈ 0.5 * (1 + erf(z/√2))
        # erf approximation (Abramowitz and Stegun 7.1.26)
        def erf_approx(x):
            sign = 1 if x >= 0 else -1
            x = abs(x)
            a1, a2, a3, a4, a5 = 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429
            p = 0.3275911
            t = 1.0 / (1.0 + p * x)
            y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * np.exp(-x * x)
            return sign * y

        p_value = 1.0 - 0.5 * (1.0 + erf_approx(z / np.sqrt(2.0)))
        p_value = max(0.0, min(1.0, p_value))  # Clamp to [0, 1]

        return chi2, p_value

    def calculate_gc_content(seq):
        \"\"\"Calculate overall GC content of a sequence.\"\"\"
        seq = seq.upper()
        gc = sum(1 for b in seq if b in 'GC')
        total = sum(1 for b in seq if b in 'ATGC')
        return gc / total if total > 0 else 0

    # Count codons from sample FFN and store sequences for GC analysis
    sample_counts = defaultdict(int)
    gene_counts = {}
    gene_sequences = {}  # Store sequences for GC content calculation
    current_gene = None
    current_seq = []

    with open("${ffn}", 'r') as f:
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
                            sample_counts[codon] += 1
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
                    sample_counts[codon] += 1
            gene_counts[current_gene] = dict(counts)

    sample_rscu = calculate_rscu(dict(sample_counts))

    # Load reference
    ref_genomes = 0
    ref_name = 'self'
    has_ref = ${has_reference ? 'True' : 'False'}

    ref_codon_counts = {}
    ref_rscu_std = {}  # Standard deviation of RSCU across reference genomes (if available)
    ref_freq_std = {}  # Standard deviation of codon frequencies (if available)
    if has_ref:
        try:
            with open("${reference}") as f:
                ref_data = json.load(f)
            reference_rscu = ref_data.get('rscu', {})
            ref_codon_counts = ref_data.get('codon_counts', {})
            ref_rscu_std = ref_data.get('rscu_std', {})  # Load std if available
            ref_freq_std = ref_data.get('freq_std', {})  # Load freq std if available
            if not reference_rscu:
                reference_rscu = calculate_rscu(ref_codon_counts)
            ref_genomes = ref_data.get('n_genomes', 0)
            ref_name = ref_data.get('name', 'unknown')
            print(f"Loaded reference: {ref_name} ({ref_genomes} genomes)")
        except Exception as e:
            print(f"Warning: Could not load reference: {e}")
            reference_rscu = sample_rscu.copy()
            ref_codon_counts = dict(sample_counts)
    else:
        reference_rscu = sample_rscu.copy()
        ref_codon_counts = dict(sample_counts)
        print("Using self-comparison (no reference available)")

    # Calculate codon frequencies for detailed comparison
    sample_freqs = calculate_codon_frequencies(dict(sample_counts))
    ref_freqs = calculate_codon_frequencies(ref_codon_counts)

    # Analyze genes with enhanced metrics
    gene_results = []
    for gene_id, counts in gene_counts.items():
        total_codons = sum(counts.values())
        if total_codons < 30:
            continue
        gene_rscu = calculate_rscu(counts)
        dist = euclidean_distance(gene_rscu, reference_rscu)
        enc = calculate_enc(counts)
        gc3s = calculate_gc3s(counts)

        # CAI - Codon Adaptation Index
        cai = calculate_cai(counts, ref_codon_counts) if ref_codon_counts else np.nan

        # Chi-square test against reference frequencies
        chi2, chi2_pval = chi_square_test(counts, ref_freqs) if ref_freqs else (np.nan, np.nan)

        # Overall GC content of the gene
        gc_content = calculate_gc_content(gene_sequences.get(gene_id, ''))

        gene_results.append({
            'gene_id': gene_id,
            'total_codons': total_codons,
            'euclidean_distance': dist,
            'enc': enc,
            'gc3s': gc3s,
            'gc_content': gc_content,
            'cai': cai,
            'chi2': chi2,
            'chi2_pval': chi2_pval,
            'contig': '_'.join(gene_id.split('_')[:-1]) if '_' in gene_id else gene_id
        })

    gene_df = pd.DataFrame(gene_results) if gene_results else pd.DataFrame()

    # Flag outliers using z-score
    zscore_threshold = ${zscore_threshold}
    if not gene_df.empty and len(gene_df) > 2:
        mean_dist = gene_df['euclidean_distance'].mean()
        std_dist = gene_df['euclidean_distance'].std()
        gene_df['distance_zscore'] = (gene_df['euclidean_distance'] - mean_dist) / std_dist if std_dist > 0 else 0
        gene_df['is_outlier'] = gene_df['distance_zscore'].abs() > zscore_threshold
    else:
        if not gene_df.empty:
            gene_df['distance_zscore'] = 0
            gene_df['is_outlier'] = False

    # Contig-level analysis with weighted aggregation
    min_contig_genes = ${min_contig_genes}
    warn_contig_outlier_frac = ${warn_contig_outlier_frac}
    if not gene_df.empty:
        # Weighted aggregation by gene length (total_codons)
        def weighted_mean(group, col):
            weights = group['total_codons']
            return np.average(group[col], weights=weights) if weights.sum() > 0 else group[col].mean()

        contig_groups = gene_df.groupby('contig')
        contig_rows = []
        for contig, group in contig_groups:
            n_genes = len(group)
            total_codons = group['total_codons'].sum()

            # Weighted means for better accuracy
            mean_distance = weighted_mean(group, 'euclidean_distance')
            mean_enc = weighted_mean(group, 'enc')
            mean_gc3s = weighted_mean(group, 'gc3s')
            mean_gc_content = weighted_mean(group, 'gc_content')
            mean_cai = weighted_mean(group, 'cai') if 'cai' in group.columns and not group['cai'].isna().all() else np.nan

            # GC correlation: correlation between gc_content and gc3s within contig
            # High correlation suggests natural variation; low correlation may indicate contamination
            if n_genes >= 5:
                gc_corr = group['gc_content'].corr(group['gc3s'])
            else:
                gc_corr = np.nan

            # Outlier analysis
            n_outlier_genes = group['is_outlier'].sum()
            outlier_fraction = n_outlier_genes / n_genes

            contig_rows.append({
                'contig': contig,
                'n_genes': n_genes,
                'total_codons': total_codons,
                'mean_distance': mean_distance,
                'std_distance': group['euclidean_distance'].std(),
                'mean_enc': mean_enc,
                'mean_gc3s': mean_gc3s,
                'mean_gc_content': mean_gc_content,
                'gc_correlation': gc_corr,
                'mean_cai': mean_cai,
                'n_outlier_genes': n_outlier_genes,
                'outlier_fraction': outlier_fraction,
                # Suspicious if high outlier fraction AND enough genes for reliable assessment
                'is_suspicious': (outlier_fraction > warn_contig_outlier_frac) and (n_genes >= min_contig_genes)
            })
        contig_df = pd.DataFrame(contig_rows)
    else:
        contig_df = pd.DataFrame()

    # Determine QC status using parameterized thresholds
    fail_outlier_pct = ${fail_outlier_pct}
    total_genes = len(gene_df)
    outlier_genes = int(gene_df['is_outlier'].sum()) if not gene_df.empty else 0
    outlier_pct = (outlier_genes / total_genes * 100) if total_genes > 0 else 0
    suspicious_contigs = int(contig_df['is_suspicious'].sum()) if not contig_df.empty else 0

    # Calculate sample-level CAI and chi-square
    sample_cai = calculate_cai(dict(sample_counts), ref_codon_counts) if ref_codon_counts else np.nan
    sample_chi2, sample_chi2_pval = chi_square_test(dict(sample_counts), ref_freqs) if ref_freqs else (np.nan, np.nan)

    # Sample-level GC correlation
    if not gene_df.empty and len(gene_df) >= 10:
        sample_gc_corr = gene_df['gc_content'].corr(gene_df['gc3s'])
    else:
        sample_gc_corr = np.nan

    if outlier_pct > fail_outlier_pct:
        qc_status = 'FAIL'
        qc_message = f'High outlier gene percentage ({outlier_pct:.1f}% > {fail_outlier_pct}%) - potential assembly issues'
    elif suspicious_contigs > 0:
        qc_status = 'WARNING'
        qc_message = f'{suspicious_contigs} suspicious contigs detected'
    else:
        qc_status = 'PASS'
        qc_message = 'Codon usage within expected range'

    # Save outputs
    gene_df.to_csv("${sample_id}.gene_codon_analysis.tsv", sep='\\t', index=False)
    if not contig_df.empty:
        contig_df.to_csv("${sample_id}.contig_codon_analysis.tsv", sep='\\t', index=False)
    else:
        open("${sample_id}.contig_codon_analysis.tsv", 'w').close()

    # Generate detailed codon frequency comparison (one row per codon)
    # Excludes M, W (single codon - no synonymous choice) and stop codons
    codon_freq_rows = []
    for aa in sorted(AA_TO_CODONS.keys()):
        if aa in ['M', 'W', '*']:  # Skip single-codon AAs and stop codons
            continue
        codons = sorted(AA_TO_CODONS[aa])
        aa_sample_total = sum(sample_counts.get(c, 0) for c in codons)
        aa_ref_total = sum(ref_codon_counts.get(c, 0) for c in codons)
        for codon in codons:
            sample_count = sample_counts.get(codon, 0)
            sample_pct = sample_freqs.get(codon, 0) * 100
            ref_pct = ref_freqs.get(codon, 0) * 100
            diff_pct = sample_pct - ref_pct

            # Reference std (if available from multi-genome reference)
            rscu_std = ref_rscu_std.get(codon, None)
            freq_std = ref_freq_std.get(codon, None)

            # Calculate z-score if we have reference std
            rscu_zscore = None
            if rscu_std and rscu_std > 0:
                rscu_zscore = (sample_rscu.get(codon, 0) - reference_rscu.get(codon, 0)) / rscu_std

            codon_freq_rows.append({
                'amino_acid': aa,
                'aa_name': AA_NAMES.get(aa, aa),
                'codon': codon,
                'n_synonymous': len(codons),
                'sample_count': sample_count,
                'sample_aa_total': aa_sample_total,
                'sample_pct': round(sample_pct, 2),
                'ref_pct': round(ref_pct, 2),
                'ref_pct_std': round(freq_std * 100, 2) if freq_std else None,
                'diff_pct': round(diff_pct, 2),
                'rscu_sample': round(sample_rscu.get(codon, 0), 3),
                'rscu_ref': round(reference_rscu.get(codon, 0), 3),
                'rscu_ref_std': round(rscu_std, 3) if rscu_std else None,
                'rscu_zscore': round(rscu_zscore, 2) if rscu_zscore is not None else None
            })
    codon_freq_df = pd.DataFrame(codon_freq_rows)
    codon_freq_df.to_csv("${sample_id}.codon_frequencies.tsv", sep='\\t', index=False)

    # Count contigs with enough genes for suspicious classification
    contigs_assessed = int((contig_df['n_genes'] >= min_contig_genes).sum()) if not contig_df.empty else 0

    summary = {
        'sample_id': '${sample_id}',
        'genus': '${genus}',
        'species': '${species}',
        'ref_level': '${ref_level}',
        'ref_name': ref_name,
        'qc_status': qc_status,
        'total_genes': total_genes,
        'outlier_genes': outlier_genes,
        'outlier_pct': round(outlier_pct, 2),
        'mean_enc': round(gene_df['enc'].mean(), 2) if not gene_df.empty else 0,
        'mean_gc3s': round(gene_df['gc3s'].mean(), 4) if not gene_df.empty else 0,
        'mean_gc_content': round(gene_df['gc_content'].mean(), 4) if not gene_df.empty else 0,
        'cai': round(sample_cai, 4) if not np.isnan(sample_cai) else None,
        'chi2': round(sample_chi2, 2) if not np.isnan(sample_chi2) else None,
        'chi2_pval': round(sample_chi2_pval, 4) if not np.isnan(sample_chi2_pval) else None,
        'gc_correlation': round(sample_gc_corr, 4) if not np.isnan(sample_gc_corr) else None,
        'suspicious_contigs': suspicious_contigs,
        'contigs_assessed': contigs_assessed,
        'ref_genomes_used': ref_genomes
    }

    pd.DataFrame([summary]).to_csv("${sample_id}.codon_qc_summary.tsv", sep='\\t', index=False)

    report = {
        'sample_id': '${sample_id}',
        'genus': '${genus}',
        'species': '${species}',
        'ref_level': '${ref_level}',
        'ref_name': ref_name,
        'qc_status': qc_status,
        'qc_message': qc_message,
        'analysis_summary': summary,
        'suspicious_contigs': contig_df[contig_df['is_suspicious']].to_dict('records') if not contig_df.empty else []
    }
    with open("${sample_id}.codon_qc_report.json", 'w') as f:
        json.dump(report, f, indent=2)

    print(f"Codon QC: ${sample_id} - {qc_status} ({outlier_genes}/{total_genes} outlier genes)")
    """
}
