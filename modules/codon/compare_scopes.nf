/*
 * Compare Internal vs External Codon Correlation Scopes
 *
 * Mann-Whitney U tests with Bonferroni correction comparing
 * the same pangenome category (core/accessory/rare) between
 * internal and external scope classifications.
 */

process COMPARE_SCOPES {
    tag "$group_id"
    label 'process_single'

    conda 'conda-forge::pandas conda-forge::numpy conda-forge::scipy'
    container null

    input:
    tuple val(group_id), val(output_subdir)
    path internal_merged_csv
    path external_merged_csv

    output:
    tuple val(group_id), path("comparison.txt"),  emit: comparison_txt,  optional: true
    tuple val(group_id), path("comparison.json"), emit: comparison_json, optional: true

    script:
    """
    #!/usr/bin/env python3
    import json
    import os
    import pandas as pd
    import numpy as np
    from scipy import stats as sp_stats

    GROUP_ID = "${group_id}"
    print(f"=== Scope comparison: {GROUP_ID} ===")

    # Load merged data from internal and external scopes
    try:
        internal = pd.read_csv("${internal_merged_csv}")
        external = pd.read_csv("${external_merged_csv}")
    except Exception as e:
        print(f"  Cannot load data: {e}")
        with open('comparison.txt', 'w') as f:
            f.write(f"Scope comparison skipped for {GROUP_ID}: {e}\\n")
        with open('comparison.json', 'w') as f:
            json.dump({'error': str(e)}, f)
        import sys
        sys.exit(0)

    metrics = ['enc', 'cai', 'gc3s', 'gc_content']
    categories = ['core', 'accessory', 'rare']
    metric_names = {'enc': 'ENC', 'cai': 'CAI', 'gc3s': 'GC3s', 'gc_content': 'GC Content'}

    lines = []
    lines.append('=' * 70)
    lines.append('INTERNAL vs EXTERNAL SCOPE COMPARISON')
    lines.append('=' * 70)
    lines.append(f'Group: {GROUP_ID}')
    lines.append(f'Test: Mann-Whitney U (two-sided)')
    lines.append(f'Correction: Bonferroni (per metric, {len(categories)} comparisons)')
    lines.append('')

    # Gene counts per category
    n_int = internal['category'].value_counts() if 'category' in internal.columns else pd.Series(dtype=int)
    n_ext = external['category'].value_counts() if 'category' in external.columns else pd.Series(dtype=int)
    lines.append('Gene counts per category:')
    lines.append(f"  {'Category':<12} {'Internal':>10} {'External':>10}")
    lines.append(f"  {'-'*12} {'-'*10} {'-'*10}")
    for cat in categories:
        ni = n_int.get(cat, 0)
        ne = n_ext.get(cat, 0)
        lines.append(f"  {cat:<12} {ni:>10,} {ne:>10,}")
    lines.append('')

    n_comparisons = len(categories)
    comparison_results = {}

    for metric in metrics:
        if metric not in internal.columns or metric not in external.columns:
            continue
        lines.append(f"\\n{metric_names[metric]} ({metric})")
        lines.append('-' * 50)
        lines.append(f"  {'Category':<12} {'Int mean':>10} {'Ext mean':>10} {'Diff':>10} {'U stat':>14} {'p-value':>12} {'p(adj)':>12} {'Sig':>5}")

        comparison_results[metric] = {}
        for cat in categories:
            int_vals = internal.loc[internal['category'] == cat, metric].dropna()
            ext_vals = external.loc[external['category'] == cat, metric].dropna()
            if len(int_vals) < 5 or len(ext_vals) < 5:
                lines.append(f"  {cat:<12} {'n/a (too few observations)':>60}")
                continue
            int_mean = int_vals.mean()
            ext_mean = ext_vals.mean()
            diff = ext_mean - int_mean
            u_stat, p_val = sp_stats.mannwhitneyu(int_vals, ext_vals, alternative='two-sided')
            p_adj = min(p_val * n_comparisons, 1.0)
            sig = '***' if p_adj < 0.001 else '**' if p_adj < 0.01 else '*' if p_adj < 0.05 else ''
            lines.append(f"  {cat:<12} {int_mean:>10.4f} {ext_mean:>10.4f} {diff:>+10.4f} {u_stat:>14.0f} {p_val:>12.2e} {p_adj:>12.2e} {sig:>5}")
            comparison_results[metric][cat] = {
                'int_mean': float(int_mean), 'ext_mean': float(ext_mean),
                'diff': float(diff), 'u_stat': float(u_stat),
                'p_value': float(p_val), 'p_adj': float(p_adj),
                'significant': bool(p_adj < 0.05)
            }

    # Effect sizes
    lines.append(f"\\n\\nEFFECT SIZES (rank-biserial correlation r)")
    lines.append('=' * 70)
    lines.append('  r ~ 0.1: small, r ~ 0.3: medium, r ~ 0.5: large')
    lines.append(f"\\n  {'Category':<12} {'ENC':>10} {'CAI':>10} {'GC3s':>10} {'GC Content':>10}")
    lines.append(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
    for cat in categories:
        row_vals = []
        for metric in metrics:
            if metric not in internal.columns or metric not in external.columns:
                row_vals.append(f"{'n/a':>10}")
                continue
            int_vals = internal.loc[internal['category'] == cat, metric].dropna()
            ext_vals = external.loc[external['category'] == cat, metric].dropna()
            if len(int_vals) < 5 or len(ext_vals) < 5:
                row_vals.append(f"{'n/a':>10}")
                continue
            u_stat, _ = sp_stats.mannwhitneyu(int_vals, ext_vals, alternative='two-sided')
            r = 1 - (2 * u_stat) / (len(int_vals) * len(ext_vals))
            row_vals.append(f"{r:>+10.4f}")
        lines.append(f"  {cat:<12} {''.join(row_vals)}")

    lines.append(f"\\n\\n{'='*70}")
    lines.append('INTERPRETATION')
    lines.append('=' * 70)
    lines.append('')
    lines.append('A significant difference means the same genes are classified')
    lines.append('differently depending on whether only internal (PAL) or only external')
    lines.append('(GCA) samples define core/accessory/rare boundaries.')
    lines.append('')
    lines.append('Large differences suggest the internal sample set is not representative')
    lines.append('of the broader genus/species diversity.')
    lines.append('')
    lines.append('Small or non-significant differences indicate the internal collection')
    lines.append('captures similar pangenomic structure as the external references.')

    with open('comparison.txt', 'w') as f:
        f.write('\\n'.join(lines))
    with open('comparison.json', 'w') as f:
        json.dump(comparison_results, f, indent=2, default=str)

    print(f"=== Done: {GROUP_ID} scope comparison ===")
    """
}
