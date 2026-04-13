/*
 * Codon Bias vs Pangenomics Correlation Analysis
 *
 * Analyzes the relationship between codon bias metrics (ENC, CAI, GC3s, GC)
 * and pangenome gene categories (core, accessory, rare).
 *
 * Runs once per (group_id, scope) combination.
 * Adapted from test/codon_pangenomics/codon_pangenomics_correlation.py
 */

process CODON_CORRELATION {
    tag "$group_id / $scope"
    label 'process_low'

    conda 'conda-forge::pandas conda-forge::numpy conda-forge::matplotlib conda-forge::seaborn conda-forge::scipy'
    container null

    input:
    tuple val(group_id), val(scope), val(output_subdir)
    path rtab_file
    path gpa_csv
    path internal_codon_files
    path reference_codon_files
    path pipeline_summary
    path taxonomy_csv
    val rare_threshold
    val core_threshold

    output:
    tuple val(group_id), val(scope), path("*.tsv"),  emit: stats,   optional: true
    tuple val(group_id), val(scope), path("*.json"), emit: tests,   optional: true
    tuple val(group_id), val(scope), path("*.csv"),  emit: data,    optional: true
    tuple val(group_id), val(scope), path("*.txt"),  emit: report,  optional: true
    tuple val(group_id), val(scope), path("*.png"),  emit: png,     optional: true
    tuple val(group_id), val(scope), path("*.pdf"),  emit: pdf,     optional: true

    script:
    """
    #!/usr/bin/env python3
    import json
    import os
    import glob
    from pathlib import Path
    from collections import defaultdict

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    from scipy import stats as sp_stats

    # Style setup (consistent with compare_assemblies.py)
    plt.rcParams['figure.figsize'] = (12, 8)
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.labelsize'] = 11
    sns.set_style("whitegrid")
    sns.set_palette("husl")

    # =========================================================================
    # Parameters from Nextflow
    # =========================================================================
    GROUP_ID = "${group_id}"
    SCOPE = "${scope}"
    RARE_THRESHOLD = ${rare_threshold}
    CORE_THRESHOLD = ${core_threshold}

    print(f"=== Codon Correlation: {GROUP_ID} / {SCOPE} scope ===")

    # =========================================================================
    # 1. Load taxonomy
    # =========================================================================
    def load_taxonomy():
        records = []
        # Internal from pipeline_summary.tsv
        ps_path = "${pipeline_summary}"
        if os.path.exists(ps_path) and os.path.getsize(ps_path) > 0:
            df_int = pd.read_csv(ps_path, sep='\\t', usecols=['sample_id', 'family', 'genus', 'species'])
            for _, row in df_int.iterrows():
                genus = row.get('genus', '')
                if pd.notna(genus) and genus not in ('NA', ''):
                    species = row.get('species', '')
                    if pd.isna(species) or species == 'NA':
                        species = ''
                    family = row.get('family', '')
                    if pd.isna(family) or family == 'NA':
                        family = ''
                    records.append({'sample_id': row['sample_id'], 'family': family,
                                    'genus': genus, 'species': species, 'source': 'internal'})
            print(f"  Internal taxonomy: {len(records)} samples")

        # External from taxonomy.csv
        n_before = len(records)
        tax_path = "${taxonomy_csv}"
        if os.path.exists(tax_path) and os.path.getsize(tax_path) > 0:
            df_ext = pd.read_csv(tax_path)
            for _, row in df_ext.iterrows():
                genus = row.get('genus', '')
                if pd.notna(genus) and genus not in ('',):
                    species = str(row.get('species', '')).replace(' ', '_')
                    if pd.isna(row.get('species', '')):
                        species = ''
                    family = row.get('family', '')
                    if pd.isna(family):
                        family = ''
                    records.append({'sample_id': row['genome_id'], 'family': family,
                                    'genus': genus, 'species': species, 'source': 'external'})
            print(f"  External taxonomy: {len(records) - n_before} samples")

        return pd.DataFrame(records) if records else pd.DataFrame(columns=['sample_id', 'family', 'genus', 'species', 'source'])

    taxonomy_df = load_taxonomy()

    # =========================================================================
    # 2. Determine sample filters based on group_id and scope
    # =========================================================================
    # Parse group_id to get taxonomy filter
    taxonomy_filter_level = None
    taxonomy_filter_name = None
    if GROUP_ID.startswith('genus_'):
        taxonomy_filter_level = 'genus'
        taxonomy_filter_name = GROUP_ID[6:]
    elif GROUP_ID.startswith('species_'):
        taxonomy_filter_level = 'species'
        taxonomy_filter_name = GROUP_ID[8:]

    # Filter taxonomy by group
    if taxonomy_filter_level and not taxonomy_df.empty:
        mask = taxonomy_df[taxonomy_filter_level] == taxonomy_filter_name
        group_taxonomy = taxonomy_df[mask]
        print(f"  Group filter: {taxonomy_filter_level}={taxonomy_filter_name} -> {len(group_taxonomy)} samples")
    else:
        group_taxonomy = taxonomy_df

    # Get sample lists by source
    internal_sample_ids = set(group_taxonomy[group_taxonomy['source'] == 'internal']['sample_id'])
    external_sample_ids = set(group_taxonomy[group_taxonomy['source'] == 'external']['sample_id'])
    all_sample_ids = internal_sample_ids | external_sample_ids

    # Use internal samples for both Rtab classification and codon data
    classification_samples = internal_sample_ids
    print(f"  Internal samples: {len(classification_samples)} classification samples")

    if len(classification_samples) < 3:
        print(f"  SKIP: Too few classification samples ({len(classification_samples)})")
        # Write empty marker
        with open('analysis_report.txt', 'w') as f:
            f.write(f"Skipped: too few samples for {GROUP_ID} / {SCOPE} ({len(classification_samples)} samples)\\n")
        import sys
        sys.exit(0)

    # =========================================================================
    # 3. Load and filter Rtab
    # =========================================================================
    print("Loading Rtab...")
    header = pd.read_csv("${rtab_file}", sep='\\t', nrows=0)
    all_rtab_cols = header.columns.tolist()

    # Filter to classification samples that exist in Rtab
    cols_to_use = [all_rtab_cols[0]] + [c for c in all_rtab_cols[1:] if c in classification_samples]
    if len(cols_to_use) <= 1:
        print("  SKIP: No matching samples in Rtab")
        with open('analysis_report.txt', 'w') as f:
            f.write(f"Skipped: no matching samples in Rtab for {GROUP_ID} / {SCOPE}\\n")
        import sys
        sys.exit(0)

    rtab = pd.read_csv("${rtab_file}", sep='\\t', index_col=0, usecols=cols_to_use).astype('int8')
    n_classification = rtab.shape[1]
    print(f"  Rtab: {rtab.shape[0]} genes x {n_classification} samples")

    # Classify genes
    n_genomes = rtab.shape[1]
    gene_counts = rtab.sum(axis=1)
    gene_freq = gene_counts / n_genomes

    gene_categories = pd.DataFrame({
        'gene_id': rtab.index,
        'n_genomes_present': gene_counts,
        'frequency': gene_freq
    }).set_index('gene_id')

    conditions = [
        gene_freq >= CORE_THRESHOLD,
        (gene_freq >= RARE_THRESHOLD) & (gene_freq < CORE_THRESHOLD),
        gene_freq < RARE_THRESHOLD
    ]
    choices = ['core', 'accessory', 'rare']
    gene_categories['category'] = np.select(conditions, choices, default='unknown')

    cat_counts = gene_categories['category'].value_counts()
    print(f"  Gene categories:")
    for cat in ['core', 'accessory', 'rare']:
        if cat in cat_counts.index:
            print(f"    {cat}: {cat_counts[cat]:,}")
    del rtab

    # =========================================================================
    # 4. Load codon data (internal + external as needed)
    # =========================================================================
    print("Loading codon data...")
    all_codon_data = []

    # Internal codon files (from CODON_QC: {sample_id}.gene_codon_analysis.tsv)
    for f in sorted(glob.glob("*.gene_codon_analysis.tsv")):
        sample_id = Path(f).stem.replace('.gene_codon_analysis', '')
        if taxonomy_filter_level and sample_id not in internal_sample_ids:
            continue
        try:
            df = pd.read_csv(f, sep='\\t')
            df['sample_id'] = sample_id
            df['data_source'] = 'internal'
            all_codon_data.append(df)
        except Exception as e:
            print(f"  Warning: Could not load {f}: {e}")

    if not all_codon_data:
        print("  SKIP: No codon data loaded")
        with open('analysis_report.txt', 'w') as f:
            f.write(f"Skipped: no codon data for {GROUP_ID} / {SCOPE}\\n")
        import sys
        sys.exit(0)

    codon_data = pd.concat(all_codon_data, ignore_index=True)
    n_codon_samples = codon_data['sample_id'].nunique()
    print(f"  Loaded internal codon data: {len(codon_data)} genes from {n_codon_samples} samples")

    # Load reference codon files (*_external_codon_metrics.tsv)
    ref_codon_data = []
    for f in sorted(glob.glob("*_external_codon_metrics.tsv")):
        try:
            df = pd.read_csv(f, sep='\\t')
            if taxonomy_filter_level and not df.empty:
                df = df[df['genome_id'].isin(external_sample_ids)]
            if not df.empty:
                df = df.rename(columns={'genome_id': 'sample_id'})
                df['data_source'] = 'reference'
                ref_codon_data.append(df)
        except Exception as e:
            print(f"  Warning: Could not load {f}: {e}")

    has_reference = len(ref_codon_data) > 0
    if has_reference:
        ref_codon_df = pd.concat(ref_codon_data, ignore_index=True)
        n_ref_samples = ref_codon_df['sample_id'].nunique()
        print(f"  Loaded reference codon data: {len(ref_codon_df)} genes from {n_ref_samples} samples")
    else:
        ref_codon_df = pd.DataFrame()
        print("  No reference codon data found (comparison plots will be skipped)")

    # =========================================================================
    # 5. Map genes to clusters via gene_presence_absence.csv
    # =========================================================================
    print("Building gene-to-cluster mapping...")
    gpa_path = "${gpa_csv}"
    locus_to_cluster = {}

    if os.path.exists(gpa_path) and os.path.getsize(gpa_path) > 0:
        # Include both internal and reference sample IDs for GPA lookup
        codon_sample_ids = set(codon_data['sample_id'].unique())
        if has_reference:
            codon_sample_ids |= set(ref_codon_df['sample_id'].unique())
        gpa_header = pd.read_csv(gpa_path, nrows=0)
        gpa_all_cols = gpa_header.columns.tolist()
        sample_cols = [c for c in gpa_all_cols if c in codon_sample_ids]
        cols_to_read = ['Gene'] + sample_cols

        if sample_cols:
            for chunk in pd.read_csv(gpa_path, usecols=cols_to_read, chunksize=50000, dtype=str):
                for _, row in chunk.iterrows():
                    cluster_id = row['Gene']
                    for sample in sample_cols:
                        gene_ids = row[sample]
                        if pd.notna(gene_ids) and gene_ids != '':
                            for gid in str(gene_ids).split('\\t'):
                                gid = gid.strip()
                                if gid:
                                    parts = gid.split('-')
                                    locus_tag = parts[-1] if len(parts) >= 3 else gid
                                    locus_to_cluster[locus_tag] = cluster_id
            print(f"  Built mapping for {len(locus_to_cluster)} locus tags from {len(sample_cols)} samples")
        else:
            print("  Warning: No matching sample columns in GPA CSV")
    else:
        print("  Warning: gene_presence_absence.csv not found")

    # Map codon data to clusters
    codon_data['cluster_id'] = codon_data['gene_id'].map(locus_to_cluster)
    n_mapped = codon_data['cluster_id'].notna().sum()
    print(f"  Mapped {n_mapped}/{len(codon_data)} genes to clusters ({n_mapped/max(len(codon_data),1)*100:.1f}%)")

    # Merge with gene categories
    merged = codon_data.merge(
        gene_categories.reset_index()[['gene_id', 'category', 'frequency']],
        left_on='cluster_id', right_on='gene_id',
        how='left', suffixes=('', '_cluster')
    )
    n_categorized = merged['category'].notna().sum()
    print(f"  Categorized {n_categorized}/{len(merged)} genes ({n_categorized/max(len(merged),1)*100:.1f}%)")

    merged_data = merged[merged['category'].isin(['core', 'accessory', 'rare'])].copy()
    if len(merged_data) == 0:
        print("  SKIP: No genes with categories after merge")
        with open('analysis_report.txt', 'w') as f:
            f.write(f"Skipped: no genes with categories for {GROUP_ID} / {SCOPE}\\n")
        import sys
        sys.exit(0)

    print(f"  Final dataset: {len(merged_data)} categorized genes")

    # =========================================================================
    # 6. Statistics
    # =========================================================================
    print("Running statistical tests...")
    metrics = ['enc', 'cai', 'gc3s', 'gc_content']
    available_metrics = [m for m in metrics if m in merged_data.columns]
    categories = ['core', 'accessory', 'rare']

    # Category stats
    stats_rows = []
    for cat in categories:
        cat_data = merged_data[merged_data['category'] == cat]
        if len(cat_data) == 0:
            continue
        row = {'category': cat, 'n_genes': len(cat_data)}
        for metric in available_metrics:
            vals = cat_data[metric].dropna()
            if len(vals) > 0:
                row[f'{metric}_mean'] = vals.mean()
                row[f'{metric}_median'] = vals.median()
                row[f'{metric}_std'] = vals.std()
                row[f'{metric}_q25'] = vals.quantile(0.25)
                row[f'{metric}_q75'] = vals.quantile(0.75)
        stats_rows.append(row)
    stats_df = pd.DataFrame(stats_rows)
    if not stats_df.empty:
        stats_df.to_csv('category_statistics.tsv', sep='\\t', index=False)

    # Statistical tests
    test_results = {}
    for metric in available_metrics:
        test_results[metric] = {}
        cat_data = {}
        for cat in categories:
            vals = merged_data[merged_data['category'] == cat][metric].dropna()
            if len(vals) > 0:
                cat_data[cat] = vals
        if len(cat_data) < 2:
            continue
        # Kruskal-Wallis
        groups = list(cat_data.values())
        if all(len(g) >= 3 for g in groups):
            try:
                h_stat, p_val = sp_stats.kruskal(*groups)
                test_results[metric]['kruskal_wallis'] = {'statistic': float(h_stat), 'p_value': float(p_val), 'significant': bool(p_val < 0.05)}
            except Exception as e:
                test_results[metric]['kruskal_wallis'] = {'error': str(e)}
        # Pairwise Mann-Whitney U
        test_results[metric]['pairwise'] = {}
        for cat1, cat2 in [('core', 'accessory'), ('core', 'rare'), ('accessory', 'rare')]:
            if cat1 in cat_data and cat2 in cat_data and len(cat_data[cat1]) >= 3 and len(cat_data[cat2]) >= 3:
                try:
                    u_stat, p_val = sp_stats.mannwhitneyu(cat_data[cat1], cat_data[cat2], alternative='two-sided')
                    test_results[metric]['pairwise'][f'{cat1}_vs_{cat2}'] = {'statistic': float(u_stat), 'p_value': float(p_val), 'significant': bool(p_val < 0.05)}
                except Exception as e:
                    test_results[metric]['pairwise'][f'{cat1}_vs_{cat2}'] = {'error': str(e)}
        # Spearman correlation with frequency
        if 'frequency' in merged_data.columns:
            valid = merged_data[[metric, 'frequency']].dropna()
            if len(valid) >= 10:
                try:
                    rho, p_val = sp_stats.spearmanr(valid[metric], valid['frequency'])
                    test_results[metric]['frequency_correlation'] = {'spearman_rho': float(rho), 'p_value': float(p_val), 'significant': bool(p_val < 0.05)}
                except Exception as e:
                    test_results[metric]['frequency_correlation'] = {'error': str(e)}

    with open('statistical_tests.json', 'w') as f:
        json.dump(test_results, f, indent=2, default=str)

    # =========================================================================
    # 7. Plots
    # =========================================================================
    print("Generating plots...")
    category_colors = {'core': '#2ecc71', 'accessory': '#3498db', 'rare': '#e74c3c'}
    category_order = ['core', 'accessory', 'rare']
    metric_labels = {
        'enc': 'ENC (Effective Number of Codons)',
        'cai': 'CAI (Codon Adaptation Index)',
        'gc3s': 'GC3s (GC content at 3rd position)',
        'gc_content': 'GC Content'
    }

    # Build title suffix
    parts = []
    if taxonomy_filter_name:
        parts.append(taxonomy_filter_name.replace('_', ' '))
    parts.append(f"{SCOPE} scope")
    title_suffix = f"\\n({', '.join(parts)})"

    present_cats = [c for c in category_order if c in merged_data['category'].values]

    # --- Violin plots ---
    if available_metrics and present_cats:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        for i, metric in enumerate(available_metrics[:4]):
            ax = axes[i]
            plot_parts = ax.violinplot(
                [merged_data[merged_data['category'] == cat][metric].dropna() for cat in present_cats],
                positions=range(len(present_cats)), showmeans=True, showmedians=True
            )
            for j, cat in enumerate(present_cats):
                if j < len(plot_parts['bodies']):
                    plot_parts['bodies'][j].set_facecolor(category_colors[cat])
                    plot_parts['bodies'][j].set_alpha(0.7)
            ax.set_xticks(range(len(present_cats)))
            ax.set_xticklabels([c.capitalize() for c in present_cats])
            ax.set_ylabel(metric_labels.get(metric, metric))
            ax.set_title(f'{metric_labels.get(metric, metric)} by Gene Category')
            ax.grid(axis='y', alpha=0.3)
        for i in range(len(available_metrics), 4):
            axes[i].set_visible(False)
        plt.suptitle(f'Codon Bias Metrics vs Pangenomics Categories{title_suffix}', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig('codon_bias_by_category_violin.png', dpi=150, bbox_inches='tight')
        plt.savefig('codon_bias_by_category_violin.pdf', bbox_inches='tight')
        plt.close()

    # --- Box plots ---
    if available_metrics and present_cats:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        box_metric_labels = {'enc': 'ENC\\n(lower = stronger bias)', 'cai': 'CAI\\n(higher = better adapted)', 'gc3s': 'GC3s', 'gc_content': 'GC Content'}
        for i, metric in enumerate(available_metrics[:4]):
            ax = axes[i]
            subset = merged_data[merged_data['category'].isin(present_cats)]
            sns.boxplot(data=subset, x='category', y=metric, hue='category', order=present_cats, hue_order=present_cats, palette=category_colors, ax=ax, showfliers=False, legend=False)
            max_points = 500
            sampled = pd.concat([g.sample(n=min(len(g), max_points), random_state=42) for _, g in subset.groupby('category')]).reset_index(drop=True)
            sns.stripplot(data=sampled, x='category', y=metric, order=present_cats, color='black', alpha=0.2, size=1.5, ax=ax)
            ax.set_xlabel('Gene Category')
            ax.set_ylabel(box_metric_labels.get(metric, metric))
            ax.grid(axis='y', alpha=0.3)
        for i in range(len(available_metrics), 4):
            axes[i].set_visible(False)
        plt.suptitle(f'Codon Bias Distribution by Gene Conservation Level{title_suffix}', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig('codon_bias_by_category_boxplot.png', dpi=150, bbox_inches='tight')
        plt.savefig('codon_bias_by_category_boxplot.pdf', bbox_inches='tight')
        plt.close()

    # --- Individual box plots ---
    metric_ylim = {'enc': (25, 65), 'cai': (0.5, 0.9), 'gc3s': (0.2, 0.7), 'gc_content': (0.2, 0.7)}
    for metric in available_metrics:
        fig, ax = plt.subplots(figsize=(5, 5))
        subset = merged_data[merged_data['category'].isin(present_cats)]
        sns.boxplot(data=subset, x='category', y=metric, hue='category', order=present_cats, hue_order=present_cats, palette=category_colors, ax=ax, showfliers=False, legend=False)
        sampled = pd.concat([g.sample(n=min(len(g), 500), random_state=42) for _, g in subset.groupby('category')]).reset_index(drop=True)
        sns.stripplot(data=sampled, x='category', y=metric, order=present_cats, color='black', alpha=0.2, size=1.5, ax=ax)
        ax.set_xlabel('Gene Category')
        ax.set_ylabel(metric_labels.get(metric, metric))
        if metric in metric_ylim:
            ax.set_ylim(metric_ylim[metric])
        ax.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'boxplot_{metric}.png', dpi=150, bbox_inches='tight')
        plt.savefig(f'boxplot_{metric}.pdf', bbox_inches='tight')
        plt.close()

    # --- Scatter: frequency vs codon metrics ---
    if 'frequency' in merged_data.columns and available_metrics:
        scatter_metrics = [m for m in ['enc', 'cai', 'gc3s'] if m in available_metrics]
        if scatter_metrics:
            fig, axes = plt.subplots(1, len(scatter_metrics), figsize=(5 * len(scatter_metrics), 5))
            if len(scatter_metrics) == 1:
                axes = [axes]
            scatter_labels = {'enc': 'ENC', 'cai': 'CAI', 'gc3s': 'GC3s'}
            for i, metric in enumerate(scatter_metrics):
                ax = axes[i]
                for cat in ['rare', 'accessory', 'core']:
                    cd = merged_data[merged_data['category'] == cat]
                    if len(cd) > 0:
                        ax.scatter(cd['frequency'], cd[metric], c=category_colors[cat], alpha=0.5, s=10, label=cat.capitalize())
                valid = merged_data[['frequency', metric]].dropna()
                if len(valid) >= 10:
                    z = np.polyfit(valid['frequency'], valid[metric], 1)
                    p = np.poly1d(z)
                    x_line = np.linspace(0, 1, 100)
                    ax.plot(x_line, p(x_line), 'k--', alpha=0.5, linewidth=2)
                    rho, pval = sp_stats.spearmanr(valid['frequency'], valid[metric])
                    ax.text(0.05, 0.95, f'rho = {rho:.3f}\\np = {pval:.2e}', transform=ax.transAxes,
                            verticalalignment='top', fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                ax.set_xlabel('Gene Frequency (proportion of genomes)')
                ax.set_ylabel(scatter_labels.get(metric, metric))
                ax.legend(loc='lower right')
                ax.grid(alpha=0.3)
            plt.suptitle(f'Codon Bias vs Gene Conservation Level{title_suffix}', fontsize=14, fontweight='bold')
            plt.tight_layout()
            plt.savefig('codon_frequency_correlation.png', dpi=150, bbox_inches='tight')
            plt.savefig('codon_frequency_correlation.pdf', bbox_inches='tight')
            plt.close()

    # --- Heatmap ---
    if not stats_df.empty:
        mean_cols = [f'{m}_mean' for m in metrics if f'{m}_mean' in stats_df.columns]
        if mean_cols:
            hmap = stats_df.set_index('category')[mean_cols]
            hmap.columns = [c.replace('_mean', '').upper() for c in hmap.columns]
            hmap_norm = (hmap - hmap.mean()) / hmap.std()
            fig, ax = plt.subplots(figsize=(8, 4))
            sns.heatmap(hmap_norm, annot=hmap.round(3), fmt='', cmap='RdYlBu_r', center=0, ax=ax,
                        cbar_kws={'label': 'Z-score (relative to mean)'})
            ax.set_ylabel('Gene Category')
            ax.set_xlabel('Codon Metric')
            ax.set_title(f'Codon Bias Metrics by Gene Category{title_suffix}\\n(colors: z-score, values: actual means)')
            plt.tight_layout()
            plt.savefig('codon_category_heatmap.png', dpi=150, bbox_inches='tight')
            plt.savefig('codon_category_heatmap.pdf', bbox_inches='tight')
            plt.close()

    # =========================================================================
    # 8. Report
    # =========================================================================
    print("Generating report...")
    lines = []
    lines.append('=' * 70)
    lines.append('CODON BIAS vs PANGENOMICS CATEGORY CORRELATION ANALYSIS')
    lines.append('=' * 70)
    lines.append('')
    lines.append('ANALYSIS PARAMETERS')
    lines.append('-' * 40)
    lines.append(f'  Group: {GROUP_ID}')
    lines.append(f'  Scope: {SCOPE}')
    lines.append(f'  Samples for gene classification: {n_classification}')
    lines.append(f'  Samples with codon data: {n_codon_samples}')
    lines.append(f'  Genes with codon + category data: {len(merged_data)}')
    lines.append('')
    lines.append('SUMMARY STATISTICS BY CATEGORY')
    lines.append('-' * 40)
    if not stats_df.empty:
        for _, row in stats_df.iterrows():
            lines.append(f"\\n{row['category'].upper()} genes (n={int(row['n_genes']):,}):")
            for metric in available_metrics:
                mc = f'{metric}_mean'
                sc = f'{metric}_std'
                if mc in row:
                    lines.append(f"  {metric.upper():12s}: mean={row[mc]:.4f} +/- {row.get(sc, 0):.4f}")
    lines.append('')
    lines.append('STATISTICAL TESTS')
    lines.append('-' * 40)
    for metric, results in test_results.items():
        lines.append(f"\\n{metric.upper()}:")
        if 'kruskal_wallis' in results:
            kw = results['kruskal_wallis']
            if 'p_value' in kw:
                sig = '***' if kw['p_value'] < 0.001 else '**' if kw['p_value'] < 0.01 else '*' if kw['p_value'] < 0.05 else ''
                lines.append(f"  Kruskal-Wallis: H={kw['statistic']:.2f}, p={kw['p_value']:.2e} {sig}")
        if 'pairwise' in results:
            for pair, pw in results['pairwise'].items():
                if 'p_value' in pw:
                    sig = '*' if pw['p_value'] < 0.05 else ''
                    lines.append(f"    {pair}: p={pw['p_value']:.2e} {sig}")
        if 'frequency_correlation' in results:
            fc = results['frequency_correlation']
            if 'spearman_rho' in fc:
                sig = '*' if fc['p_value'] < 0.05 else ''
                lines.append(f"  Correlation with frequency: rho={fc['spearman_rho']:.3f}, p={fc['p_value']:.2e} {sig}")
    lines.append('')
    lines.append('INTERPRETATION GUIDE')
    lines.append('-' * 40)
    lines.append('* p < 0.05, ** p < 0.01, *** p < 0.001')
    lines.append('')
    lines.append('ENC (Effective Number of Codons):')
    lines.append('  - Range: 20-61')
    lines.append('  - Lower values = stronger codon bias')
    lines.append('  - Core genes with lower ENC suggest selection for translational efficiency')
    lines.append('')
    lines.append('CAI (Codon Adaptation Index):')
    lines.append('  - Range: 0-1')
    lines.append('  - Higher values = better adapted to reference optimal codons')
    lines.append('  - Core genes with higher CAI suggest adaptation to host expression machinery')
    lines.append('')
    lines.append('GC3s (GC at 3rd codon position):')
    lines.append('  - Reflects compositional bias')
    lines.append('  - Differences may indicate different origins (HGT vs vertical inheritance)')
    lines.append('')
    lines.append('=' * 70)

    with open('analysis_report.txt', 'w') as f:
        f.write('\\n'.join(lines))

    # Save merged data
    merged_data.to_csv('merged_data.csv', index=False)

    # =========================================================================
    # 9. Internal vs Reference comparison (if reference data available)
    # =========================================================================
    if has_reference and len(ref_codon_df) > 0:
        print("\\nProcessing reference data for comparison...")

        # Map reference genes to clusters
        ref_codon_df['cluster_id'] = ref_codon_df['gene_id'].map(locus_to_cluster)
        n_ref_mapped = ref_codon_df['cluster_id'].notna().sum()
        print(f"  Mapped {n_ref_mapped}/{len(ref_codon_df)} reference genes to clusters ({n_ref_mapped/max(len(ref_codon_df),1)*100:.1f}%)")

        ref_merged = ref_codon_df.merge(
            gene_categories.reset_index()[['gene_id', 'category', 'frequency']],
            left_on='cluster_id', right_on='gene_id',
            how='left', suffixes=('', '_cluster')
        )
        ref_merged_data = ref_merged[ref_merged['category'].isin(['core', 'accessory', 'rare'])].copy()
        print(f"  Reference categorized genes: {len(ref_merged_data)}")

        if len(ref_merged_data) > 0:
            ref_merged_data.to_csv('reference_merged_data.csv', index=False)

            # Build combined dataset with source labels
            merged_data['source'] = 'Internal'
            ref_merged_data['source'] = 'Reference'
            combined = pd.concat([
                merged_data[['category', 'source'] + available_metrics],
                ref_merged_data[['category', 'source'] + [m for m in available_metrics if m in ref_merged_data.columns]]
            ], ignore_index=True)

            ref_present_cats = [c for c in category_order if c in ref_merged_data['category'].values]
            cmp_cats = [c for c in category_order if c in present_cats and c in ref_present_cats]

            # --- Grouped comparison box plots per metric ---
            print("  Generating comparison plots...")
            source_palette = {'Internal': '#3498db', 'Reference': '#e74c3c'}
            cmp_metric_labels = {
                'enc': 'ENC (Effective Number of Codons)',
                'cai': 'CAI (Codon Adaptation Index)',
                'gc3s': 'GC3s (GC at 3rd position)',
                'gc_content': 'GC Content'
            }

            if cmp_cats:
                for metric in available_metrics:
                    if metric not in ref_merged_data.columns:
                        continue
                    fig, ax = plt.subplots(figsize=(6, 5))
                    cmp_data = combined[combined['category'].isin(cmp_cats)].dropna(subset=[metric])

                    sns.boxplot(data=cmp_data, x='category', y=metric, hue='source',
                                order=cmp_cats, hue_order=['Internal', 'Reference'],
                                palette=source_palette, ax=ax, width=0.6,
                                boxprops=dict(alpha=0.9), fliersize=0)

                    # Add p-value annotations per category
                    for k, cat in enumerate(cmp_cats):
                        int_v = merged_data.loc[merged_data['category'] == cat, metric].dropna()
                        ref_v = ref_merged_data.loc[ref_merged_data['category'] == cat, metric].dropna()
                        if len(int_v) >= 5 and len(ref_v) >= 5:
                            _, pval = sp_stats.mannwhitneyu(int_v, ref_v, alternative='two-sided')
                            sig_text = f'p={pval:.3f}' if pval >= 0.001 else 'p<0.001'
                            if pval < 0.05:
                                sig_text += '*'
                            y_max = cmp_data.loc[cmp_data['category'] == cat, metric].max()
                            ax.text(k, y_max, sig_text, ha='center', va='bottom',
                                    fontsize=8, style='italic')

                    ax.set_xlabel('')
                    ax.set_ylabel(cmp_metric_labels.get(metric, metric))
                    handles, labels = ax.get_legend_handles_labels()
                    ax.legend(handles[:2], labels[:2], fontsize=8)
                    plt.tight_layout()
                    plt.savefig(f'comparison_boxplot_{metric}.png', dpi=150, bbox_inches='tight')
                    plt.savefig(f'comparison_boxplot_{metric}.pdf', bbox_inches='tight')
                    plt.close()

                # --- Combined 2x2 comparison plot ---
                n_metrics = len([m for m in available_metrics if m in ref_merged_data.columns])
                n_cols = min(2, n_metrics)
                n_rows = (n_metrics + n_cols - 1) // n_cols
                fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
                if n_metrics == 1:
                    axes = [axes]
                else:
                    axes = axes.flatten()
                plot_idx = 0
                for metric in available_metrics:
                    if metric not in ref_merged_data.columns:
                        continue
                    ax = axes[plot_idx]
                    cmp_data = combined[combined['category'].isin(cmp_cats)].dropna(subset=[metric])

                    sns.boxplot(data=cmp_data, x='category', y=metric, hue='source',
                                order=cmp_cats, hue_order=['Internal', 'Reference'],
                                palette=source_palette, ax=ax, width=0.6,
                                boxprops=dict(alpha=0.9), fliersize=0)

                    ax.set_xlabel('')
                    ax.set_ylabel(cmp_metric_labels.get(metric, metric))
                    handles, labels = ax.get_legend_handles_labels()
                    ax.legend(handles[:2], labels[:2], fontsize=8)
                    plot_idx += 1

                for i in range(plot_idx, len(axes)):
                    axes[i].set_visible(False)
                plt.suptitle(f'Internal vs Reference Codon Bias Comparison{title_suffix}', fontsize=12, fontweight='bold')
                plt.tight_layout()
                plt.savefig('comparison_all_metrics.png', dpi=150, bbox_inches='tight')
                plt.savefig('comparison_all_metrics.pdf', bbox_inches='tight')
                plt.close()

            # --- Comparison statistics ---
            print("  Running comparison statistics...")
            comparison_results = {}
            cmp_lines = []
            cmp_lines.append('=' * 70)
            cmp_lines.append('INTERNAL vs REFERENCE CODON BIAS COMPARISON')
            cmp_lines.append('=' * 70)
            cmp_lines.append(f'Group: {GROUP_ID}')
            cmp_lines.append(f'Internal genes: {len(merged_data):,}')
            cmp_lines.append(f'Reference genes: {len(ref_merged_data):,}')
            cmp_lines.append(f'Test: Mann-Whitney U (two-sided)')
            cmp_lines.append(f'Correction: Bonferroni (per metric, {len(cmp_cats)} comparisons)')
            cmp_lines.append('')

            n_int_by_cat = merged_data['category'].value_counts()
            n_ref_by_cat = ref_merged_data['category'].value_counts()
            cmp_lines.append(f"  {'Category':<12} {'Internal':>10} {'Reference':>10}")
            cmp_lines.append(f"  {'-'*12} {'-'*10} {'-'*10}")
            for cat in categories:
                cmp_lines.append(f"  {cat:<12} {n_int_by_cat.get(cat, 0):>10,} {n_ref_by_cat.get(cat, 0):>10,}")
            cmp_lines.append('')

            for metric in available_metrics:
                if metric not in ref_merged_data.columns:
                    continue
                cmp_lines.append(f"\\n{metric_labels.get(metric, metric)}")
                cmp_lines.append('-' * 50)
                cmp_lines.append(f"  {'Category':<12} {'Int mean':>10} {'Ref mean':>10} {'Diff':>10} {'U stat':>14} {'p-value':>12} {'p(adj)':>12} {'Sig':>5}")
                comparison_results[metric] = {}

                for cat in cmp_cats:
                    int_vals = merged_data.loc[merged_data['category'] == cat, metric].dropna()
                    ref_vals = ref_merged_data.loc[ref_merged_data['category'] == cat, metric].dropna()
                    if len(int_vals) < 5 or len(ref_vals) < 5:
                        cmp_lines.append(f"  {cat:<12} {'n/a (too few observations)':>60}")
                        continue
                    int_mean = int_vals.mean()
                    ref_mean = ref_vals.mean()
                    diff = ref_mean - int_mean
                    u_stat, p_val = sp_stats.mannwhitneyu(int_vals, ref_vals, alternative='two-sided')
                    p_adj = min(p_val * len(cmp_cats), 1.0)
                    sig = '***' if p_adj < 0.001 else '**' if p_adj < 0.01 else '*' if p_adj < 0.05 else ''
                    cmp_lines.append(f"  {cat:<12} {int_mean:>10.4f} {ref_mean:>10.4f} {diff:>+10.4f} {u_stat:>14.0f} {p_val:>12.2e} {p_adj:>12.2e} {sig:>5}")
                    comparison_results[metric][cat] = {
                        'int_mean': float(int_mean), 'ref_mean': float(ref_mean),
                        'diff': float(diff), 'u_stat': float(u_stat),
                        'p_value': float(p_val), 'p_adj': float(p_adj),
                        'significant': bool(p_adj < 0.05)
                    }

            cmp_lines.append(f"\\n\\n{'='*70}")
            cmp_lines.append('INTERPRETATION')
            cmp_lines.append('=' * 70)
            cmp_lines.append('A significant difference means genes in the same pangenome category')
            cmp_lines.append('show different codon bias between internal (PAL) and reference (GCA) genomes.')
            cmp_lines.append('')
            cmp_lines.append('Large differences may indicate:')
            cmp_lines.append('  - Different selective pressures on codon usage')
            cmp_lines.append('  - Different genomic GC content between sample sets')
            cmp_lines.append('  - Compositional differences due to different ecological niches')

            with open('comparison_report.txt', 'w') as f:
                f.write('\\n'.join(cmp_lines))
            with open('comparison_statistics.json', 'w') as f:
                json.dump(comparison_results, f, indent=2, default=str)

            print("  Comparison analysis complete")

    print(f"=== Done: {GROUP_ID} / {SCOPE} ===")
    """
}
