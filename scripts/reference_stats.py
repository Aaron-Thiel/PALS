#!/usr/bin/env python3
"""
Reference Selector Statistics Script
Extracts statistics from selected_references.csv files and generates visualizations.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import sys

def load_all_references(base_path: str = "results") -> pd.DataFrame:
    """Load all selected_references.csv files into a single DataFrame."""
    pattern = f"{base_path}/PAL*/classification/references/selected_references.csv"
    files = glob.glob(pattern)

    if not files:
        print(f"No files found matching pattern: {pattern}")
        sys.exit(1)

    all_data = []
    for f in files:
        # Extract sample ID from path
        sample_id = Path(f).parts[-4]  # PAL###
        try:
            df = pd.read_csv(f)
            df['sample_id'] = sample_id
            all_data.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}")

    combined = pd.concat(all_data, ignore_index=True)
    print(f"Loaded {len(files)} files with {len(combined)} total references")
    return combined


def compute_sample_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-sample statistics for ANI, AF, and quality_score."""
    stats = df.groupby('sample_id').agg(
        n_references=('accession', 'count'),
        # ANI stats
        ani_min=('ani', 'min'),
        ani_max=('ani', 'max'),
        ani_mean=('ani', 'mean'),
        ani_median=('ani', 'median'),
        ani_std=('ani', 'std'),
        # AF stats
        af_min=('af_query', 'min'),
        af_max=('af_query', 'max'),
        af_mean=('af_query', 'mean'),
        af_median=('af_query', 'median'),
        af_std=('af_query', 'std'),
        # Quality score stats
        quality_min=('quality_score', 'min'),
        quality_max=('quality_score', 'max'),
        quality_mean=('quality_score', 'mean'),
        quality_median=('quality_score', 'median'),
        quality_std=('quality_score', 'std'),
    ).reset_index()

    return stats


def print_summary_stats(df: pd.DataFrame, sample_stats: pd.DataFrame):
    """Print summary statistics to console."""
    print("\n" + "="*80)
    print("REFERENCE SELECTOR STATISTICS SUMMARY")
    print("="*80)

    print(f"\nTotal samples analyzed: {sample_stats['sample_id'].nunique()}")
    print(f"Total references selected: {len(df)}")

    print("\n--- References per Sample ---")
    print(f"  Min:    {sample_stats['n_references'].min()}")
    print(f"  Max:    {sample_stats['n_references'].max()}")
    print(f"  Mean:   {sample_stats['n_references'].mean():.1f}")
    print(f"  Median: {sample_stats['n_references'].median():.1f}")

    metrics = [("ANI", "ani"), ("AF (Alignment Fraction)", "af_query"), ("Quality Score", "quality_score")]
    # Add composite_score if it exists in the dataframe
    if "composite_score" in df.columns:
        metrics.append(("Composite Score", "composite_score"))

    for metric, col in metrics:
        if col in df.columns:
            print(f"\n--- {metric} (across all references) ---")
            print(f"  Min:    {df[col].min():.2f}")
            print(f"  Max:    {df[col].max():.2f}")
            print(f"  Mean:   {df[col].mean():.2f}")
            print(f"  Median: {df[col].median():.2f}")
            print(f"  Std:    {df[col].std():.2f}")

    # Selection reason breakdown
    print("\n--- Selection Reason Breakdown ---")
    reason_counts = df['selection_reason'].value_counts()
    for reason, count in reason_counts.items():
        pct = count / len(df) * 100
        print(f"  {reason}: {count} ({pct:.1f}%)")


def create_histograms(df: pd.DataFrame, output_dir: Path):
    """Create histograms for ANI, AF, and quality_score."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    metrics = [
        ('ani', 'ANI (%)', 'steelblue'),
        ('af_query', 'Alignment Fraction (%)', 'darkorange'),
        ('quality_score', 'Quality Score', 'forestgreen')
    ]

    for ax, (col, label, color) in zip(axes, metrics):
        data = df[col].dropna()
        ax.hist(data, bins=30, color=color, edgecolor='black', alpha=0.7)
        ax.set_xlabel(label, fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_title(f'Distribution of {label}', fontsize=14)
        ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {data.mean():.2f}')
        ax.axvline(data.median(), color='darkred', linestyle=':', linewidth=2, label=f'Median: {data.median():.2f}')
        ax.legend()

    plt.tight_layout()
    outfile = output_dir / 'reference_histograms.png'
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved histogram: {outfile}")


def create_boxplots(df: pd.DataFrame, output_dir: Path):
    """Create boxplots for ANI, AF, and quality_score."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 6))

    metrics = [
        ('ani', 'ANI (%)', 'steelblue'),
        ('af_query', 'Alignment Fraction (%)', 'darkorange'),
        ('quality_score', 'Quality Score', 'forestgreen')
    ]

    for ax, (col, label, color) in zip(axes, metrics):
        data = df[col].dropna()
        bp = ax.boxplot(data, patch_artist=True)
        bp['boxes'][0].set_facecolor(color)
        bp['boxes'][0].set_alpha(0.7)
        ax.set_ylabel(label, fontsize=12)
        ax.set_title(f'{label} Distribution', fontsize=14)
        ax.set_xticklabels(['All References'])

        # Add stats annotation
        stats_text = f'n={len(data)}\nMean={data.mean():.2f}\nMedian={data.median():.2f}'
        ax.text(1.15, data.median(), stats_text, fontsize=10, verticalalignment='center')

    plt.tight_layout()
    outfile = output_dir / 'reference_boxplots.png'
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved boxplot: {outfile}")


def create_per_sample_boxplots(df: pd.DataFrame, output_dir: Path):
    """Create per-sample boxplots for each metric."""
    # Sort samples by median ANI
    sample_order = df.groupby('sample_id')['ani'].median().sort_values(ascending=False).index.tolist()

    fig, axes = plt.subplots(3, 1, figsize=(16, 14))

    metrics = [
        ('ani', 'ANI (%)', 'steelblue'),
        ('af_query', 'Alignment Fraction (%)', 'darkorange'),
        ('quality_score', 'Quality Score', 'forestgreen')
    ]

    for ax, (col, label, color) in zip(axes, metrics):
        data_by_sample = [df[df['sample_id'] == s][col].dropna().values for s in sample_order]
        bp = ax.boxplot(data_by_sample, patch_artist=True)
        for box in bp['boxes']:
            box.set_facecolor(color)
            box.set_alpha(0.7)
        ax.set_xticklabels(sample_order, rotation=90, fontsize=8)
        ax.set_ylabel(label, fontsize=12)
        ax.set_title(f'{label} by Sample', fontsize=14)
        ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    outfile = output_dir / 'reference_boxplots_per_sample.png'
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved per-sample boxplot: {outfile}")


def create_selection_reason_plot(df: pd.DataFrame, output_dir: Path):
    """Create a bar chart showing selection reason distribution."""
    fig, ax = plt.subplots(figsize=(10, 6))

    reason_counts = df['selection_reason'].value_counts()
    colors = plt.cm.Set2(np.linspace(0, 1, len(reason_counts)))

    bars = ax.bar(reason_counts.index, reason_counts.values, color=colors, edgecolor='black')
    ax.set_xlabel('Selection Reason', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Reference Selection Reasons', fontsize=14)
    plt.xticks(rotation=45, ha='right')

    # Add value labels on bars
    for bar, val in zip(bars, reason_counts.values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(val), ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    outfile = output_dir / 'reference_selection_reasons.png'
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved selection reasons plot: {outfile}")


def main():
    # Setup paths
    base_path = Path("results")
    output_dir = Path("results/reference_stats")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print("Loading reference data...")
    df = load_all_references(str(base_path))

    # Compute per-sample statistics
    sample_stats = compute_sample_stats(df)

    # Save sample stats to CSV
    stats_file = output_dir / 'sample_statistics.csv'
    sample_stats.to_csv(stats_file, index=False)
    print(f"Saved sample statistics: {stats_file}")

    # Print summary
    print_summary_stats(df, sample_stats)

    # Create visualizations
    print("\nGenerating visualizations...")
    create_histograms(df, output_dir)
    create_boxplots(df, output_dir)
    create_per_sample_boxplots(df, output_dir)
    create_selection_reason_plot(df, output_dir)

    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
