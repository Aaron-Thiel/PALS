/*
 * Aggregate Codon QC Results
 *
 * Combines individual sample QC results into batch summary
 * for easy review and reporting.
 */

process AGGREGATE_CODON_QC {
    label 'process_single'

    publishDir "${params.outdir}/codon", mode: 'copy'

    container 'quay.io/biocontainers/pandas:2.2.1'
    conda 'conda-forge::pandas'

    input:
    path summary_files

    output:
    path "codon_qc_batch_summary.tsv", emit: batch_summary
    path "codon_qc_batch_report.json", emit: batch_report

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import json
    import glob

    files = glob.glob("*.codon_qc_summary.tsv")

    if files:
        dfs = [pd.read_csv(f, sep='\\t') for f in files]
        combined = pd.concat(dfs, ignore_index=True)
        combined = combined.sort_values('sample_id').reset_index(drop=True)
        combined.to_csv("codon_qc_batch_summary.tsv", sep='\\t', index=False)

        # Build summary report
        report = {
            'total_samples': len(combined),
            'samples_passed': int((combined['qc_status'] == 'PASS').sum()),
            'samples_warning': int((combined['qc_status'] == 'WARNING').sum()),
            'samples_failed': int((combined['qc_status'] == 'FAIL').sum()),
            'mean_outlier_pct': float(combined['outlier_pct'].mean()),
            'failed_samples': combined[combined['qc_status'] == 'FAIL']['sample_id'].tolist(),
            'warning_samples': combined[combined['qc_status'] == 'WARNING']['sample_id'].tolist(),
            'by_ref_level': {
                'species': int((combined['ref_level'] == 'species').sum()),
                'genus': int((combined['ref_level'] == 'genus').sum()),
                'none': int((combined['ref_level'] == 'none').sum())
            },
            'by_genus': combined.groupby('genus').agg({
                'sample_id': 'count',
                'qc_status': lambda x: (x == 'PASS').sum(),
                'ref_genomes_used': 'first'
            }).rename(columns={
                'sample_id': 'n_samples',
                'qc_status': 'n_passed'
            }).to_dict('index')
        }
    else:
        pd.DataFrame().to_csv("codon_qc_batch_summary.tsv", sep='\\t', index=False)
        report = {'error': 'No summary files found'}

    with open("codon_qc_batch_report.json", 'w') as f:
        json.dump(report, f, indent=2)

    print(f"Aggregated {len(files)} codon QC results")
    if files:
        passed = report['samples_passed']
        warned = report['samples_warning']
        failed = report['samples_failed']
        print(f"  PASS: {passed}, WARNING: {warned}, FAIL: {failed}")
    """
}
