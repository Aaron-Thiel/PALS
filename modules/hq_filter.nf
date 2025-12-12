/*
 * HQ_FILTER - High Quality Sample Filter
 * Filters samples based on CheckM2 quality metrics before classification
 * Only samples meeting quality thresholds proceed to taxonomic classification
 */

process HQ_FILTER {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/hq_filter", mode: 'copy'

    input:
    tuple val(sample_id), path(genome_fasta), path(checkm2_json)

    output:
    tuple val(sample_id), path("filtered_genome.fasta"), emit: hq_genome, optional: true
    tuple val(sample_id), path("hq_filter_summary.txt"), emit: summary
    tuple val(sample_id), path("hq_filter_result.json"), emit: json_report

    script:
    """
    #!/usr/bin/env python3
    import json
    import shutil
    import os

    # Read CheckM2 results
    checkm2_file = "${checkm2_json}"
    sample_id = "${sample_id}"
    genome_fasta = "${genome_fasta}"

    # Quality thresholds from params
    min_completeness = ${params.filter_min_completeness}
    max_contamination = ${params.filter_max_contamination}

    # Parse CheckM2 JSON
    completeness = None
    contamination = None
    try:
        with open(checkm2_file, 'r') as f:
            data = json.load(f)
            completeness = data.get('completeness')
            contamination = data.get('contamination')
    except Exception as e:
        print(f"Error reading CheckM2 JSON: {e}")

    # Determine if sample is high quality
    is_hq = False
    reasons = []

    if completeness is None or contamination is None:
        reasons.append("CheckM2 results unavailable or invalid")
    else:
        # Check completeness
        if completeness >= min_completeness:
            reasons.append(f"Completeness OK: {completeness:.1f}% >= {min_completeness}%")
            completeness_ok = True
        else:
            reasons.append(f"Completeness FAIL: {completeness:.1f}% < {min_completeness}%")
            completeness_ok = False

        # Check contamination
        if contamination <= max_contamination:
            reasons.append(f"Contamination OK: {contamination:.1f}% <= {max_contamination}%")
            contamination_ok = True
        else:
            reasons.append(f"Contamination FAIL: {contamination:.1f}% > {max_contamination}%")
            contamination_ok = False

        is_hq = completeness_ok and contamination_ok

    # Create output files
    result = {
        "sample_id": sample_id,
        "is_high_quality": is_hq,
        "completeness": completeness,
        "contamination": contamination,
        "min_completeness_threshold": min_completeness,
        "max_contamination_threshold": max_contamination,
        "reasons": reasons
    }

    with open("hq_filter_result.json", 'w') as f:
        json.dump(result, f, indent=2)

    # Create summary
    status = "PASS - HIGH QUALITY" if is_hq else "FAIL - LOW QUALITY"
    with open("hq_filter_summary.txt", 'w') as f:
        f.write(f"HQ Filter Results for {sample_id}\\n")
        f.write("=" * 50 + "\\n\\n")
        f.write(f"Status: {status}\\n\\n")
        f.write(f"Quality Metrics:\\n")
        f.write(f"  Completeness: {completeness}%\\n")
        f.write(f"  Contamination: {contamination}%\\n\\n")
        f.write(f"Thresholds:\\n")
        f.write(f"  Min completeness: {min_completeness}%\\n")
        f.write(f"  Max contamination: {max_contamination}%\\n\\n")
        f.write(f"Details:\\n")
        for reason in reasons:
            f.write(f"  - {reason}\\n")
        f.write(f"\\nNext step: {'Proceed to classification' if is_hq else 'Excluded from downstream analysis'}\\n")

    # Copy genome only if high quality
    if is_hq:
        shutil.copy(genome_fasta, "filtered_genome.fasta")
        print(f"Sample {sample_id}: HIGH QUALITY - proceeding to classification")
    else:
        print(f"Sample {sample_id}: LOW QUALITY - excluded from classification")
    """
}
