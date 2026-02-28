/*
 * PASA scaffold selection filter
 * Selects the best PASA scaffold (standard vs sensitive) based on QC results
 * Uses CheckM2 completeness (with heavy penalty <90%) + L90 contiguity from QUAST
 */

process PASA_FILTER {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/pasa_filter", mode: 'copy'

    input:
    tuple val(sample_id), path(pasa_standard), path(pasa_sensitive), path(checkm2_standard_json), path(quast_standard_tsv), path(checkm2_sensitive_json), path(quast_sensitive_tsv)

    output:
    tuple val(sample_id), path("selected_scaffold.fasta"), emit: selected_scaffold
    tuple val(sample_id), path("pasa_selection_summary.txt"), emit: summary
    tuple val(sample_id), path("pasa_selection_summary.json"), emit: json_reports

    script:
    """
    python3 << 'PYEOF'
import json
import shutil
import os
import csv

def get_checkm2_completeness(json_file):
    """Extract completeness from CheckM2 summary JSON."""
    if not json_file or str(json_file).startswith('NO_QC_FILE') or not os.path.exists(json_file):
        return None
    try:
        if os.path.getsize(json_file) == 0:
            return None
        with open(json_file, 'r') as f:
            data = json.load(f)
        return data.get('completeness', 0.0)
    except Exception as e:
        print(f"Warning: Could not parse CheckM2 JSON {json_file}: {e}")
        return None

def get_l90_from_quast(tsv_file):
    """Extract L90 from QUAST transposed report TSV."""
    if not tsv_file or str(tsv_file).startswith('NO_QC_FILE') or not os.path.exists(tsv_file):
        return None
    try:
        if os.path.getsize(tsv_file) == 0:
            return None
        with open(tsv_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\\t')
            for row in reader:
                val = row.get('L90')
                if val is not None:
                    return int(val)
        return None
    except Exception as e:
        print(f"Warning: Could not parse QUAST TSV {tsv_file}: {e}")
        return None

def is_valid_scaffold(scaffold_path):
    try:
        if not scaffold_path or scaffold_path == 'NO_FILE':
            return False
        if not os.path.exists(scaffold_path):
            return False
        if os.path.getsize(scaffold_path) == 0:
            return False
        return True
    except Exception:
        return False

# Check availability of each scaffold
standard_available = is_valid_scaffold("${pasa_standard}")
sensitive_available = is_valid_scaffold("${pasa_sensitive}")

print(f"Standard PASA available: {standard_available}")
print(f"Sensitive PASA available: {sensitive_available}")

# Handle cases where only one scaffold is available
if standard_available and not sensitive_available:
    selected_scaffold = "${pasa_standard}"
    scaffold_type = "standard"
    std_completeness = get_checkm2_completeness("${checkm2_standard_json}")
    std_l90 = get_l90_from_quast("${quast_standard_tsv}")
    sens_completeness = None
    sens_l90 = None
    selection_reason = f"Standard PASA selected (only available option, completeness: {std_completeness}%, L90: {std_l90})"
elif sensitive_available and not standard_available:
    selected_scaffold = "${pasa_sensitive}"
    scaffold_type = "sensitive"
    sens_completeness = get_checkm2_completeness("${checkm2_sensitive_json}")
    sens_l90 = get_l90_from_quast("${quast_sensitive_tsv}")
    std_completeness = None
    std_l90 = None
    selection_reason = f"Sensitive PASA selected (only available option, completeness: {sens_completeness}%, L90: {sens_l90})"
elif not standard_available and not sensitive_available:
    raise RuntimeError("Neither standard nor sensitive PASA scaffold is available")
else:
    # Both available - compare using CheckM2 completeness + L90
    std_completeness = get_checkm2_completeness("${checkm2_standard_json}") or 0.0
    sens_completeness = get_checkm2_completeness("${checkm2_sensitive_json}") or 0.0
    std_l90 = get_l90_from_quast("${quast_standard_tsv}")
    sens_l90 = get_l90_from_quast("${quast_sensitive_tsv}")

    print(f"Standard  - CheckM2 completeness: {std_completeness}%, L90: {std_l90}")
    print(f"Sensitive - CheckM2 completeness: {sens_completeness}%, L90: {sens_l90}")

    # Tier 1: Heavily penalize completeness < 90%
    std_above_threshold = std_completeness >= 90.0
    sens_above_threshold = sens_completeness >= 90.0

    if std_above_threshold and not sens_above_threshold:
        selected_scaffold = "${pasa_standard}"
        scaffold_type = "standard"
        selection_reason = (
            f"Standard PASA selected (completeness {std_completeness}% >= 90% vs "
            f"sensitive {sens_completeness}% < 90%)"
        )
    elif sens_above_threshold and not std_above_threshold:
        selected_scaffold = "${pasa_sensitive}"
        scaffold_type = "sensitive"
        selection_reason = (
            f"Sensitive PASA selected (completeness {sens_completeness}% >= 90% vs "
            f"standard {std_completeness}% < 90%)"
        )
    else:
        # Tier 2: Both above or both below 90% - compare L90 (lower is better)
        if std_l90 is not None and sens_l90 is not None:
            if std_l90 < sens_l90:
                selected_scaffold = "${pasa_standard}"
                scaffold_type = "standard"
                selection_reason = (
                    f"Standard PASA selected (better contiguity: L90={std_l90} vs {sens_l90}, "
                    f"completeness: {std_completeness}% vs {sens_completeness}%)"
                )
            elif sens_l90 < std_l90:
                selected_scaffold = "${pasa_sensitive}"
                scaffold_type = "sensitive"
                selection_reason = (
                    f"Sensitive PASA selected (better contiguity: L90={sens_l90} vs {std_l90}, "
                    f"completeness: {sens_completeness}% vs {std_completeness}%)"
                )
            else:
                # L90 tied - pick higher completeness
                if std_completeness >= sens_completeness:
                    selected_scaffold = "${pasa_standard}"
                    scaffold_type = "standard"
                    selection_reason = (
                        f"Standard PASA selected (tied L90={std_l90}, "
                        f"higher completeness: {std_completeness}% vs {sens_completeness}%)"
                    )
                else:
                    selected_scaffold = "${pasa_sensitive}"
                    scaffold_type = "sensitive"
                    selection_reason = (
                        f"Sensitive PASA selected (tied L90={sens_l90}, "
                        f"higher completeness: {sens_completeness}% vs {std_completeness}%)"
                    )
        else:
            # L90 unavailable for one/both - fall back to completeness
            if std_completeness >= sens_completeness:
                selected_scaffold = "${pasa_standard}"
                scaffold_type = "standard"
                selection_reason = (
                    f"Standard PASA selected (L90 unavailable, "
                    f"completeness: {std_completeness}% vs {sens_completeness}%)"
                )
            else:
                selected_scaffold = "${pasa_sensitive}"
                scaffold_type = "sensitive"
                selection_reason = (
                    f"Sensitive PASA selected (L90 unavailable, "
                    f"completeness: {sens_completeness}% vs {std_completeness}%)"
                )

# Copy selected scaffold
shutil.copy(selected_scaffold, "selected_scaffold.fasta")

# Create summary report
sample_id = "${sample_id}"

summary_lines = [
    f"PASA Scaffold Selection Results for {sample_id}",
    "==============================================",
    "",
    f"Standard PASA: ${pasa_standard}",
    f"  Available: {standard_available}",
    f"  CheckM2 Completeness: {std_completeness}%" if std_completeness is not None else "  CheckM2 Completeness: N/A",
    f"  L90: {std_l90}" if std_l90 is not None else "  L90: N/A",
    "",
    f"Sensitive PASA: ${pasa_sensitive}",
    f"  Available: {sensitive_available}",
    f"  CheckM2 Completeness: {sens_completeness}%" if sens_completeness is not None else "  CheckM2 Completeness: N/A",
    f"  L90: {sens_l90}" if sens_l90 is not None else "  L90: N/A",
    "",
    f"Selection Decision: {selection_reason}",
    f"Selected Scaffold: {scaffold_type}",
    "",
    "Scoring: Completeness >= 90% required (hard gate), then lowest L90 wins"
]

with open("pasa_selection_summary.txt", 'w') as f:
    f.write("\\n".join(summary_lines))

# Create JSON summary
selection_data = {
    "sample_id": sample_id,
    "tool": "pasa_filter",
    "selected_scaffold": scaffold_type,
    "standard_available": standard_available,
    "sensitive_available": sensitive_available,
    "standard_completeness": std_completeness,
    "sensitive_completeness": sens_completeness,
    "standard_l90": std_l90,
    "sensitive_l90": sens_l90,
    "selection_reason": selection_reason
}

with open("pasa_selection_summary.json", 'w') as f:
    json.dump(selection_data, f, indent=2)

print(f"PASA scaffold selection completed for {sample_id}")
print(f"Selected: {scaffold_type} ({selection_reason})")
PYEOF
    """
}