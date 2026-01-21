/*
 * PASA scaffold selection filter
 * Selects the best PASA scaffold (standard vs sensitive) based on QC results
 * Uses BUSCO completeness as quality measure from the comprehensive QC JSON
 */

process PASA_FILTER {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/pasa_filter", mode: 'copy'

    input:
    tuple val(sample_id), path(pasa_standard), path(pasa_sensitive), path(qc_standard_json), path(qc_sensitive_json)

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

# Parse QC JSON files to extract BUSCO completeness scores
def get_busco_score(json_file):
    # Extract BUSCO complete_percent from QC comprehensive summary JSON
    # Handle placeholder files (NO_QC_FILE*) gracefully
    if not json_file or str(json_file).startswith('NO_QC_FILE') or not os.path.exists(json_file):
        return None
    try:
        if os.path.getsize(json_file) == 0:
            return None
        with open(json_file, 'r') as f:
            data = json.load(f)
        # The JSON contains a list of tool results
        for entry in data:
            if entry.get('tool') == 'busco':
                return entry.get('complete_percent', 0.0)
        return 0.0
    except Exception as e:
        print(f"Warning: Could not parse {json_file}: {e}")
        return None

# Check if a scaffold file is valid (exists and has content)
def is_valid_scaffold(scaffold_path):
    try:
        # Handle placeholder files (NO_FILE) gracefully
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
    # Only standard available - use it by default
    selected_scaffold = "${pasa_standard}"
    scaffold_type = "standard"
    standard_score = get_busco_score("${qc_standard_json}") or 0.0
    selection_reason = f"Standard PASA selected (only available option, BUSCO: {standard_score}%)"
    sensitive_score = None
    print(f"Only standard PASA available, selecting it by default")
elif sensitive_available and not standard_available:
    # Only sensitive available - use it by default
    selected_scaffold = "${pasa_sensitive}"
    scaffold_type = "sensitive"
    sensitive_score = get_busco_score("${qc_sensitive_json}") or 0.0
    selection_reason = f"Sensitive PASA selected (only available option, BUSCO: {sensitive_score}%)"
    standard_score = None
    print(f"Only sensitive PASA available, selecting it by default")
elif not standard_available and not sensitive_available:
    # Neither available - error
    raise RuntimeError("Neither standard nor sensitive PASA scaffold is available")
else:
    # Both available - compare BUSCO scores
    standard_score = get_busco_score("${qc_standard_json}") or 0.0
    sensitive_score = get_busco_score("${qc_sensitive_json}") or 0.0

    print(f"Standard PASA BUSCO score: {standard_score}")
    print(f"Sensitive PASA BUSCO score: {sensitive_score}")

    # Select the best scaffold based on BUSCO completeness
    if standard_score > sensitive_score:
        selected_scaffold = "${pasa_standard}"
        scaffold_type = "standard"
        selection_reason = f"Standard PASA selected (BUSCO: {standard_score}% vs {sensitive_score}%)"
    elif sensitive_score > standard_score:
        selected_scaffold = "${pasa_sensitive}"
        scaffold_type = "sensitive"
        selection_reason = f"Sensitive PASA selected (BUSCO: {sensitive_score}% vs {standard_score}%)"
    else:
        # Tie - default to sensitive
        selected_scaffold = "${pasa_sensitive}"
        scaffold_type = "sensitive"
        selection_reason = f"Sensitive PASA selected (default - equal BUSCO: {standard_score}%)"

# Copy selected scaffold
shutil.copy(selected_scaffold, "selected_scaffold.fasta")

# Create summary report
sample_id = "${sample_id}"
pasa_std = "${pasa_standard}"
pasa_sens = "${pasa_sensitive}"

summary_lines = [
    f"PASA Scaffold Selection Results for {sample_id}",
    "==============================================",
    "",
    f"Standard PASA: {pasa_std}",
    f"  Available: {standard_available}",
    f"  BUSCO Completeness: {standard_score}%" if standard_score is not None else "  BUSCO Completeness: N/A",
    "",
    f"Sensitive PASA: {pasa_sens}",
    f"  Available: {sensitive_available}",
    f"  BUSCO Completeness: {sensitive_score}%" if sensitive_score is not None else "  BUSCO Completeness: N/A",
    "",
    f"Selection Decision: {selection_reason}",
    f"Selected Scaffold: {scaffold_type}",
    "",
    "Quality Score: BUSCO completeness percentage (higher is better)"
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
    "standard_busco": standard_score,
    "sensitive_busco": sensitive_score,
    "selection_reason": selection_reason
}

with open("pasa_selection_summary.json", 'w') as f:
    json.dump(selection_data, f, indent=2)

print(f"PASA scaffold selection completed for {sample_id}")
print(f"Selected: {scaffold_type} (BUSCO: {standard_score if scaffold_type == 'standard' else sensitive_score}%)")
PYEOF
    """
}
