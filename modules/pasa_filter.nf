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
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        # The JSON contains a list of tool results
        for entry in data:
            if entry.get('tool') == 'busco':
                return entry.get('complete_percent', 0.0)
        return 0.0
    except Exception as e:
        print(f"Warning: Could not parse {json_file}: {e}")
        return 0.0

# Get scores for both scaffolds
standard_score = get_busco_score("${qc_standard_json}")
sensitive_score = get_busco_score("${qc_sensitive_json}")

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
    f"BUSCO Completeness: {standard_score}%",
    "",
    f"Sensitive PASA: {pasa_sens}",
    f"BUSCO Completeness: {sensitive_score}%",
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
