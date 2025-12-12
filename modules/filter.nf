/*
 * Comprehensive QC filtering process
 * Filters genomes based on multiple criteria:
 * - GTDB-Tk taxonomic classification (Lactobacillaceae)
 * - CheckM2 contamination threshold (<5%)
 * - BUSCO completeness threshold (>90%)
 * All criteria are optional and user-configurable
 */

process FILTER {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/filter", mode: 'copy'
    
    input:
    tuple val(sample_id), path(genome_fasta), path(gtdbtk_classification), path(checkm2_results), path(busco_results)
    
    output:
    tuple val(sample_id), path("filter_summary.txt"), emit: summary
    tuple val(sample_id), path("*.json"), emit: json_reports
    tuple val(sample_id), path("filtered_genome.fasta"), emit: filtered_genome, optional: true
    tuple val(sample_id), path("passed_filter.flag"), emit: passed_filter, optional: true
    
    script:
    """
    # Initialize filter variables
    filter_passed=true
    filter_reasons=()
    
    # Read GTDB-Tk classification if provided
    if [ -f "${gtdbtk_classification}" ] && [ "${params.filter_family}" = "true" ]; then
        classification=\$(cat ${gtdbtk_classification})
        if [ "\$classification" = "FAILED" ]; then
            filter_passed=false
            filter_reasons+=("GTDB-Tk classification failed")
        elif ! echo "\$classification" | grep -q "${params.filter_target_family}"; then
            filter_passed=false
            filter_reasons+=("Not ${params.filter_target_family} family (found: \$classification)")
        else
            filter_reasons+=("✓ ${params.filter_target_family} family confirmed")
        fi
    else
        classification="Not analyzed"
        if [ "${params.filter_family}" = "true" ]; then
            filter_reasons+=("GTDB-Tk classification file missing")
        fi
    fi
    
    # Extract CheckM2 results for both contamination and completeness
    contamination="Not analyzed"
    completeness="Not analyzed"
    
    # Check if CheckM2 JSON results are provided
    if [ -f "${checkm2_results}" ] && [ "${checkm2_results}" != "no_checkm2" ] && [ "${checkm2_results}" != "*/empty_checkm2.txt" ]; then
        echo "Processing CheckM2 results from: ${checkm2_results}"
        
        # Parse CheckM2 JSON results
        python3 << 'CHECKM2_EOF'
import json
import sys

try:
    with open("${checkm2_results}", 'r') as f:
        data = json.load(f)
    
    # Extract contamination and completeness from CheckM2 JSON
    contamination_val = data.get('contamination', 'N/A')
    completeness_val = data.get('completeness', 'N/A')
    
    # Write results to temporary files
    with open('checkm2_contamination.tmp', 'w') as f:
        f.write(str(contamination_val))
    with open('checkm2_completeness.tmp', 'w') as f:
        f.write(str(completeness_val))
    
    print(f"Extracted CheckM2: contamination={contamination_val}, completeness={completeness_val}")
    
except Exception as e:
    print(f"Error parsing CheckM2 JSON: {e}")
    with open('checkm2_contamination.tmp', 'w') as f:
        f.write('N/A')
    with open('checkm2_completeness.tmp', 'w') as f:
        f.write('N/A')

CHECKM2_EOF
        
        # Read parsed values
        if [ -f "checkm2_contamination.tmp" ]; then
            contamination=\$(cat checkm2_contamination.tmp)
        fi
        if [ -f "checkm2_completeness.tmp" ]; then
            completeness=\$(cat checkm2_completeness.tmp)
        fi
        
        # Clean up temp files
        rm -f checkm2_contamination.tmp checkm2_completeness.tmp
    fi
    
    # Apply contamination filter if enabled
    if [ "${params.filter_contamination}" = "true" ]; then
        if [ "\$contamination" != "Not analyzed" ] && [ "\$contamination" != "N/A" ]; then
            # Check if contamination is above threshold
            if awk "BEGIN {exit !(\$contamination > ${params.filter_max_contamination})}"; then
                filter_passed=false
                filter_reasons+=("Contamination too high: \${contamination}% > ${params.filter_max_contamination}%")
            else
                filter_reasons+=("✓ Contamination acceptable: \${contamination}% ≤ ${params.filter_max_contamination}%")
            fi
        else
            # Don't fail if CheckM2 is completely unavailable - just warn
            if [ "${checkm2_results}" = "no_checkm2" ]; then
                filter_reasons+=("⚠ CheckM2 results unavailable - contamination filter skipped")
            else
                filter_passed=false
                filter_reasons+=("CheckM2 contamination value missing or invalid")
            fi
        fi
    else
        filter_reasons+=("Contamination filter disabled")
    fi
    
    # Apply completeness filter if enabled (using CheckM2 instead of BUSCO)
    if [ "${params.filter_completeness}" = "true" ]; then
        if [ "\$completeness" != "Not analyzed" ] && [ "\$completeness" != "N/A" ]; then
            # Check if completeness is below threshold
            if awk "BEGIN {exit !(\$completeness < ${params.filter_min_completeness})}"; then
                filter_passed=false
                filter_reasons+=("Completeness too low: \${completeness}% < ${params.filter_min_completeness}%")
            else
                filter_reasons+=("✓ Completeness acceptable: \${completeness}% ≥ ${params.filter_min_completeness}%")
            fi
        else
            # Don't fail if CheckM2 is completely unavailable - just warn
            if [ "${checkm2_results}" = "no_checkm2" ]; then
                filter_reasons+=("⚠ CheckM2 results unavailable - completeness filter skipped")
            else
                filter_passed=false
                filter_reasons+=("CheckM2 completeness value missing or invalid")
            fi
        fi
    else
        filter_reasons+=("Completeness filter disabled")
    fi
    
    # Note: Completeness filtering is now handled above using CheckM2 data
    # BUSCO results are no longer used for filtering, only CheckM2
    
    # Create output files if all filters passed
    if [ "\$filter_passed" = "true" ]; then
        cp ${genome_fasta} filtered_genome.fasta
        touch passed_filter.flag
        echo "✅ PASS: Sample ${sample_id} passed all quality filters"
    else
        echo "❌ FILTERED: Sample ${sample_id} failed quality filters"
    fi
    
    # Create detailed summary report
    cat > filter_summary.txt << EOF
Comprehensive QC Filter Results for ${sample_id}
===============================================

Genome: ${genome_fasta}

Filter Criteria Applied:
- Family filter: ${params.filter_family} (target: ${params.filter_target_family})
- CheckM2 contamination filter: ${params.filter_contamination} (max: ${params.filter_max_contamination}%)
- CheckM2 completeness filter: ${params.filter_completeness} (min: ${params.filter_min_completeness}%)

Results:
- GTDB-Tk Classification: \$classification
- CheckM2 Contamination: \$contamination
- CheckM2 Completeness: \$completeness

Filter Status: \$(if [ "\$filter_passed" = "true" ]; then echo "PASS"; else echo "FILTERED"; fi)

Filter Details:
\$(printf '%s\\n' "\${filter_reasons[@]}")

Next Step: \$(if [ "\$filter_passed" = "true" ]; then echo "Proceed to annotation (Bakta)"; else echo "Excluded from downstream analysis"; fi)
EOF
    
    # Convert filtering results to JSON for MultiQC
    python3 << 'EOF'
import json
import os
import glob

# Create filtering summary data
filter_data = {
    "sample_id": "${sample_id}",
    "tool": "comprehensive_filter",
    "filter_passed": "\$filter_passed" == "true",
    "filter_criteria": {
        "family_filter": "${params.filter_family}" == "true",
        "contamination_filter": "${params.filter_contamination}" == "true", 
        "completeness_filter": "${params.filter_completeness}" == "true"
    },
    "filter_thresholds": {
        "target_family": "${params.filter_target_family}",
        "max_contamination": ${params.filter_max_contamination},
        "min_completeness": ${params.filter_min_completeness}
    },
    "filter_status": "PASS" if "\$filter_passed" == "true" else "FILTERED"
}

# Read actual values from files
try:
    # GTDB-Tk classification
    if os.path.exists("${gtdbtk_classification}"):
        with open("${gtdbtk_classification}", 'r') as f:
            classification = f.read().strip()
            filter_data["gtdbtk_classification"] = classification
            filter_data["is_target_family"] = "${params.filter_target_family}" in classification
    
    # CheckM2 contamination
    checkm2_files = glob.glob("${checkm2_results}/**/quality_report.tsv", recursive=True)
    if checkm2_files:
        with open(checkm2_files[0], 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:
                contamination = lines[1].split('\t')[2]
                filter_data["checkm2_contamination"] = float(contamination) if contamination.replace('.', '').isdigit() else None
    
    # BUSCO completeness  
    busco_files = glob.glob("${busco_results}/**/short_summary*.txt", recursive=True)
    if busco_files:
        with open(busco_files[0], 'r') as f:
            content = f.read()
            import re
            match = re.search(r'C:([0-9.]+)%', content)
            if match:
                filter_data["busco_completeness"] = float(match.group(1))

except Exception as e:
    print(f"Error parsing QC files: {e}")

# Save as JSON file for MultiQC
with open("comprehensive_filter_summary.json", 'w') as f:
    json.dump(filter_data, f, indent=2)

print("Comprehensive filter JSON summary created for MultiQC")
EOF
    
    echo "Comprehensive QC filtering completed for ${sample_id}"
    echo "Filter result: \$(if [ "\$filter_passed" = "true" ]; then echo "PASS"; else echo "FILTERED"; fi)"
    """
}