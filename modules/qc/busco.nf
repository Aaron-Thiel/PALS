/*
 * BUSCO QC module - general quality control assessment
 * Can analyze one or more assemblies for completeness
 */

process BUSCO {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/qc_${stage_name}/busco", mode: 'copy'
    
    container 'staphb/busco:latest'
    
    input:
    tuple val(sample_id), path(assemblies)
    val stage_name
    
    output:
    tuple val(sample_id), path("busco_results"), emit: results
    tuple val(sample_id), path("busco_summary.txt"), emit: summary
    tuple val(sample_id), path("*.json"), emit: json_reports
    
    script:
    def assembly_list = assemblies instanceof List ? assemblies : [assemblies]
    def _assembly_count = assembly_list.size()
    
    """
    # Set up BUSCO database path
    export BUSCO_CONFIG_FILE=busco_config.ini
    
    # Create config file for BUSCO
    cat > busco_config.ini << EOF
[busco]
offline = True
download_path = /databases/busco
EOF
    
    mkdir -p busco_results
    
    # Process each assembly
    assembly_index=1
    for assembly in ${assembly_list.join(' ')}; do
        if [ "\$assembly" != "NO_FILE" ] && [ -f "\$assembly" ]; then
            echo "Running BUSCO on \$assembly (assembly \$assembly_index)..."
            busco \\
                -i \$assembly \\
                -o busco_results/assembly_\${assembly_index} \\
                -l /databases/busco/lineages/lactobacillaceae_odb12 \\
                -m genome \\
                -c ${task.cpus} \\
                -f
        else
            echo "Skipping invalid assembly: \$assembly"
        fi
        assembly_index=\$((assembly_index + 1))
    done
    
    # Create summary report
    cat > busco_summary.txt << EOF
BUSCO QC Analysis Report for ${sample_id}
Stage: ${stage_name}
=========================================

EOF
    
    # Extract completeness for each assembly
    assembly_index=1
    for assembly in ${assembly_list.join(' ')}; do
        if [ "\$assembly" != "NO_FILE" ] && [ -f "\$assembly" ] && [ -d "busco_results/assembly_\${assembly_index}" ]; then
            summary_file=\$(find busco_results/assembly_\${assembly_index} -name "short_summary*.txt" | head -1)
            if [ -f "\$summary_file" ]; then
                completeness=\$(grep -E "^\\s*C:" "\$summary_file" | sed 's/.*C:\\([0-9.]*\\)%.*/\\1/')
                single_copy=\$(grep -E "^\\s*S:" "\$summary_file" | sed 's/.*S:\\([0-9.]*\\)%.*/\\1/')
                duplicated=\$(grep -E "^\\s*D:" "\$summary_file" | sed 's/.*D:\\([0-9.]*\\)%.*/\\1/')
                fragmented=\$(grep -E "^\\s*F:" "\$summary_file" | sed 's/.*F:\\([0-9.]*\\)%.*/\\1/')
                missing=\$(grep -E "^\\s*M:" "\$summary_file" | sed 's/.*M:\\([0-9.]*\\)%.*/\\1/')
                
                cat >> busco_summary.txt << INNER_EOF
Assembly \${assembly_index} (\$assembly):
  Complete BUSCOs: \${completeness}%
  Single-copy: \${single_copy}%
  Duplicated: \${duplicated}%
  Fragmented: \${fragmented}%
  Missing: \${missing}%

INNER_EOF
            fi
        fi
        assembly_index=\$((assembly_index + 1))
    done
    
    echo "Analysis based on lactobacillaceae_odb12 lineage" >> busco_summary.txt
    
    # Convert to JSON for MultiQC
    python3 << 'EOF'
import json
import re
import os
import glob

def parse_busco_summary(file_path, assembly_name):
    data = {"assembly_name": assembly_name, "stage": "${stage_name}"}
    try:
        with open(file_path, 'r') as f:
            content = f.read()
            
        # Extract percentages
        complete_match = re.search(r'C:([0-9.]+)%', content)
        single_match = re.search(r'S:([0-9.]+)%', content)
        duplicated_match = re.search(r'D:([0-9.]+)%', content)
        fragmented_match = re.search(r'F:([0-9.]+)%', content)
        missing_match = re.search(r'M:([0-9.]+)%', content)
        
        if complete_match:
            data['complete_percent'] = float(complete_match.group(1))
        if single_match:
            data['single_copy_percent'] = float(single_match.group(1))
        if duplicated_match:
            data['duplicated_percent'] = float(duplicated_match.group(1))
        if fragmented_match:
            data['fragmented_percent'] = float(fragmented_match.group(1))
        if missing_match:
            data['missing_percent'] = float(missing_match.group(1))
            
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
    
    return data

# Find all summary files and process them
summary_files = glob.glob("busco_results/*/short_summary*.txt")
all_data = []

for i, summary_file in enumerate(summary_files, 1):
    assembly_name = f"assembly_{i}"
    data = parse_busco_summary(summary_file, assembly_name)
    all_data.append(data)

# Save as JSON file
with open("busco_qc_summary.json", 'w') as f:
    json.dump(all_data, f, indent=2)

print("BUSCO QC JSON summary created for MultiQC")
EOF
    
    echo "BUSCO QC analysis completed for ${sample_id} at stage ${stage_name}"
    """
}