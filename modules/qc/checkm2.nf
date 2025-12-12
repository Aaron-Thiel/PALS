/*
 * CheckM2 quality assessment process
 * Evaluates genome completeness and contamination
 * Compares PASA sensitive vs standard scaffolds
 */

process CHECKM2 {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/qc_${stage_name}/checkm2", mode: 'copy'

    container 'staphb/checkm2:latest'

    input:
    tuple val(sample_id), path(assembly)
    val stage_name

    output:
    tuple val(sample_id), path("${sample_id}_checkm2"), emit: checkm2_results
    tuple val(sample_id), path("checkm2_summary.txt"), emit: summary
    tuple val(sample_id), path("*.json"), emit: json_reports
    
    script:
    """
    echo "Running CheckM2 quality assessment on ${assembly}..."
    
    # Set CheckM2 database path
    export CHECKM2_DB=/databases/checkm2
    
    # Create input directory structure that CheckM2 expects
    mkdir -p input_genomes
    cp ${assembly} input_genomes/${sample_id}.fasta
    
    # Run CheckM2 with explicit database file path
    checkm2 predict \\
        --threads ${task.cpus} \\
        --input input_genomes \\
        --output-directory ${sample_id}_checkm2 \\
        --extension .fasta \\
        --database_path /databases/checkm2/uniref100.KO.1.dmnd \\
        --force
    
    # Extract quality metrics
    if [ -f "${sample_id}_checkm2/quality_report.tsv" ]; then
        completeness=\$(tail -n +2 ${sample_id}_checkm2/quality_report.tsv | cut -f2 | head -1)
        contamination=\$(tail -n +2 ${sample_id}_checkm2/quality_report.tsv | cut -f3 | head -1)
        genome_size=\$(tail -n +2 ${sample_id}_checkm2/quality_report.tsv | cut -f4 | head -1)
    else
        completeness="N/A"
        contamination="N/A" 
        genome_size="N/A"
    fi
    
    # Create summary report
    cat > checkm2_summary.txt << EOF
CheckM2 Quality Assessment for ${sample_id}
==========================================

Assembly: ${assembly}
Completeness: \${completeness}%
Contamination: \${contamination}%
Genome Size: \${genome_size} bp

Quality Score: \$(awk 'BEGIN {printf "%.2f", '\$completeness' - 5 * '\$contamination'}')
EOF
    
    # Convert CheckM2 output to JSON for MultiQC
    python3 << 'EOF'
import json
import csv
import os

def parse_checkm2_tsv(file_path, sample_id):
    data = {"sample_id": sample_id}
    try:
        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    data['completeness'] = float(row.get('Completeness', 0))
                    data['contamination'] = float(row.get('Contamination', 0))
                    data['genome_size'] = int(row.get('Genome_Size', 0))
                    data['quality_score'] = data['completeness'] - 5 * data['contamination']
                    break  # Only process first row
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        data['completeness'] = 0
        data['contamination'] = 0
        data['quality_score'] = 0
    
    return data

# Parse CheckM2 results
checkm2_data = parse_checkm2_tsv("${sample_id}_checkm2/quality_report.tsv", "${sample_id}")

# Save as JSON file for MultiQC
with open("checkm2_summary.json", 'w') as f:
    json.dump(checkm2_data, f, indent=2)

print("CheckM2 JSON summary created for MultiQC")
EOF
    
    echo "CheckM2 analysis completed for ${sample_id}"
    echo "Completeness=\${completeness}%, Contamination=\${contamination}%"
    """
}