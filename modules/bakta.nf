/*
 * Bakta annotation process
 * Annotates bacterial genomes with Bakta
 */

process BAKTA {
    tag "$sample_id"

    publishDir "${params.outdir}/bakta/${sample_id}", mode: 'copy'

    container 'staphb/bakta:latest'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}.gff3"), emit: gff
    tuple val(sample_id), path("${sample_id}.gbff"), emit: genbank
    tuple val(sample_id), path("${sample_id}.faa"), emit: proteins
    tuple val(sample_id), path("${sample_id}.ffn"), emit: nucleotides
    tuple val(sample_id), path("${sample_id}.tsv"), emit: tsv
    tuple val(sample_id), path("${sample_id}.fna"), emit: fna
    tuple val(sample_id), path("${sample_id}_bakta_summary.json"), emit: json_report
    tuple val(sample_id), path("${sample_id}.txt"), emit: summary, optional: true
    
    script:
    """
    bakta \\
        --db /databases/bakta/db \\
        --prefix ${sample_id} \\
        --output . \\
        --threads ${task.cpus} \\
        --verbose \\
        --force \\
        ${contigs}
    
    # Create JSON summary for MultiQC
    python3 << 'EOF'
import json
import os

def parse_bakta_tsv(tsv_file):
    data = {
        "sample_id": "${sample_id}",
        "tool": "bakta",
        "annotation_stats": {}
    }

    if not os.path.exists(tsv_file):
        return data

    # Count features by type
    # Bakta TSV columns: Sequence Id, Type, Start, Stop, Strand, Locus Tag, Gene, Product, DbXrefs
    feature_counts = {}
    total_features = 0

    try:
        with open(tsv_file, 'r') as f:
            for line in f:
                # Skip comment and header lines (start with #)
                if line.startswith('#'):
                    continue
                if not line.strip():
                    continue

                fields = line.strip().split('\\t')
                if len(fields) > 1:
                    feature_type = fields[1]  # Type is column index 1
                    feature_counts[feature_type] = feature_counts.get(feature_type, 0) + 1
                    total_features += 1

        data["annotation_stats"] = {
            "total_features": total_features,
            "feature_counts": feature_counts
        }

    except Exception as e:
        print(f"Error parsing Bakta TSV: {e}")

    return data

# Parse the TSV file
tsv_file = "${sample_id}.tsv"
summary_data = parse_bakta_tsv(tsv_file)

# Save as JSON
with open("${sample_id}_bakta_summary.json", 'w') as f:
    json.dump(summary_data, f, indent=2)

print("Bakta JSON summary created for MultiQC")
EOF
    
    echo "Bakta annotation completed for ${sample_id}"
    """
}
