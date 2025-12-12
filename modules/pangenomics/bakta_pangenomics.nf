/*
 * Bakta annotation process for pangenomics
 * Re-annotates PASA scaffolds to generate consistent annotation files
 * Outputs: .gff3 (for Panta), .faa (for EggNOG), .fna + .gff3 (for antiSMASH)
 */

process BAKTA_PANGENOMICS {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/pangenomics/bakta", mode: 'copy'

    container 'staphb/bakta:latest'

    input:
    tuple val(sample_id), path(scaffold)

    output:
    tuple val(sample_id), path("${sample_id}.gff3"), emit: gff
    tuple val(sample_id), path("${sample_id}.faa"), emit: faa
    tuple val(sample_id), path("${sample_id}.fna"), emit: fna
    tuple val(sample_id), path("${sample_id}.ffn"), emit: ffn
    tuple val(sample_id), path("${sample_id}.gbff"), emit: genbank
    tuple val(sample_id), path("${sample_id}.tsv"), emit: tsv
    tuple val(sample_id), path("${sample_id}_bakta_pangenomics_summary.json"), emit: json_report

    script:
    """
    bakta \\
        --db /databases/bakta/db \\
        --prefix ${sample_id} \\
        --output . \\
        --threads ${task.cpus} \\
        --verbose \\
        --force \\
        ${scaffold}

    # Create JSON summary
    python3 << 'EOF'
import json
import os

def parse_bakta_tsv(tsv_file):
    data = {
        "sample_id": "${sample_id}",
        "tool": "bakta_pangenomics",
        "annotation_stats": {}
    }

    if not os.path.exists(tsv_file):
        return data

    feature_counts = {}
    total_features = 0

    try:
        with open(tsv_file, 'r') as f:
            lines = f.readlines()

        for line in lines[1:]:
            if line.strip():
                fields = line.split('\\t')
                if len(fields) > 2:
                    feature_type = fields[2]
                    feature_counts[feature_type] = feature_counts.get(feature_type, 0) + 1
                    total_features += 1

        data["annotation_stats"] = {
            "total_features": total_features,
            "feature_counts": feature_counts
        }

    except Exception as e:
        print(f"Error parsing Bakta TSV: {e}")

    return data

tsv_file = "${sample_id}.tsv"
summary_data = parse_bakta_tsv(tsv_file)

with open("${sample_id}_bakta_pangenomics_summary.json", 'w') as f:
    json.dump(summary_data, f, indent=2)

print("Bakta pangenomics JSON summary created")
EOF

    echo "Bakta annotation for pangenomics completed for ${sample_id}"
    """
}
