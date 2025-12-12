/*
 * Sourmash genus-level taxonomic classification process
 * Classifies bacterial genomes using sourmash with GTDB RS226 database
 * Extracts genus-level classification for downstream genome downloads
 */

process SOURMASH {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/classification/sourmash", mode: 'copy'
    
    container 'nanozoo/sourmash:4.8.14--4df0447'
    
    input:
    tuple val(sample_id), path(genome_fasta)
    
    output:
    tuple val(sample_id), path("sourmash_output"), emit: sourmash_results
    tuple val(sample_id), path("sourmash_summary.txt"), emit: summary
    tuple val(sample_id), path("*.json"), emit: json_reports
    tuple val(sample_id), path("classification.txt"), emit: classification
    tuple val(sample_id), path("genus.txt"), emit: genus
    tuple val(sample_id), path("family.txt"), emit: family
    
    script:
    """
    # Create sourmash signature from input genome (k=31 to match database)
    sourmash sketch dna \\
        -p k=31,scaled=1000 \\
        --name ${sample_id} \\
        -o ${sample_id}.sig \\
        ${genome_fasta}
    
    # Set up database paths (databases should be pre-installed)
    SOURMASH_DB_DIR="/databases/sourmash"
    
    # Create output directory
    mkdir -p sourmash_output
    
    # Run sourmash gather for classification using k=31 database
    # Note: --output-dir not available in sourmash 2.3.0, using -o only
    sourmash gather \\
        ${sample_id}.sig \\
        "\$SOURMASH_DB_DIR/gtdb-reps-rs226-k31.dna.zip" \\
        --threshold-bp 5000 \\
        -o sourmash_gather.csv
    
    # Skip taxonomy command as it may not be available in this version
    # We'll extract taxonomy directly from the GTDB database names
    
    # Extract classification results
    if [ -f "sourmash_gather.csv" ] && [ -s "sourmash_gather.csv" ]; then
        # Get best match (first line after header)  
        best_match=\$(tail -n +2 sourmash_gather.csv | head -1)
        
        if [ -n "\$best_match" ]; then
            # Extract taxon information from taxonomy file using best match name
            # Get full match name and extract just the accession ID (before first space)
            full_match_name=\$(echo "\$best_match" | cut -d',' -f10 | tr -d '"')
            match_id=\$(echo "\$full_match_name" | cut -d' ' -f1)
            
            echo "Full match name: \$full_match_name"
            echo "Extracted ID: \$match_id"
            
            # Get taxonomic information from lineages file (CSV format)
            tax_line=\$(grep "^\$match_id," "\$SOURMASH_DB_DIR/gtdb-rs226.lineages.csv" | head -1)
            
            if [ -n "\$tax_line" ]; then
                # Extract genus and family from CSV columns (f5=family, f6=genus)
                genus=\$(echo "\$tax_line" | cut -d',' -f7 | tr -d '"' | sed 's/g__//')
                family=\$(echo "\$tax_line" | cut -d',' -f6 | tr -d '"' | sed 's/f__//')
                
                # Reconstruct full taxonomy lineage for display
                superkingdom=\$(echo "\$tax_line" | cut -d',' -f2 | tr -d '"')
                phylum=\$(echo "\$tax_line" | cut -d',' -f3 | tr -d '"')
                class=\$(echo "\$tax_line" | cut -d',' -f4 | tr -d '"')
                order=\$(echo "\$tax_line" | cut -d',' -f5 | tr -d '"')
                species=\$(echo "\$tax_line" | cut -d',' -f8 | tr -d '"')
                full_taxonomy="\$superkingdom;\$phylum;\$class;\$order;\$family;\$genus;\$species"
                
                # Clean up empty or malformed entries
                if [[ "\$genus" == "" ]] || [[ "\$genus" == "g__" ]]; then
                    genus="Unclassified"
                fi
                if [[ "\$family" == "" ]] || [[ "\$family" == "f__" ]]; then
                    family="Unclassified"
                fi
                
                echo "✓ Sample ${sample_id} classified: \$full_taxonomy"
                echo "  Genus: \$genus"
                echo "  Family: \$family"
                
                # Save classification outputs
                echo "\$full_taxonomy" > classification.txt
                echo "\$genus" > genus.txt
                echo "\$family" > family.txt
                
                classification="\$full_taxonomy"
                classification_status="SUCCESS"
            else
                echo "✗ Sample ${sample_id}: No taxonomy found for best match"
                classification="No taxonomy found"
                classification_status="FAILED"
                echo "FAILED" > classification.txt
                echo "Unclassified" > genus.txt
                echo "Unclassified" > family.txt
            fi
        else
            echo "✗ Sample ${sample_id}: No significant matches found"
            classification="No matches found"
            classification_status="FAILED"
            echo "FAILED" > classification.txt
            echo "Unclassified" > genus.txt
            echo "Unclassified" > family.txt
        fi
    else
        echo "✗ Sample ${sample_id}: Sourmash gather failed or produced no output"
        classification="Sourmash failed"
        classification_status="FAILED"
        echo "FAILED" > classification.txt
        echo "Unclassified" > genus.txt
        echo "Unclassified" > family.txt
    fi
    
    # Create summary report
    cat > sourmash_summary.txt << EOF
Sourmash Taxonomic Classification for ${sample_id}
==================================================

Genome: ${genome_fasta}
Classification: \${classification}
Genus: \$(cat genus.txt 2>/dev/null || echo "Unknown")
Family: \$(cat family.txt 2>/dev/null || echo "Unknown")
Status: \${classification_status}

Database: GTDB RS226 Representatives
k-mer size: 31
Method: sourmash gather + taxonomy
Note: Memory-efficient alternative to GTDB-Tk for genus-level classification
EOF
    
    # Convert sourmash output to JSON for MultiQC
    python3 << 'EOF'
import json
import csv
import os

def parse_sourmash_results(sample_id):
    data = {"sample_id": sample_id, "tool": "sourmash_classify"}
    
    try:
        # Read genus and family
        if os.path.exists("genus.txt"):
            with open("genus.txt", 'r') as f:
                data['genus'] = f.read().strip()
        else:
            data['genus'] = 'Unknown'
            
        if os.path.exists("family.txt"):
            with open("family.txt", 'r') as f:
                data['family'] = f.read().strip()
        else:
            data['family'] = 'Unknown'
            
        if os.path.exists("classification.txt"):
            with open("classification.txt", 'r') as f:
                data['classification'] = f.read().strip()
        else:
            data['classification'] = 'Unknown'
        
        # Parse gather results if available
        gather_file = "sourmash_gather.csv"
        if os.path.exists(gather_file):
            with open(gather_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    data['intersect_bp'] = int(row.get('intersect_bp', 0))
                    data['f_match'] = float(row.get('f_match', 0))
                    data['f_unique_to_query'] = float(row.get('f_unique_to_query', 0))
                    data['f_unique_weighted'] = float(row.get('f_unique_weighted', 0))
                    break  # Only take best match
        else:
            data['intersect_bp'] = 0
            data['f_match'] = 0
            data['f_unique_to_query'] = 0
            data['f_unique_weighted'] = 0
            
    except Exception as e:
        print(f"Error parsing sourmash results: {e}")
        data['genus'] = 'Error'
        data['family'] = 'Error'
        data['classification'] = 'Error'
        data['intersect_bp'] = 0
        data['f_match'] = 0
        data['f_unique_to_query'] = 0
        data['f_unique_weighted'] = 0
    
    return data

# Parse sourmash results
sourmash_data = parse_sourmash_results("${sample_id}")

# Save as JSON file for MultiQC
with open("sourmash_classify_summary.json", 'w') as f:
    json.dump(sourmash_data, f, indent=2)

print("Sourmash classification JSON summary created for MultiQC")
EOF
    
    echo "Sourmash classification completed for ${sample_id}"
    echo "Classification: \${classification}"
    """
}