/*
 * REFERENCE_SELECTOR - Select best reference genomes for a sample
 * 
 * Combines ANI results with quality metadata to select optimal references
 * No external Python dependencies - uses pure bash/awk
 * 
 * Compatible with skani v0.3.x output format (7 columns)
 */
// GCA/GCF deduplication handled upstream in GENUS_DOWNLOAD
process REFERENCE_SELECTOR {
    tag "${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}/classification/references", mode: 'copy'

    container 'ubuntu:22.04'

    input:
    tuple val(sample_id), val(genus), path(ani_results), path(quality_metadata)

    output:
    tuple val(sample_id), path("selected_references.csv"), emit: references
    tuple val(sample_id), path("selection_summary.json"), emit: summary

    script:
    """
    #!/bin/bash
    set -eu
    # Note: removed pipefail to avoid SIGPIPE issues

    echo "Reference selection for ${sample_id} (${genus})"

    # Thresholds
    MIN_ANI=80
    MIN_AF=0.5
    HIGH_ANI=95
    MED_ANI=90
    HIGH_COMP=95
    HIGH_CONT=5
    MED_COMP=90
    MED_CONT=10
    MIN_COMP=90
    MAX_CONT=10
    MAX_REFS=1000

    # Create output header
    echo "accession,species,ani,af_query,completeness,contamination,quality_score,quality_source,selection_reason,ref_file,gff_file" > selected_references.csv

    # Check if we have ANI results
    ani_count=\$(tail -n +2 ${ani_results} 2>/dev/null | wc -l || echo 0)
    
    if [ "\$ani_count" -eq 0 ]; then
        echo "No ANI results available"
        cat > selection_summary.json << EOF
{
    "sample_id": "${sample_id}",
    "genus": "${genus}",
    "status": "no_candidates",
    "total_candidates": 0,
    "selected_count": 0
}
EOF
        exit 0
    fi

    # Process ANI results and join with quality data
    # ANI columns (skani v0.3.x): Ref_file, Query_file, ANI, Align_fraction_ref, Align_fraction_query, Ref_name, Query_name
    # Quality columns: accession,species,family,genus,completeness,contamination,quality_source,assembly_level,file_path

    # Process ANI results directly (no deduplication needed - handled in GENUS_DOWNLOAD)
    echo "Processing ANI results and quality metadata..."
    tail -n +2 ${ani_results} > ani_data.tmp
    
    # Create temp file with merged data
    rm -f all_candidates.tmp
    touch all_candidates.tmp
    
    # Read all 7 columns from skani output (fixed from original 5)
    while IFS=\$'\\t' read -r ref_file _query_file ani _align_fraction_ref af_query _ref_name _query_name; do
        # Extract accession from ref_file path
        accession=\$(basename "\$ref_file" .fna)
        # Also handle .fna.gz files
        accession=\$(echo "\$accession" | sed 's/\\.gz\$//')
        
        # Skip if ANI or AF below threshold
        if awk -v ani="\$ani" -v min_ani="\$MIN_ANI" -v af="\$af_query" -v min_af="\$MIN_AF" \\
            'BEGIN {exit !(ani < min_ani || af < min_af)}'; then
            continue
        fi
        
        # Look up quality data
        quality_line=\$(awk -F',' -v acc="\$accession" '\$1 == acc {print; exit}' ${quality_metadata})
        
        if [ -n "\$quality_line" ]; then
            species=\$(echo "\$quality_line" | cut -d',' -f2)
            completeness=\$(echo "\$quality_line" | cut -d',' -f5)
            contamination=\$(echo "\$quality_line" | cut -d',' -f6)
            quality_source=\$(echo "\$quality_line" | cut -d',' -f7)
        else
            species="unknown"
            completeness=""
            contamination=""
            quality_source="unknown"
        fi
        
        # Calculate quality score if we have data
        if [ -n "\$completeness" ] && [ -n "\$contamination" ]; then
            quality_score=\$(awk -v comp="\$completeness" -v cont="\$contamination" 'BEGIN {printf "%.1f", comp - 5 * cont}')
        else
            quality_score="0"
            completeness="0"
            contamination="100"
        fi
        
        # Determine selection tier using awk with variables
        high_qual=false
        med_qual=false
        high_ani_match=false
        med_ani_match=false
        
        awk -v comp="\$completeness" -v hc="\$HIGH_COMP" -v cont="\$contamination" -v hcont="\$HIGH_CONT" \\
            'BEGIN {exit !(comp >= hc && cont <= hcont)}' && high_qual=true
        awk -v comp="\$completeness" -v mc="\$MED_COMP" -v cont="\$contamination" -v mcont="\$MED_CONT" \\
            'BEGIN {exit !(comp >= mc && cont <= mcont)}' && med_qual=true
        awk -v ani="\$ani" -v ha="\$HIGH_ANI" 'BEGIN {exit !(ani >= ha)}' && high_ani_match=true
        awk -v ani="\$ani" -v ma="\$MED_ANI" 'BEGIN {exit !(ani >= ma)}' && med_ani_match=true
        
        if [ "\$high_qual" = true ] && [ "\$high_ani_match" = true ]; then
            tier=1
            reason="high_quality_high_ani"
        elif [ "\$med_qual" = true ] && [ "\$high_ani_match" = true ]; then
            tier=2
            reason="medium_quality_high_ani"
        elif [ "\$high_qual" = true ] && [ "\$med_ani_match" = true ]; then
            tier=3
            reason="high_quality_medium_ani"
        elif [ "\$med_qual" = true ] && [ "\$med_ani_match" = true ]; then
            tier=4
            reason="medium_quality_medium_ani"
        else
            # Check minimum quality requirements for fallback
            min_qual_ok=false
            awk -v comp="\$completeness" -v minc="\$MIN_COMP" -v cont="\$contamination" -v maxcont="\$MAX_CONT" \\
                'BEGIN {exit !(comp >= minc && cont <= maxcont)}' && min_qual_ok=true

            if [ "\$min_qual_ok" = true ]; then
                tier=5
                reason="fallback"
            else
                # Skip this candidate - doesn't meet minimum quality
                continue
            fi
        fi

        # Determine GFF3 file path from FASTA path
        gff_file=\$(echo "\$ref_file" | sed 's/\\.fna\$/.gff3/' | sed 's/\\.fna\\.gz\$/.gff3/')
        
        # Ensure GFF3 file has correct permissions if it exists
        if [ -f "\$gff_file" ]; then
            chmod 644 "\$gff_file" 2>/dev/null || true
        fi
        
        # Output: tier, ani, quality_score, accession, species, af_query, completeness, contamination, quality_source, reason, ref_file, gff_file
        echo "\$tier,\$ani,\$quality_score,\$accession,\$species,\$af_query,\$completeness,\$contamination,\$quality_source,\$reason,\$ref_file,\$gff_file" >> all_candidates.tmp
        
    done < ani_data.tmp
    
    # Check if we have any candidates
    if [ ! -s all_candidates.tmp ]; then
        echo "No candidates passed filtering thresholds"
        cat > selection_summary.json << EOF
{
    "sample_id": "${sample_id}",
    "genus": "${genus}",
    "status": "no_candidates_after_filter",
    "total_candidates": \$ani_count,
    "selected_count": 0
}
EOF
        rm -f ani_data.tmp all_candidates.tmp
        exit 0
    fi
    
    # Sort all candidates by tier, then by ANI, then by quality
    sort -t',' -k1,1n -k2,2nr -k3,3nr all_candidates.tmp > sorted_candidates.tmp
    
    # Implement tiered selection strategy
    echo "Implementing tiered selection strategy..."
    
    # Count tier 1 references
    tier1_count=\$(awk -F',' '\$1 == 1' sorted_candidates.tmp | wc -l)
    echo "Available tier 1 references: \$tier1_count"
    
    rm -f candidates.tmp
    
    if [ "\$tier1_count" -ge 10 ]; then
        # Use all tier 1 references if we have at least 10
        echo "Using all \$tier1_count tier 1 references (high quality + high ANI)"
        awk -F',' '\$1 == 1' sorted_candidates.tmp > candidates.tmp
    else
        # Use all tier 1 and fill with best from other tiers
        echo "Using all \$tier1_count tier 1 references, filling with best from other tiers"
        
        # Add all tier 1 references
        awk -F',' '\$1 == 1' sorted_candidates.tmp > candidates.tmp
        
        # Calculate how many more we need, up to MAX_REFS total
        needed=\$((MAX_REFS - tier1_count))
        
        if [ "\$needed" -gt 0 ]; then
            # Add best references from tiers 2-5 to reach MAX_REFS
            awk -F',' '\$1 > 1' sorted_candidates.tmp | head -n \$needed >> candidates.tmp
            total_selected=\$(wc -l < candidates.tmp)
            echo "Added \$needed references from lower tiers (total: \$total_selected)"
        fi
    fi
    
    rm -f ani_data.tmp all_candidates.tmp sorted_candidates.tmp

    # Convert to final output format
    selected_count=0
    while IFS=',' read -r tier ani quality_score accession species af_query completeness contamination quality_source reason ref_file gff_file; do
        echo "\$accession,\$species,\$ani,\$af_query,\$completeness,\$contamination,\$quality_score,\$quality_source,\$reason,\$ref_file,\$gff_file" >> selected_references.csv
        selected_count=\$((selected_count + 1))
    done < candidates.tmp

    # Get top reference info for summary
    if [ "\$selected_count" -gt 0 ]; then
        top_line=\$(head -1 candidates.tmp)
        top_accession=\$(echo "\$top_line" | cut -d',' -f4)
        top_species=\$(echo "\$top_line" | cut -d',' -f5)
        top_ani=\$(echo "\$top_line" | cut -d',' -f2)
        top_quality=\$(echo "\$top_line" | cut -d',' -f3)
        
        # Calculate mean ANI
        mean_ani=\$(awk -F',' '{sum+=\$2; count++} END {if(count>0) printf "%.2f", sum/count; else print "0"}' candidates.tmp)
    else
        top_accession="none"
        top_species="none"
        top_ani="0"
        top_quality="0"
        mean_ani="0"
    fi

    # Count by selection reason (ensure clean integers)
    tier1=\$(grep -c "high_quality_high_ani" candidates.tmp 2>/dev/null || true)
    tier1=\${tier1:-0}
    [ -z "\$tier1" ] && tier1=0
    tier2=\$(grep -c "medium_quality_high_ani" candidates.tmp 2>/dev/null || true)
    tier2=\${tier2:-0}
    [ -z "\$tier2" ] && tier2=0
    tier3=\$(grep -c "high_quality_medium_ani" candidates.tmp 2>/dev/null || true)
    tier3=\${tier3:-0}
    [ -z "\$tier3" ] && tier3=0
    tier4=\$(grep -c "medium_quality_medium_ani" candidates.tmp 2>/dev/null || true)
    tier4=\${tier4:-0}
    [ -z "\$tier4" ] && tier4=0
    tier5=\$(grep -c "fallback" candidates.tmp 2>/dev/null || true)
    tier5=\${tier5:-0}
    [ -z "\$tier5" ] && tier5=0

    # Create summary JSON
    cat > selection_summary.json << EOF
{
    "sample_id": "${sample_id}",
    "genus": "${genus}",
    "status": "\$([ \$selected_count -gt 0 ] && echo 'success' || echo 'no_selection')",
    "total_candidates": \$ani_count,
    "selected_count": \$selected_count,
    "mean_ani": \$mean_ani,
    "top_reference": {
        "accession": "\$top_accession",
        "species": "\$top_species",
        "ani": \$top_ani,
        "quality_score": \$top_quality
    },
    "selection_breakdown": {
        "high_quality_high_ani": \$tier1,
        "medium_quality_high_ani": \$tier2,
        "high_quality_medium_ani": \$tier3,
        "medium_quality_medium_ani": \$tier4,
        "fallback": \$tier5
    }
}
EOF

    rm -f candidates.tmp

    echo "Selected \$selected_count references"
    [ "\$selected_count" -gt 0 ] && echo "Top: \$top_accession (\$top_species) - ANI: \$top_ani%"
    """
}
