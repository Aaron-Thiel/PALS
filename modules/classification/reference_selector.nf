/*
 * REFERENCE_SELECTOR - Select best reference genomes for a sample
 *
 * Uses composite scoring: score = (ANI × 0.8) + (AF³ × 0.1) + (Quality³ × 0.1)
 * ANI is primary factor; AF and quality penalize only when low (cubic)
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

    # Hard filter thresholds
    MIN_ANI=80
    MIN_AF=0.5
    MIN_COMP=90
    MAX_CONT=10

    # Selection parameters
    MIN_REFS=10
    MAX_REFS=100
    SCORE_CUTOFF=0.88

    # Create output header (include sample_id for downstream matching)
    echo "sample_id,accession,species,ani,af_query,completeness,contamination,quality_score,composite_score,quality_source,selection_reason,ref_file,gff_file" > selected_references.csv

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

    echo "Processing ANI results and quality metadata..."
    tail -n +2 ${ani_results} > ani_data.tmp

    # Create temp file with merged data
    rm -f all_candidates.tmp
    touch all_candidates.tmp

    # Read all 7 columns from skani output
    while IFS=\$'\\t' read -r ref_file _query_file ani _align_fraction_ref af_query _ref_name _query_name; do
        # Extract accession from ref_file path
        accession=\$(basename "\$ref_file" .fna)
        accession=\$(echo "\$accession" | sed 's/\\.gz\$//')

        # Skip if ANI or AF below hard thresholds
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

        # Skip if missing quality data or below minimum quality thresholds
        if [ -z "\$completeness" ] || [ -z "\$contamination" ]; then
            continue
        fi

        # Check minimum quality thresholds
        if awk -v comp="\$completeness" -v minc="\$MIN_COMP" -v cont="\$contamination" -v maxcont="\$MAX_CONT" \\
            'BEGIN {exit !(comp < minc || cont > maxcont)}'; then
            continue
        fi

        # Calculate quality score (completeness - 5 * contamination)
        quality_score=\$(awk -v comp="\$completeness" -v cont="\$contamination" 'BEGIN {printf "%.1f", comp - 5 * cont}')

        # Calculate composite score: (ANI/100 * 0.8) + ((AF/100)^3 * 0.1) + ((Quality/100)^3 * 0.1)
        composite_score=\$(awk -v ani="\$ani" -v af="\$af_query" -v qual="\$quality_score" 'BEGIN {
            ani_norm = ani / 100
            af_norm = af / 100
            qual_norm = qual / 100
            if (qual_norm < 0) qual_norm = 0
            score = (ani_norm * 0.8) + ((af_norm ^ 3) * 0.1) + ((qual_norm ^ 3) * 0.1)
            printf "%.4f", score
        }')

        # Determine GFF3 file path from FASTA path
        gff_file=\$(echo "\$ref_file" | sed 's/\\.fna\$/.gff3/' | sed 's/\\.fna\\.gz\$/.gff3/')

        # Ensure GFF3 file has correct permissions if it exists
        if [ -f "\$gff_file" ]; then
            chmod 644 "\$gff_file" 2>/dev/null || true
        fi

        # Output: composite_score, ani, quality_score, af_query, accession, species, completeness, contamination, quality_source, ref_file, gff_file
        echo "\$composite_score,\$ani,\$quality_score,\$af_query,\$accession,\$species,\$completeness,\$contamination,\$quality_source,\$ref_file,\$gff_file" >> all_candidates.tmp

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

    # Sort all candidates by composite score (descending)
    sort -t',' -k1,1nr all_candidates.tmp > sorted_candidates.tmp

    echo "Implementing score-based selection..."

    # Count candidates above score cutoff
    above_cutoff=\$(awk -F',' -v cutoff="\$SCORE_CUTOFF" '\$1 >= cutoff' sorted_candidates.tmp | wc -l)
    total_candidates=\$(wc -l < sorted_candidates.tmp)
    echo "Candidates above score cutoff (\$SCORE_CUTOFF): \$above_cutoff / \$total_candidates"

    rm -f candidates.tmp

    # Selection logic:
    # 1. Take all refs with score >= SCORE_CUTOFF (up to MAX_REFS)
    # 2. If count < MIN_REFS, fill from lower scores until MIN_REFS reached
    if [ "\$above_cutoff" -ge "\$MIN_REFS" ]; then
        # Take all above cutoff, up to MAX_REFS
        awk -F',' -v cutoff="\$SCORE_CUTOFF" '\$1 >= cutoff' sorted_candidates.tmp | head -n \$MAX_REFS > candidates.tmp
        echo "Selected \$(wc -l < candidates.tmp) references above score cutoff"
    else
        # Take all above cutoff + fill to MIN_REFS from below cutoff
        awk -F',' -v cutoff="\$SCORE_CUTOFF" '\$1 >= cutoff' sorted_candidates.tmp > candidates.tmp
        current=\$(wc -l < candidates.tmp)
        needed=\$((MIN_REFS - current))
        if [ "\$needed" -gt 0 ]; then
            awk -F',' -v cutoff="\$SCORE_CUTOFF" '\$1 < cutoff' sorted_candidates.tmp | head -n \$needed >> candidates.tmp
            echo "Selected \$current refs above cutoff + \$needed fallbacks to reach MIN_REFS"
        fi
    fi

    rm -f ani_data.tmp all_candidates.tmp sorted_candidates.tmp

    # Convert to final output format
    # Input format: composite_score, ani, quality_score, af_query, accession, species, completeness, contamination, quality_source, ref_file, gff_file
    selected_count=0
    good_ref_count=0
    fallback_count=0

    while IFS=',' read -r composite_score ani quality_score af_query accession species completeness contamination quality_source ref_file gff_file; do
        # Determine selection reason based on score
        if awk -v score="\$composite_score" -v cutoff="\$SCORE_CUTOFF" 'BEGIN {exit !(score >= cutoff)}'; then
            reason="good_reference"
            good_ref_count=\$((good_ref_count + 1))
        else
            reason="fallback"
            fallback_count=\$((fallback_count + 1))
        fi

        echo "${sample_id},\$accession,\$species,\$ani,\$af_query,\$completeness,\$contamination,\$quality_score,\$composite_score,\$quality_source,\$reason,\$ref_file,\$gff_file" >> selected_references.csv
        selected_count=\$((selected_count + 1))
    done < candidates.tmp

    # Get top reference info for summary
    if [ "\$selected_count" -gt 0 ]; then
        top_line=\$(head -1 candidates.tmp)
        top_accession=\$(echo "\$top_line" | cut -d',' -f5)
        top_species=\$(echo "\$top_line" | cut -d',' -f6)
        top_ani=\$(echo "\$top_line" | cut -d',' -f2)
        top_quality=\$(echo "\$top_line" | cut -d',' -f3)
        top_composite=\$(echo "\$top_line" | cut -d',' -f1)

        # Calculate mean ANI and mean composite score
        mean_ani=\$(awk -F',' '{sum+=\$2; count++} END {if(count>0) printf "%.2f", sum/count; else print "0"}' candidates.tmp)
        mean_score=\$(awk -F',' '{sum+=\$1; count++} END {if(count>0) printf "%.4f", sum/count; else print "0"}' candidates.tmp)
        min_score=\$(tail -1 candidates.tmp | cut -d',' -f1)
    else
        top_accession="none"
        top_species="none"
        top_ani="0"
        top_quality="0"
        top_composite="0"
        mean_ani="0"
        mean_score="0"
        min_score="0"
    fi

    # Create summary JSON
    cat > selection_summary.json << EOF
{
    "sample_id": "${sample_id}",
    "genus": "${genus}",
    "status": "\$([ \$selected_count -gt 0 ] && echo 'success' || echo 'no_selection')",
    "total_candidates": \$ani_count,
    "selected_count": \$selected_count,
    "score_cutoff": \$SCORE_CUTOFF,
    "mean_ani": \$mean_ani,
    "mean_composite_score": \$mean_score,
    "min_composite_score": \$min_score,
    "top_reference": {
        "accession": "\$top_accession",
        "species": "\$top_species",
        "ani": \$top_ani,
        "quality_score": \$top_quality,
        "composite_score": \$top_composite
    },
    "selection_breakdown": {
        "good_reference": \$good_ref_count,
        "fallback": \$fallback_count
    }
}
EOF

    rm -f candidates.tmp

    echo "Selected \$selected_count references (\$good_ref_count good, \$fallback_count fallback)"
    [ "\$selected_count" -gt 0 ] && echo "Top: \$top_accession (\$top_species) - ANI: \$top_ani%, Score: \$top_composite"
    """
}
