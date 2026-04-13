/*
 * BAKTA_GENUS - Annotate reference genomes missing GFF3 files
 *
 * Only runs Bakta on genomes where gff_path is empty in metadata
 * Updates metadata with annotation source (NCBI or Bakta)
 * Uses storeDir for caching at genus level
 *
 * INTERMEDIATE CACHING: GFF3 files are written directly to the database
 * directory as they complete. On timeout/restart, already-annotated
 * genomes are skipped, preserving progress across multiple runs.
 */

process BAKTA_GENUS {
    tag "${genus}"

    // Use storeDir for caching - stores outputs in the genus database directory
    storeDir "${params.genome_database}/${family}/${genus}"

    container 'staphb/bakta:latest'

    input:
    tuple val(genus), val(family), path("genome_metadata_checkm2.csv")

    output:
    tuple val(genus), val(family), path("genome_metadata.csv"), emit: metadata
    tuple val(genus), val(family), path("bakta_summary.json"), emit: summary

    script:
    """
    #!/bin/bash
    set -eu

    echo "════════════════════════════════════════════════════════════════"
    echo "BAKTA ANNOTATION: ${genus}"
    echo "════════════════════════════════════════════════════════════════"

    GENUS_DIR="${params.genome_database}/${family}/${genus}"
    BAKTA_DB=\${BAKTA_DB:-/databases/bakta/db}

    # Step 1: Analyze input metadata - count genomes missing GFF3
    echo -e "\\n[1/3] Analyzing metadata..."

    total_count=\$(tail -n +2 genome_metadata_checkm2.csv | wc -l)

    # Count genomes with empty gff_path (column 10)
    missing_gff=\$(awk -F',' 'NR>1 && (\$10 == "" || \$10 ~ /^[[:space:]]*\$/) {count++} END {print count+0}' genome_metadata_checkm2.csv)
    has_gff=\$(awk -F',' 'NR>1 && \$10 != "" && \$10 !~ /^[[:space:]]*\$/ {count++} END {print count+0}' genome_metadata_checkm2.csv)
    
    echo "Total genomes: \$total_count"
    echo "With GFF3 (NCBI): \$has_gff"
    echo "Missing GFF3: \$missing_gff"

    # Handle case where no genomes need annotation
    if [ "\$total_count" -eq 0 ]; then
        echo "No genomes in metadata"
        # Add gff_source column to header if not already present
        header=\$(head -1 genome_metadata_checkm2.csv)
        if echo "\$header" | grep -q "gff_source"; then
            echo "\$header" > genome_metadata.csv
        else
            echo "\$header,gff_source" > genome_metadata.csv
        fi
        # Create JSON using printf to avoid heredoc variable expansion issues
        printf '{
    "genus": "%s",
    "family": "%s",
    "status": "empty",
    "reason": "no_genomes_in_metadata",
    "total_genomes": 0,
    "bakta_annotated": 0
}\n' "${genus}" "${family}" > bakta_summary.json
        exit 0
    fi

    if [ "\$missing_gff" -eq 0 ]; then
        echo "✓ All genomes have GFF3 annotations - skipping Bakta"
        # Add gff_source column - only set to NCBI if gff_path is not empty
        # Use awk for robust CSV handling and to avoid carriage return issues
        awk -F',' 'BEGIN {OFS=","}
            NR==1 {
                # Header: add gff_source if not present
                if (\$0 ~ /gff_source/) {
                    gsub(/\\r/, ""); print
                } else {
                    gsub(/\\r/, ""); print \$0 ",gff_source"
                }
                next
            }
            {
                # Remove any carriage returns
                gsub(/\\r/, "")
                # Check if gff_path (column 10) is empty
                if (\$10 == "" || \$10 ~ /^[[:space:]]*\$/) {
                    print \$0 ",none"
                } else {
                    print \$0 ",NCBI"
                }
            }
        ' genome_metadata_checkm2.csv > genome_metadata.csv

        # Count NCBI vs none
        final_ncbi=\$(grep -c ",NCBI\$" genome_metadata.csv 2>/dev/null || echo "0")
        final_none=\$(grep -c ",none\$" genome_metadata.csv 2>/dev/null || echo "0")

        # Create JSON using printf to avoid heredoc variable expansion issues
        printf '{
    "genus": "%s",
    "family": "%s",
    "status": "skipped",
    "reason": "all_have_gff",
    "total_genomes": %d,
    "ncbi_annotated": %d,
    "bakta_annotated": 0,
    "no_annotation": %d
}\\n' "${genus}" "${family}" \$total_count \$final_ncbi \$final_none > bakta_summary.json
        exit 0
    fi

    # Step 2: Run Bakta on genomes missing GFF3
    # Quality thresholds - only annotate genomes meeting minimum quality
    MIN_COMP=90
    MAX_CONT=10

    # Get list of accessions missing GFF3 that meet quality thresholds
    # Columns: 1=accession, 5=completeness, 6=contamination, 9=file_path, 10=gff_path
    awk -F',' -v min_comp="\$MIN_COMP" -v max_cont="\$MAX_CONT" '
        NR>1 && (\$10 == "" || \$10 ~ /^[[:space:]]*\$/) {
            # Check quality thresholds
            comp = \$5 + 0
            cont = \$6 + 0
            if (comp >= min_comp && cont <= max_cont) {
                print \$1, \$9
            }
        }
    ' genome_metadata_checkm2.csv > missing_list.txt

    qualified_count=\$(wc -l < missing_list.txt | tr -d ' ')
    skipped_low_quality=\$((missing_gff - qualified_count))

    echo -e "\\n[2/3] Running Bakta on \$qualified_count genomes (skipping \$skipped_low_quality low-quality genomes)..."
    echo "  Quality thresholds: completeness >= \$MIN_COMP%, contamination <= \$MAX_CONT%"

    mkdir -p bakta_work
    bakta_success=0

    # Run Bakta on each genome
    # Track genomes skipped because GFF3 already exists in database (from previous partial run)
    skipped_existing=0

    if [ -d "\$BAKTA_DB" ]; then
        while read -r accession file_path; do
            # Get species directory from file path
            species_dir=\$(echo "\$file_path" | sed 's|^[.]/||' | xargs dirname)

            # CHECKPOINT: Skip if GFF3 already exists in database (from previous partial run)
            if [ -f "\$GENUS_DIR/\$species_dir/\${accession}.gff3" ]; then
                echo "  ✓ Skipping \$accession - GFF3 already exists in database"
                skipped_existing=\$((skipped_existing + 1))
                continue
            fi

            # Find the genome file in the genus directory
            fasta=\$(find "\$GENUS_DIR" -name "\${accession}.fna" -type f 2>/dev/null | head -1)

            if [ ! -f "\$fasta" ]; then
                echo "  WARNING: Could not find \$accession.fna"
                continue
            fi

            echo "  Annotating: \$accession (\$species_dir)"

            # Create temp output directory for this accession
            outdir="bakta_work/\${accession}"
            mkdir -p "\$outdir"

            if bakta \\
                --db "\$BAKTA_DB" \\
                --output "\$outdir" \\
                --prefix "\$accession" \\
                --threads ${task.cpus} \\
                --skip-plot \\
                --force \\
                "\$fasta" 2>&1; then

                # Write GFF3 directly to database directory (survives timeout)
                if [ -f "\$outdir/\${accession}.gff3" ]; then
                    mkdir -p "\$GENUS_DIR/\$species_dir"
                    cp "\$outdir/\${accession}.gff3" "\$GENUS_DIR/\$species_dir/"
                    bakta_success=\$((bakta_success + 1))
                    echo "    ✓ Saved to database: \$GENUS_DIR/\$species_dir/\${accession}.gff3"
                fi
            else
                echo "    ✗ Failed to annotate \$accession"
            fi

            # Clean up bakta output immediately (only keep GFF3)
            rm -rf "\$outdir"

        done < missing_list.txt
    else
        echo "WARNING: Bakta database not found at \$BAKTA_DB"
        echo "Skipping annotation - genomes will remain without GFF3"
    fi

    # Cleanup work directory
    rm -rf bakta_work missing_list.txt

    echo "Successfully annotated: \$bakta_success genomes"
    if [ "\$skipped_existing" -gt 0 ]; then
        echo "Skipped (already in database): \$skipped_existing genomes"
    fi

    # Step 3: Update metadata with gff_path and gff_source
    echo -e "\\n[3/3] Updating metadata..."

    # Create a lookup file of all GFF3 files in the database directory
    # This includes both NCBI GFF3s (from merge) and Bakta-generated ones
    find "\$GENUS_DIR" -name "*.gff3" -type f 2>/dev/null | while read -r gff; do
        acc=\$(basename "\$gff" .gff3)
        # Store relative path from GENUS_DIR
        rel_path=\$(echo "\$gff" | sed "s|^\$GENUS_DIR/||")
        echo "\$acc \$rel_path"
    done > bakta_gff_lookup.txt

    # Use awk for robust CSV processing - handles carriage returns properly
    awk -F',' -v lookup="bakta_gff_lookup.txt" '
    BEGIN {
        OFS=","
        # Load Bakta GFF lookup
        while ((getline line < lookup) > 0) {
            split(line, parts, " ")
            bakta_gff[parts[1]] = parts[2]
        }
        close(lookup)
    }
    NR==1 {
        # Header: add gff_source if not present
        gsub(/\\r/, "")
        if (\$0 ~ /gff_source/) {
            print
        } else {
            print \$0 ",gff_source"
        }
        next
    }
    {
        # Remove carriage returns
        gsub(/\\r/, "")

        accession = \$1
        gff_path = \$10

        # Determine gff_source and potentially update gff_path
        if (gff_path != "" && gff_path !~ /^[[:space:]]*\$/) {
            # Already has GFF from NCBI
            print \$0 ",NCBI"
        } else if (accession in bakta_gff) {
            # Bakta created a GFF3 - update gff_path
            \$10 = bakta_gff[accession]
            print \$0 ",Bakta"
        } else {
            # No annotation available
            print \$0 ",none"
        }
    }
    ' genome_metadata_checkm2.csv > genome_metadata.csv

    rm -f bakta_gff_lookup.txt

    # Count final statistics (ensure clean integers)
    final_ncbi=\$(grep -c ",NCBI\$" genome_metadata.csv 2>/dev/null || true)
    final_ncbi=\${final_ncbi:-0}
    [ -z "\$final_ncbi" ] && final_ncbi=0
    final_bakta=\$(grep -c ",Bakta\$" genome_metadata.csv 2>/dev/null || true)
    final_bakta=\${final_bakta:-0}
    [ -z "\$final_bakta" ] && final_bakta=0
    final_none=\$(grep -c ",none\$" genome_metadata.csv 2>/dev/null || true)
    final_none=\${final_none:-0}
    [ -z "\$final_none" ] && final_none=0

    # Create summary JSON using printf to avoid heredoc variable expansion issues
    TIMESTAMP=\$(date -Iseconds)
    printf '{
    "genus": "%s",
    "family": "%s",
    "status": "success",
    "total_genomes": %d,
    "ncbi_annotated": %d,
    "bakta_annotated": %d,
    "no_annotation": %d,
    "timestamp": "%s"
}\n' "${genus}" "${family}" \$total_count \$final_ncbi \$final_bakta \$final_none "\$TIMESTAMP" > bakta_summary.json

    echo -e "\\n════════════════════════════════════════════════════════════════"
    echo "COMPLETED: Bakta annotation for ${genus}"
    echo "  - NCBI GFF3: \$final_ncbi"
    echo "  - Bakta GFF3: \$final_bakta"
    echo "  - No annotation: \$final_none"
    echo "════════════════════════════════════════════════════════════════"
    """
}
