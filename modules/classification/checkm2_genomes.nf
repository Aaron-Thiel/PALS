/*
 * CHECKM2_GENUS - Quality assessment for genomes missing NCBI quality data
 *
 * Only runs CheckM2 on genomes with quality_source="pending"
 * If all genomes have quality data, completes quickly
 *
 * INTERMEDIATE CACHING: Results are written to a checkpoint file
 * (checkm2_progress.tsv) in the database directory after each batch.
 * On timeout/restart, already-assessed genomes are skipped, preserving
 * progress across multiple runs. Checkpoint is cleaned up on completion.
 */

process CHECKM2_GENUS {
    tag "${genus}"
    
    storeDir "${params.genome_database}/${family}/${genus}"

    container 'staphb/checkm2:latest'

    input:
    tuple val(genus), val(family), path("genome_metadata_merged.csv")

    output:
    tuple val(genus), val(family), path("genome_metadata_checkm2.csv"), emit: metadata
    tuple val(genus), val(family), path("checkm2_summary.json"), emit: summary

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "════════════════════════════════════════════════════════════════"
    echo "CHECKM2 QUALITY ASSESSMENT: ${genus}"
    echo "════════════════════════════════════════════════════════════════"

    GENUS_DIR="${params.genome_database}/${family}/${genus}"

    # Step 1: Analyze input metadata
    echo -e "\\n[1/3] Analyzing metadata..."

    total_count=\$(tail -n +2 genome_metadata_merged.csv | wc -l)
    pending_count=\$(awk -F',' '\$7 == "pending" {count++} END {print count+0}' genome_metadata_merged.csv)
    
    echo "Total genomes: \$total_count"
    echo "Pending quality assessment: \$pending_count"

    # Handle case where no genomes exist or all have quality data
    if [ "\$total_count" -eq 0 ]; then
        echo "No genomes in metadata - copying input as-is"
        cp genome_metadata_merged.csv genome_metadata_checkm2.csv
        cat > checkm2_summary.json << EOF
{
    "genus": "${genus}",
    "family": "${family}",
    "status": "empty",
    "reason": "no_genomes_in_metadata",
    "total_genomes": 0,
    "checkm2_analyzed": 0
}
EOF
        exit 0
    fi

    if [ "\$pending_count" -eq 0 ]; then
        echo "✓ All genomes have quality data - skipping CheckM2"
        cp genome_metadata_merged.csv genome_metadata_checkm2.csv
        cat > checkm2_summary.json << EOF
{
    "genus": "${genus}",
    "family": "${family}",
    "status": "skipped",
    "reason": "no_pending_genomes",
    "total_genomes": \$total_count,
    "checkm2_analyzed": 0
}
EOF
        exit 0
    fi

    # Step 2: Run CheckM2 on pending genomes
    echo -e "\\n[2/3] Running CheckM2 on \$pending_count genomes..."

    # CHECKPOINT: Load existing progress from previous partial runs
    PROGRESS_FILE="\$GENUS_DIR/checkm2_progress.tsv"
    declare -A completed_accessions
    skipped_from_checkpoint=0

    if [ -f "\$PROGRESS_FILE" ]; then
        echo "Found checkpoint file - loading previous progress..."
        while IFS=\$'\\t' read -r name comp cont rest; do
            # Skip header
            [[ "\$name" == "Name" ]] && continue
            # Extract accession (remove .fna suffix if present)
            acc=\$(echo "\$name" | sed 's/.fna\$//')
            completed_accessions["\$acc"]=1
        done < "\$PROGRESS_FILE"
        echo "  Loaded \${#completed_accessions[@]} previously completed genomes"
    fi

    mkdir -p checkm2_input checkm2_output

    # Copy pending genomes to input directory (skip already-completed ones)
    while IFS=',' read -r accession species family_col genus_col completeness contamination quality_source assembly_level file_path gff_path rest; do
        # Skip header
        [ "\$accession" = "accession" ] && continue

        # Only process pending genomes
        [ "\$quality_source" != "pending" ] && continue

        # CHECKPOINT: Skip if already completed in previous run
        if [[ -v "completed_accessions[\$accession]" ]]; then
            echo "  ✓ Skipping \$accession - already in checkpoint"
            skipped_from_checkpoint=\$((skipped_from_checkpoint + 1))
            continue
        fi

        # Find the genome file
        fasta=\$(find "\$GENUS_DIR" -name "\${accession}.fna" -type f 2>/dev/null | head -1)
        if [ -f "\$fasta" ]; then
            cp "\$fasta" checkm2_input/
            echo "  Queued: \$accession"
        else
            echo "  WARNING: Could not find \$accession.fna in \$GENUS_DIR"
        fi
    done < genome_metadata_merged.csv

    if [ "\$skipped_from_checkpoint" -gt 0 ]; then
        echo "Skipped \$skipped_from_checkpoint genomes (already in checkpoint)"
    fi

    copied_count=\$(find checkm2_input -name "*.fna" 2>/dev/null | wc -l)
    echo "Copied \$copied_count genomes for analysis"

    checkm2_success=false
    if [ "\$copied_count" -gt 0 ]; then
        # Set database path
        export CHECKM2DB=\${CHECKM2DB:-/databases/checkm2}

        # Adaptive batch processing: try all genomes first, then reduce batch size exponentially
        attempt=-1
        max_attempts=6  # Minimum batch size will be 100/2^5 = 3

        while [ "\$checkm2_success" = "false" ] && [ "\$attempt" -lt "\$max_attempts" ]; do
            attempt=\$((attempt + 1))

            if [ "\$attempt" -eq 0 ]; then
                batch_size=\$copied_count  # Try all genomes first
                echo "Attempt \$((attempt + 1)): Processing all \$copied_count genomes..."
            else
                batch_size=\$((100 / (2 ** (attempt - 1))))  # 100, 50, 25, 12, 6, 3
                [ "\$batch_size" -lt 1 ] && batch_size=1
                echo "Attempt \$((attempt + 1)): Reducing batch size to \$batch_size genomes..."
            fi

            # Clear previous output
            rm -rf checkm2_output/*
            mkdir -p checkm2_output

            # Get list of all genome files
            genome_files=(\$(find checkm2_input -name "*.fna" -type f))
            total_files=\${#genome_files[@]}

            if [ "\$batch_size" -ge "\$total_files" ]; then
                # Process all at once
                if checkm2 predict \\
                    --threads ${task.cpus} \\
                    --input checkm2_input \\
                    --output-directory checkm2_output \\
                    --extension .fna \\
                    --database_path \${CHECKM2DB}/uniref100.KO.1.dmnd \\
                    --force 2>&1; then
                    checkm2_success=true
                    echo "✓ CheckM2 completed successfully"

                    # CHECKPOINT: Write results to persistent progress file
                    if [ -f "checkm2_output/quality_report.tsv" ]; then
                        if [ ! -f "\$PROGRESS_FILE" ]; then
                            cp checkm2_output/quality_report.tsv "\$PROGRESS_FILE"
                        else
                            tail -n +2 checkm2_output/quality_report.tsv >> "\$PROGRESS_FILE"
                        fi
                        echo "Checkpoint: Saved results to database"
                    fi
                else
                    echo "WARNING: CheckM2 failed with batch size \$batch_size"
                fi
            else
                # Process in batches
                batch_num=0
                batch_success=true
                mkdir -p checkm2_batch_input

                for ((i=0; i<total_files; i+=batch_size)); do
                    batch_num=\$((batch_num + 1))
                    rm -rf checkm2_batch_input/*

                    # Copy batch of genomes
                    for ((j=i; j<i+batch_size && j<total_files; j++)); do
                        cp "\${genome_files[j]}" checkm2_batch_input/
                    done

                    batch_count=\$(find checkm2_batch_input -name "*.fna" | wc -l)
                    echo "  Processing batch \$batch_num (\$batch_count genomes)..."

                    mkdir -p checkm2_batch_output
                    if checkm2 predict \\
                        --threads ${task.cpus} \\
                        --input checkm2_batch_input \\
                        --output-directory checkm2_batch_output \\
                        --extension .fna \\
                        --database_path \${CHECKM2DB}/uniref100.KO.1.dmnd \\
                        --force 2>&1; then
                        # Merge results into local output
                        if [ -f "checkm2_batch_output/quality_report.tsv" ]; then
                            if [ ! -f "checkm2_output/quality_report.tsv" ]; then
                                cp checkm2_batch_output/quality_report.tsv checkm2_output/
                            else
                                tail -n +2 checkm2_batch_output/quality_report.tsv >> checkm2_output/quality_report.tsv
                            fi

                            # CHECKPOINT: Also write to persistent progress file in database
                            if [ ! -f "\$PROGRESS_FILE" ]; then
                                cp checkm2_batch_output/quality_report.tsv "\$PROGRESS_FILE"
                                echo "    Checkpoint: Created progress file"
                            else
                                tail -n +2 checkm2_batch_output/quality_report.tsv >> "\$PROGRESS_FILE"
                                echo "    Checkpoint: Appended batch results"
                            fi
                        fi
                        echo "  ✓ Batch \$batch_num completed"
                    else
                        echo "  ✗ Batch \$batch_num failed"
                        batch_success=false
                        break
                    fi
                    rm -rf checkm2_batch_output
                done

                rm -rf checkm2_batch_input

                if [ "\$batch_success" = "true" ]; then
                    checkm2_success=true
                    echo "✓ CheckM2 completed successfully with batch size \$batch_size"
                fi
            fi
        done

        if [ "\$checkm2_success" = "false" ]; then
            echo "WARNING: CheckM2 failed after \$max_attempts attempts, will keep pending status"
        fi
    fi

    # Step 3: Update metadata with CheckM2 results
    echo -e "\\n[3/3] Updating metadata..."

    # Use checkpoint file as the primary source of results (includes previous + current runs)
    # Fall back to local output if checkpoint doesn't exist
    RESULTS_FILE="\$PROGRESS_FILE"
    if [ ! -f "\$RESULTS_FILE" ]; then
        RESULTS_FILE="checkm2_output/quality_report.tsv"
    fi

    # Create updated metadata using bash (no Python dependency)
    {
        # Header
        head -1 genome_metadata_merged.csv

        # Process each line - read all fields including gff_path (10th column)
        tail -n +2 genome_metadata_merged.csv | while IFS=',' read -r accession species family_col genus_col completeness contamination quality_source assembly_level file_path gff_path rest; do
            if [ "\$quality_source" = "pending" ]; then
                # Try to get CheckM2 results from checkpoint or local output
                if [ -f "\$RESULTS_FILE" ]; then
                    checkm2_line=\$(awk -F'\\t' -v acc="\$accession" '\$1 == acc".fna" || \$1 == acc {print; exit}' "\$RESULTS_FILE")
                    if [ -n "\$checkm2_line" ]; then
                        new_completeness=\$(echo "\$checkm2_line" | cut -f2)
                        new_contamination=\$(echo "\$checkm2_line" | cut -f3)
                        echo "\$accession,\$species,\$family_col,\$genus_col,\$new_completeness,\$new_contamination,CheckM2,\$assembly_level,\$file_path,\$gff_path"
                        continue
                    fi
                fi
            fi
            # Keep original line if no update
            echo "\$accession,\$species,\$family_col,\$genus_col,\$completeness,\$contamination,\$quality_source,\$assembly_level,\$file_path,\$gff_path"
        done
    } > genome_metadata_checkm2.csv

    # Count final statistics (ensure clean integers)
    final_total=\$(tail -n +2 genome_metadata_checkm2.csv | wc -l | tr -d '[:space:]')
    final_ncbi=\$(grep -c ",NCBI," genome_metadata_checkm2.csv 2>/dev/null || true)
    final_ncbi=\${final_ncbi:-0}
    [ -z "\$final_ncbi" ] && final_ncbi=0
    final_checkm2=\$(grep -c ",CheckM2," genome_metadata_checkm2.csv 2>/dev/null || true)
    final_checkm2=\${final_checkm2:-0}
    [ -z "\$final_checkm2" ] && final_checkm2=0
    final_pending=\$(grep -c ",pending," genome_metadata_checkm2.csv 2>/dev/null || true)
    final_pending=\${final_pending:-0}
    [ -z "\$final_pending" ] && final_pending=0

    # Create summary JSON
    cat > checkm2_summary.json << EOF
{
    "genus": "${genus}",
    "family": "${family}",
    "status": "success",
    "total_genomes": \$final_total,
    "ncbi_quality": \$final_ncbi,
    "checkm2_analyzed": \$final_checkm2,
    "still_pending": \$final_pending,
    "timestamp": "\$(date -Iseconds)"
}
EOF

    # Cleanup
    rm -rf checkm2_input

    # Clean up checkpoint file only if all pending genomes have been processed
    # (no remaining pending means the checkpoint served its purpose)
    if [ "\$final_pending" -eq 0 ] && [ -f "\$PROGRESS_FILE" ]; then
        rm -f "\$PROGRESS_FILE"
        echo "Checkpoint file cleaned up (all genomes processed)"
    elif [ -f "\$PROGRESS_FILE" ]; then
        echo "Checkpoint file retained at \$PROGRESS_FILE (\$final_pending genomes still pending)"
    fi

    echo -e "\\n════════════════════════════════════════════════════════════════"
    echo "COMPLETED: CheckM2 analysis for ${genus}"
    echo "  - NCBI quality: \$final_ncbi"
    echo "  - CheckM2 analyzed: \$final_checkm2"
    echo "  - Still pending: \$final_pending"
    if [ "\$skipped_from_checkpoint" -gt 0 ]; then
        echo "  - From checkpoint: \$skipped_from_checkpoint"
    fi
    echo "════════════════════════════════════════════════════════════════"
    """
}
