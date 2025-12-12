/*
 * GENUS_DOWNLOAD - Download reference genomes from NCBI
 * 
 * Uses storeDir for caching - if genus already downloaded, process is SKIPPED
 * Downloads both FASTA and GFF3 annotation files
 */
process GENUS_DOWNLOAD {
    tag "${genus}"

    storeDir "${params.genome_database}/${family}/${genus}"

    container 'staphb/ncbi-datasets:latest'

    input:
    tuple val(sample_id), val(genus), val(family)

    output:
    tuple val(genus), val(family), path("genome_metadata_download.csv"), emit: metadata
    tuple val(genus), val(family), path("download_info.json"), emit: info
    tuple val(genus), val(family), path("*/*.fna"), emit: genomes, optional: true
    tuple val(genus), val(family), path("*/*.gff"), emit: annotations, optional: true
    tuple val(genus), val(family), path("problematic_accessions.txt"), emit: problematic, optional: true
    tuple val(genus), val(family), path("annotations.txt"), emit: annotations_list, optional: true

    when:
    genus && genus != "Unclassified" && genus != "Unknown"

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "════════════════════════════════════════════════════════════════"
    echo "GENUS DOWNLOAD: ${genus}"
    echo "Family: ${family}"
    echo "════════════════════════════════════════════════════════════════"

    # Initialize counters
    ASSEMBLY_COUNT=0
    DOWNLOAD_COUNT=0
    DOWNLOADED_COUNT=0
    SPECIES_COUNT=0
    GFF_COUNT=0

    # Step 1: Query NCBI for assemblies
    echo -e "\\n[1/4] Querying NCBI for ${genus} genomes..."
    
    if ! datasets summary genome taxon "${genus}" --exclude-atypical > search_results.json 2>query_error.log; then
        echo "Warning: NCBI query failed"
        cat query_error.log || true
        echo '{"reports": []}' > search_results.json
    fi
    
    echo "Query result size: \$(wc -c < search_results.json) bytes"
    
    # Extract accessions
    if command -v jq >/dev/null 2>&1; then
        echo "Using jq for JSON parsing"
        ASSEMBLY_COUNT=\$(jq '.reports | length' search_results.json 2>/dev/null || echo "0")
        jq -r '.reports[]? | .accession' search_results.json > accessions.txt 2>/dev/null || touch accessions.txt
    else
        echo "jq not available, using grep fallback"
        grep -oE '"accession"[[:space:]]*:[[:space:]]*"GC[FA]_[0-9]+[.][0-9]+"' search_results.json 2>/dev/null | \\
            grep -oE 'GC[FA]_[0-9]+[.][0-9]+' > accessions.txt || touch accessions.txt
        ASSEMBLY_COUNT=\$(wc -l < accessions.txt)
    fi
    
    echo "First 3 accessions (before deduplication):"
    head -3 accessions.txt || echo "(none)"

    # Step 1.5: Deduplicate GCA/GCF pairs, prioritizing GCF
    echo -e "\\n[1.5/4] Deduplicating GCA/GCF pairs (prioritizing GCF)..."
    
    if [ -s accessions.txt ]; then
        # Create temporary file for deduplicated accessions
        > accessions_dedup.tmp
        
        # Process each accession to find GCA/GCF pairs
        while read -r accession; do
            # Extract numerical part (e.g., GCA_001434575.1 -> 001434575.1)
            numerical=\$(echo "\$accession" | sed 's/^GC[AF]_//')
            echo "\$numerical \$accession" >> accession_pairs.tmp
        done < accessions.txt
        
        # Sort by numerical part, then by accession (GCF comes before GCA lexicographically)
        sort accession_pairs.tmp | awk '{
            numerical = \$1
            accession = \$2
            
            if (numerical != prev_numerical) {
                # New numerical ID - always include the first one
                print accession
                prev_numerical = numerical
            } else {
                # Same numerical ID - prefer GCF over GCA
                if (index(accession, "GCF_") == 1) {
                    # This is GCF, replace the previous GCA entry
                    # Remove the last line and add GCF
                    "tail -n 1 accessions_dedup.tmp" | getline last_line
                    if (index(last_line, "GCA_") == 1) {
                        # Previous was GCA, replace with GCF
                        system("head -n -1 accessions_dedup.tmp > temp_dedup.txt && mv temp_dedup.txt accessions_dedup.tmp")
                        print accession
                    } else {
                        # Previous was also GCF, keep the existing one
                    }
                }
                # If current is GCA and we already have this numerical ID, skip it
            }
        }' > accessions_dedup.tmp
        
        # Use a simpler approach that definitely works
        sort accession_pairs.tmp | awk '{
            numerical = \$1
            accession = \$2
            
            if (!(numerical in seen)) {
                # First time seeing this numerical ID
                best[numerical] = accession
                seen[numerical] = 1
            } else {
                # Already seen this numerical ID - prefer GCF
                if (index(accession, "GCF_") == 1) {
                    best[numerical] = accession
                }
            }
        } END {
            for (num in best) {
                print best[num]
            }
        }' | sort > accessions_dedup.txt
        
        # Replace original with deduplicated
        mv accessions_dedup.txt accessions.txt
        rm -f accession_pairs.tmp accessions_dedup.tmp temp_dedup.txt
        
        DEDUP_COUNT=\$(wc -l < accessions.txt | tr -d ' ')
        echo "Deduplicated: \$ASSEMBLY_COUNT -> \$DEDUP_COUNT accessions (prioritized GCF over GCA)"
    else
        DEDUP_COUNT=0
    fi

    echo "First 3 accessions (after deduplication):"
    head -3 accessions.txt || echo "(none)"

    DOWNLOAD_COUNT=\$DEDUP_COUNT
    echo "Found \$ASSEMBLY_COUNT assemblies, \$DOWNLOAD_COUNT accessions to download"
    
    # Create annotations.txt - find which accessions have GFF3 annotations
    echo "Checking which accessions have GFF3 annotations available..."
    touch annotations.txt
    
    # Use dataformat to query annotation status (available in ncbi-datasets container)
    # Try annotinfo-release-date - if it has a date, annotation exists
    if command -v dataformat >/dev/null 2>&1; then
        datasets summary genome accession --inputfile accessions.txt --as-json-lines 2>/dev/null | \\
            dataformat tsv genome --fields accession,annotinfo-release-date 2>/dev/null > annotations_check.tmp || true
        
        # Debug: show first few lines
        echo "  Annotation release date check (first 10 lines):"
        head -10 annotations_check.tmp 2>/dev/null || echo "  (empty)"
        
        # Filter to only accessions with a release date (not empty, not "na")
        awk -F'\\t' 'NR>1 && \$2 != "" && \$2 != "na" {print \$1}' annotations_check.tmp | sort -u > annotations.txt
        rm -f annotations_check.tmp
    elif command -v jq >/dev/null 2>&1; then
        jq -r '.reports[]? | select(.annotation_info != null) | .accession' search_results.json 2>/dev/null | sort -u > annotations.txt || true
    else
        echo "  Warning: Cannot determine annotations pre-download"
    fi
    
    annotations_expected=\$(wc -l < annotations.txt)
    echo "Created annotations.txt (\$annotations_expected accessions have GFF3 annotations available)"

    # Handle empty genus case
    if [ "\$DOWNLOAD_COUNT" -eq 0 ]; then
        echo "WARNING: No genomes available for ${genus} in NCBI"
        
        echo "accession,species,family,genus,completeness,contamination,quality_source,assembly_level,file_path,gff_path" > genome_metadata_download.csv
        
        printf '{
    "genus": "%s",
    "family": "%s",
    "status": "empty",
    "message": "No genomes found in NCBI for this genus",
    "assemblies_found": 0,
    "genomes_downloaded": 0,
    "species_count": 0,
    "timestamp": "%s"
}\\n' "${genus}" "${family}" "\$(date -Iseconds)" > download_info.json
        
        rm -f search_results.json accessions.txt query_error.log
        echo "════════════════════════════════════════════════════════════════"
        exit 0
    fi

    # Step 2: Create accession -> species mapping
    echo -e "\\n[2/4] Creating species mapping..."
    
    if command -v dataformat >/dev/null 2>&1; then
        datasets summary genome accession --inputfile accessions.txt --as-json-lines 2>/dev/null | \\
            dataformat tsv genome --fields accession,organism-name,assminfo-level,checkm-completeness,checkm-contamination \\
            > accession_info.tsv 2>/dev/null || touch accession_info.tsv
        echo "Created accession_info.tsv with \$(wc -l < accession_info.tsv) lines"
    else
        touch accession_info.tsv
    fi

    # Step 3: Download genomes with GFF3 annotations
    echo -e "\\n[3/4] Downloading \$DOWNLOAD_COUNT genomes (with GFF3 annotations)..."
    
    # Strategy: Binary search for problematic accessions
    # 1. Try all at once (fastest when it works)
    # 2. If fails, use binary search to find problematic accession(s)
    # 3. Binary search: O(log n) downloads to find culprit vs O(n) individual
    # 4. Exclude problematic, download rest in one batch
    
    cp accessions.txt current_accessions.txt
    touch problematic_accessions.txt
    mkdir -p ncbi_dataset/data
    
    # Function to test if a batch downloads successfully
    test_batch() {
        local batch_file="\$1"
        local zip_name="\$2"
        
        rm -f "\$zip_name"
        datasets download genome accession --inputfile "\$batch_file" --include gff3,genome --filename "\$zip_name" 2>/dev/null || true
        
        if [ -f "\$zip_name" ] && unzip -t "\$zip_name" >/dev/null 2>&1; then
            return 0
        fi
        rm -f "\$zip_name"
        return 1
    }
    
    # Function to find problematic accession(s) using binary search
    # Writes problematic accessions to stdout
    binary_search_problematic() {
        local input_file="\$1"
        local depth="\$2"
        local indent=""
        for i in \$(seq 1 \$depth); do indent="  \$indent"; done
        
        local count=\$(wc -l < "\$input_file")
        
        # Base case: single accession - it must be the problem
        if [ \$count -eq 1 ]; then
            local acc=\$(cat "\$input_file")
            echo "\$acc"
            echo "\${indent}FOUND: \$acc" >&2
            return
        fi
        
        # Base case: very small list - test individually
        if [ \$count -le 3 ]; then
            while read -r acc; do
                [ -z "\$acc" ] && continue
                echo "\$acc" > "test_single_\$depth.txt"
                if ! test_batch "test_single_\$depth.txt" "test_single_\$depth.zip"; then
                    echo "\$acc"
                    echo "\${indent}FOUND: \$acc" >&2
                fi
                rm -f "test_single_\$depth.txt" "test_single_\$depth.zip"
            done < "\$input_file"
            return
        fi
        
        echo "\${indent}Testing \$count accessions..." >&2
        
        # Split in half
        local half=\$((count / 2))
        head -n \$half "\$input_file" > "first_half_\$depth.txt"
        tail -n +\$((half + 1)) "\$input_file" > "second_half_\$depth.txt"
        
        # Test first half
        if ! test_batch "first_half_\$depth.txt" "test_first_\$depth.zip"; then
            echo "\${indent}Problem in first half (\$half accessions)" >&2
            binary_search_problematic "first_half_\$depth.txt" \$((depth + 1))
        fi
        
        # Test second half
        local second_count=\$((count - half))
        if ! test_batch "second_half_\$depth.txt" "test_second_\$depth.zip"; then
            echo "\${indent}Problem in second half (\$second_count accessions)" >&2
            binary_search_problematic "second_half_\$depth.txt" \$((depth + 1))
        fi
        
        rm -f "first_half_\$depth.txt" "second_half_\$depth.txt" "test_first_\$depth.zip" "test_second_\$depth.zip"
    }
    
    # Try downloading everything at once first
    echo "Attempting to download all \$DOWNLOAD_COUNT genomes at once..."
    
    if test_batch current_accessions.txt genomes.zip; then
        echo "Download successful! Extracting..."
        unzip -o -q genomes.zip
        rm -f genomes.zip
    else
        echo "Download failed - starting binary search for problematic accession(s)..."
        echo "(This requires ~log2(\$DOWNLOAD_COUNT) test downloads)"
        echo ""
        
        # Find problematic accessions
        problematic=\$(binary_search_problematic current_accessions.txt 0)
        
        if [ -n "\$problematic" ]; then
            echo "\$problematic" >> problematic_accessions.txt
            sort -u problematic_accessions.txt -o problematic_accessions.txt
            
            problematic_count=\$(wc -l < problematic_accessions.txt)
            echo -e "\\nFound \$problematic_count problematic accession(s):"
            cat problematic_accessions.txt
            
            # Remove problematic from list and retry
            sort current_accessions.txt > sorted_current.txt
            sort problematic_accessions.txt > sorted_problematic.txt
            comm -23 sorted_current.txt sorted_problematic.txt > clean_accessions.txt
            rm -f sorted_current.txt sorted_problematic.txt
            
            clean_count=\$(wc -l < clean_accessions.txt)
            echo -e "\\nDownloading \$clean_count clean accessions..."
            
            if test_batch clean_accessions.txt genomes.zip; then
                echo "Clean download successful! Extracting..."
                unzip -o -q genomes.zip
                rm -f genomes.zip
            else
                echo "ERROR: Download still failing after excluding problematic accessions"
                echo "This may indicate multiple problematic accessions or a different issue"
            fi
            
            rm -f clean_accessions.txt
        else
            echo "ERROR: Could not identify any problematic accessions"
            echo "The issue may be transient - consider retrying later"
        fi
    fi
    
    rm -f current_accessions.txt
    
    # Report results
    downloaded_count=\$(find ncbi_dataset/data -maxdepth 1 -type d -name "GC*" 2>/dev/null | wc -l)
    problematic_count=\$(wc -l < problematic_accessions.txt 2>/dev/null || echo 0)
    echo -e "\\nDownload summary: \$downloaded_count successful, \$problematic_count problematic"
    
    # Organize genomes by species
    if [ -d "ncbi_dataset/data" ]; then
        echo "Organizing genomes by species..."
        
        for acc_dir in ncbi_dataset/data/GC*/; do
            [ -d "\$acc_dir" ] || continue
            
            ACCESSION=\$(basename "\$acc_dir")
            GENOMIC_FILE=\$(find "\$acc_dir" -name "*.fna" -type f | head -1)
            [ -f "\$GENOMIC_FILE" ] || continue
            
            # Get species from mapping
            ORGANISM=""
            if [ -s accession_info.tsv ]; then
                ORGANISM=\$(awk -F'\\t' -v acc="\$ACCESSION" 'NR>1 && \$1 == acc {print \$2; exit}' accession_info.tsv)
                if [ -z "\$ORGANISM" ]; then
                    base_acc=\${ACCESSION%.*}
                    ORGANISM=\$(awk -F'\\t' -v acc="\$base_acc" 'NR>1 && index(\$1, acc)==1 {print \$2; exit}' accession_info.tsv)
                fi
            fi
            
            # Fallback to jq
            if [ -z "\$ORGANISM" ] && command -v jq >/dev/null 2>&1; then
                ORGANISM=\$(jq -r ".reports[]? | select(.accession == \\"\$ACCESSION\\") | .organism.organism_name // empty" search_results.json 2>/dev/null | head -1)
            fi
            
            [ -z "\$ORGANISM" ] && ORGANISM="${genus} sp."
            
            # Clean species name
            SPECIES=\$(echo "\$ORGANISM" | awk '{print \$1"_"\$2}' | sed 's/[^a-zA-Z0-9_]/_/g; s/__*/_/g; s/_\$//')
            [ -z "\$SPECIES" ] || [ "\$SPECIES" = "_" ] && SPECIES="${genus}_sp"
            
            mkdir -p "\$SPECIES"
            cp "\$GENOMIC_FILE" "\$SPECIES/\${ACCESSION}.fna"
            
            # Copy GFF if exists (NCBI downloads as .gff)
            GFF_FILE=\$(find "\$acc_dir" -type f -name "*.gff" | head -1)
            if [ -f "\$GFF_FILE" ]; then
                cp "\$GFF_FILE" "\$SPECIES/\${ACCESSION}.gff"
                GFF_COUNT=\$((GFF_COUNT + 1))
            fi
            
            DOWNLOADED_COUNT=\$((DOWNLOADED_COUNT + 1))
            
            if [ \$((DOWNLOADED_COUNT % 100)) -eq 0 ]; then
                echo "  Processed \$DOWNLOADED_COUNT genomes..."
            fi
        done
        
        rm -rf ncbi_dataset genomes.zip
    else
        echo "ERROR: No genome data available"
    fi

    echo "Organized \$DOWNLOADED_COUNT genomes with \$GFF_COUNT GFF annotations"

    # Update annotations.txt with actual GFF files found (overwrite the pre-download estimate)
    find . -maxdepth 2 -name "*.gff" -type f -printf "%f\\n" | sed 's/\\.gff\$//' | sort > annotations.txt
    echo "Updated annotations.txt with \$(wc -l < annotations.txt) actual GFF files"

    # Step 4: Create metadata CSV
    echo -e "\\n[4/4] Creating metadata..."

    SPECIES_COUNT=\$(find . -maxdepth 1 -type d ! -name "." 2>/dev/null | wc -l)

    # Create metadata header
    echo "accession,species,family,genus,completeness,contamination,quality_source,assembly_level,file_path,gff_path" > genome_metadata_download.csv

    # Find all fasta files and create metadata entries
    # Use a temp file to avoid subshell issues with the while loop
    find . -maxdepth 2 -name "*.fna" -type f > fasta_list.txt
    
    while read -r fasta; do
        acc=\$(basename "\$fasta" .fna)
        species_dir=\$(dirname "\$fasta" | sed 's|^[.]/||')
        
        # Check for GFF file
        gff_path=""
        if [ -f "\${species_dir}/\${acc}.gff" ]; then
            gff_path="\${species_dir}/\${acc}.gff"
        fi
        
        # Get quality info from accession_info.tsv
        completeness=""
        contamination=""
        assembly_level=""
        quality_source="pending"
        species_name="\$species_dir"
        
        if [ -s accession_info.tsv ]; then
            # Try exact match first
            info_line=\$(awk -F'\\t' -v acc="\$acc" 'NR>1 && \$1 == acc {print; exit}' accession_info.tsv)
            # Try without version number
            if [ -z "\$info_line" ]; then
                base_acc=\${acc%.*}
                info_line=\$(awk -F'\\t' -v acc="\$base_acc" 'NR>1 && index(\$1, acc)==1 {print; exit}' accession_info.tsv)
            fi
            
            if [ -n "\$info_line" ]; then
                species_name=\$(echo "\$info_line" | cut -f2 | sed 's/,/_/g')
                assembly_level=\$(echo "\$info_line" | cut -f3)
                completeness=\$(echo "\$info_line" | cut -f4)
                contamination=\$(echo "\$info_line" | cut -f5)
                
                # Check if we have valid completeness from NCBI
                if [ -n "\$completeness" ] && [ "\$completeness" != "na" ] && [ "\$completeness" != "" ]; then
                    quality_source="NCBI"
                    # If contamination is empty but completeness exists, treat as 0
                    # NCBI web shows "0%" but API returns empty for zero contamination
                    if [ -z "\$contamination" ] || [ "\$contamination" = "na" ]; then
                        contamination="0"
                    fi
                else
                    completeness=""
                    contamination=""
                fi
            fi
        fi
        
        echo "\$acc,\$species_name,${family},${genus},\$completeness,\$contamination,\$quality_source,\$assembly_level,\$fasta,\$gff_path"
    done < fasta_list.txt >> genome_metadata_download.csv
    
    rm -f fasta_list.txt

    # Count quality sources (ensure clean integers)
    genome_count=\$(tail -n +2 genome_metadata_download.csv | wc -l | tr -d '[:space:]')
    ncbi_quality=\$(grep -c ",NCBI," genome_metadata_download.csv 2>/dev/null || true)
    ncbi_quality=\${ncbi_quality:-0}
    [ -z "\$ncbi_quality" ] && ncbi_quality=0
    pending_quality=\$(grep -c ",pending," genome_metadata_download.csv 2>/dev/null || true)
    pending_quality=\${pending_quality:-0}
    [ -z "\$pending_quality" ] && pending_quality=0
    gff_in_metadata=\$(awk -F',' 'NR>1 && \$10 != "" {count++} END {print count+0}' genome_metadata_download.csv)
    problematic_count=0
    [ -f problematic_accessions.txt ] && problematic_count=\$(wc -l < problematic_accessions.txt | tr -d '[:space:]')

    # Determine status
    if [ \$DOWNLOADED_COUNT -gt 0 ]; then
        status="success"
    else
        status="failed"
    fi

    # Create info JSON using printf for robustness
    printf '{
    "genus": "%s",
    "family": "%s",
    "status": "%s",
    "assemblies_found": %d,
    "genomes_downloaded": %d,
    "gff_annotations": %d,
    "species_count": %d,
    "ncbi_quality_available": %d,
    "pending_checkm2": %d,
    "problematic_accessions": %d,
    "timestamp": "%s"
}\\n' \\
        "${genus}" \\
        "${family}" \\
        "\$status" \\
        "\$ASSEMBLY_COUNT" \\
        "\$DOWNLOADED_COUNT" \\
        "\$GFF_COUNT" \\
        "\$SPECIES_COUNT" \\
        "\$ncbi_quality" \\
        "\$pending_quality" \\
        "\$problematic_count" \\
        "\$(date -Iseconds)" > download_info.json

    # Cleanup
    rm -f search_results.json accessions.txt accession_info.tsv query_error.log

    echo -e "\\n════════════════════════════════════════════════════════════════"
    echo "COMPLETED: \$DOWNLOADED_COUNT genomes in \$SPECIES_COUNT species"
    echo "  - With GFF3 annotations: \$GFF_COUNT"
    echo "  - With NCBI quality: \$ncbi_quality"
    echo "  - Pending CheckM2: \$pending_quality"
    [ \$problematic_count -gt 0 ] && echo "  - Problematic accessions: \$problematic_count (see problematic_accessions.txt)"
    echo "════════════════════════════════════════════════════════════════"
    """
}
