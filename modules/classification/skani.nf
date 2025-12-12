/*
 * SKANI Module for Nextflow - v0.3.x compatible
 * 
 * SKANI_BUILD_DB - Build sketch database for a genus (cached via storeDir)
 * SKANI_SEARCH - Compare sample against genus database (per-sample)
 * 
 * Workflow: skani sketch → skani search
 * Compatible with classification.nf and reference_selector.nf
 */

process SKANI_BUILD_DB {
    tag "${genus}"
    
    // Cache the entire database directory
    storeDir "${params.genome_database}/${family}/${genus}/skani_db"

    container 'staphb/skani:latest'  // v0.3.x

    input:
    tuple val(genus), val(family), path(metadata)

    output:
    tuple val(genus), val(family), path("sketch_db"), emit: database
    tuple val(genus), val(family), path("genome_count.txt"), emit: stats

    script:
    // Presets for fragmented/contig references:
    // --medium: -c 70, good for ANI ≤95% and N50 ≤10 kb
    // --slow: -c 30, for very fragmented genomes with N50 ~3kb  
    // --small-genomes: -c 30 -m 200 --faster-small, for small genomes/contigs
    def sketch_preset = params.skani_preset ?: '--medium'
    """
    #!/bin/bash
    set -euo pipefail

    echo "════════════════════════════════════════════════════════════════"
    echo "SKANI DATABASE BUILD: ${genus}"
    echo "skani version: \$(skani --version 2>&1 || echo 'unknown')"
    echo "════════════════════════════════════════════════════════════════"

    GENUS_DIR="${params.genome_database}/${family}/${genus}"

    # Find all reference genomes (support multiple extensions)
    find "\$GENUS_DIR" \\( -name "*.fna" -o -name "*.fna.gz" -o -name "*.fasta" -o -name "*.fa" \\) -type f > genome_list.txt 2>/dev/null || touch genome_list.txt
    genome_count=\$(wc -l < genome_list.txt | tr -d ' ')
    
    echo "Found \$genome_count genomes in \$GENUS_DIR"
    echo "\$genome_count" > genome_count.txt

    if [ "\$genome_count" -eq 0 ]; then
        echo "WARNING: No genomes found - creating empty database marker"
        mkdir -p sketch_db
        echo "empty_db" > sketch_db/empty_marker
        exit 0
    fi

    # Build sketch database using skani sketch
    # In v0.3.x, this creates a unified database
    echo "Building sketch database with preset: ${sketch_preset}"
    
    skani sketch \\
        -l genome_list.txt \\
        -o sketch_db \\
        -t ${task.cpus} \\
        ${sketch_preset}

    # Verify database was created
    if [ ! -f "sketch_db/markers.bin" ]; then
        echo "ERROR: Database creation failed - markers.bin not found"
        ls -la sketch_db/ 2>/dev/null || echo "sketch_db directory doesn't exist"
        exit 1
    fi

    # Store metadata
    cat > sketch_db/db_info.txt << EOF
genus: ${genus}
family: ${family}
genome_count: \$genome_count
created: \$(date -Iseconds)
skani_version: \$(skani --version 2>&1 || echo "unknown")
preset: ${sketch_preset}
EOF

    echo "✓ Database created successfully"
    echo "  - Genomes indexed: \$genome_count"
    ls -lh sketch_db/
    
    rm -f genome_list.txt

    echo "════════════════════════════════════════════════════════════════"
    """
}


process SKANI_SEARCH {
    tag "${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}/classification/skani", mode: 'copy'

    container 'staphb/skani:latest'  // v0.3.x

    input:
    tuple val(sample_id), path(query_fasta), val(genus), val(family), path(sketch_db)

    output:
    tuple val(sample_id), val(genus), path("ani_results.tsv"), emit: ani_results
    tuple val(sample_id), val(genus), path("species_assignment.txt"), emit: species
    tuple val(sample_id), path("skani_summary.txt"), emit: summary

    script:
    // Note: search inherits parameters from the sketch database
    // Presets like --medium are only valid for sketch/dist, not search
    """
    #!/bin/bash
    set -eu

    echo "════════════════════════════════════════════════════════════════"
    echo "SKANI SEARCH: ${sample_id} vs ${genus}"
    echo "════════════════════════════════════════════════════════════════"

    # Validate database
    if [ ! -d "${sketch_db}" ]; then
        echo "ERROR: Database directory not found: ${sketch_db}"
        echo "sp_nov_${genus}" > species_assignment.txt
        printf "Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name\\n" > ani_results.tsv
        echo "ERROR: Database not found" > skani_summary.txt
        exit 0
    fi

    # Check for empty database marker
    if [ -f "${sketch_db}/empty_marker" ]; then
        echo "WARNING: Empty database (no reference genomes for ${genus})"
        echo "sp_nov_${genus}" > species_assignment.txt
        printf "Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name\\n" > ani_results.tsv
        cat > skani_summary.txt << EOF
Sample: ${sample_id}
Genus: ${genus}
Status: No reference genomes available
Assignment: sp_nov_${genus} (novel species - no references)
EOF
        exit 0
    fi

    # Validate database structure
    if [ ! -f "${sketch_db}/markers.bin" ]; then
        echo "ERROR: Invalid database - markers.bin not found"
        ls -la "${sketch_db}/" || true
        echo "sp_nov_${genus}" > species_assignment.txt
        printf "Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name\\n" > ani_results.tsv
        echo "ERROR: Invalid database structure" > skani_summary.txt
        exit 0
    fi

    echo "Database: ${sketch_db}"
    echo "Query: ${query_fasta}"

    # Count query contigs
    query_contigs=\$(grep -c "^>" ${query_fasta} || echo "0")
    echo "Query contains \$query_contigs contigs"

    # Run skani search
    # Note: parameters like -c are inherited from the sketch database
    # Use --min-af 0 to return ALL genus-level matches, letting REFERENCE_SELECTOR do filtering
    echo "Running skani search..."

    if ! skani search \\
        ${query_fasta} \\
        -d "${sketch_db}" \\
        -o ani_results.tsv \\
        -t ${task.cpus} \\
        -n 1000 \\
        --min-af 0 \\
        2>&1; then
        
        echo "WARNING: skani search failed"
        echo "sp_nov_${genus}" > species_assignment.txt
        printf "Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name\\n" > ani_results.tsv
        echo "skani search failed" > skani_summary.txt
        exit 0
    fi

    echo "Search completed"

    # Check if we have results
    if [ ! -s ani_results.tsv ]; then
        echo "No results returned"
        echo "sp_nov_${genus}" > species_assignment.txt
        printf "Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name\\n" > ani_results.tsv
        echo "No matches found" > skani_summary.txt
        exit 0
    fi

    result_lines=\$(wc -l < ani_results.tsv)
    data_lines=\$((result_lines - 1))
    
    echo "Found \$data_lines hits"

    if [ "\$data_lines" -gt 0 ]; then
        # Sort by ANI (column 3) descending, get best match
        tail -n +2 ani_results.tsv | sort -t\$'\\t' -k3 -nr > sorted_results.tmp
        best=\$(head -1 sorted_results.tmp)
        
        # Parse best match fields (7 columns: Ref_file, Query_file, ANI, AF_ref, AF_query, Ref_name, Query_name)
        ref_file=\$(echo "\$best" | cut -f1)
        ani=\$(echo "\$best" | cut -f3)
        af_ref=\$(echo "\$best" | cut -f4)
        af_query=\$(echo "\$best" | cut -f5)
        
        # Extract species from path: .../Species/accession.fna
        species=\$(echo "\$ref_file" | awk -F'/' '{print \$(NF-1)}')
        [ -z "\$species" ] || [ "\$species" = "\$ref_file" ] && species="${genus}_sp"

        echo "Best: \$ref_file (ANI: \${ani}%, AF: \$af_query)"

        # Species assignment based on ANI thresholds
        assigned=""
        confidence=""
        
        if awk "BEGIN {exit !(\$ani >= 95.0 && \$af_query >= 0.6)}"; then
            assigned="\$species"
            confidence="high"
        elif awk "BEGIN {exit !(\$ani >= 95.0 && \$af_query >= 0.3)}"; then
            assigned="\$species"
            confidence="medium"
        elif awk "BEGIN {exit !(\$ani >= 90.0 && \$af_query >= 0.4)}"; then
            assigned="cf_\$species"
            confidence="medium"
        elif awk "BEGIN {exit !(\$ani >= 85.0)}"; then
            assigned="aff_\$species"
            confidence="low"
        else
            assigned="sp_nov_${genus}"
            confidence="novel"
        fi

        echo "\$assigned" > species_assignment.txt

        # Build summary
        cat > skani_summary.txt << EOF
═══════════════════════════════════════════════════════════════════
SKANI ANALYSIS: ${sample_id}
═══════════════════════════════════════════════════════════════════
Query: ${query_fasta} (\$query_contigs contigs)
Target genus: ${genus}

BEST MATCH
───────────────────────────────────────────────────────────────────
Reference: \$ref_file
Species: \$species
ANI: \${ani}%
Aligned fraction (ref): \$af_ref
Aligned fraction (query): \$af_query

ASSIGNMENT
───────────────────────────────────────────────────────────────────
Species: \$assigned
Confidence: \$confidence

TOP HITS (\$data_lines total)
───────────────────────────────────────────────────────────────────
EOF
        head -5 sorted_results.tmp | while IFS=\$'\\t' read -r rf qf a afr afq rn qn; do
            echo "  \$(basename \$rf): \${a}% ANI, \$afq AF" >> skani_summary.txt
        done
        echo "═══════════════════════════════════════════════════════════════════" >> skani_summary.txt

        rm -f sorted_results.tmp
        echo "✓ Assigned: \$assigned (\$confidence)"

    else
        echo "No significant matches"
        echo "sp_nov_${genus}" > species_assignment.txt
        cat > skani_summary.txt << EOF
═══════════════════════════════════════════════════════════════════
SKANI ANALYSIS: ${sample_id}
═══════════════════════════════════════════════════════════════════
Query: ${query_fasta}
Target genus: ${genus}

RESULT: No significant matches found
Assignment: sp_nov_${genus} (potential novel species)
═══════════════════════════════════════════════════════════════════
EOF
    fi

    echo "════════════════════════════════════════════════════════════════"
    """
}
