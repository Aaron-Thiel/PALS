/*
 * SUMMARY - Generate pipeline summary TSV
 *
 * Collects key metrics from all pipeline stages into a single file:
 * - Sample name
 * - Taxonomy (family/genus/species)
 * - SPAdes QC: completeness, contamination, genome_size
 * - PASA QC: completeness, contamination, genome_size (for comparison)
 * - Gene counts from Bakta
 * - Best reference genome info
 */

process SUMMARY {
    publishDir "${params.outdir}", mode: 'copy'

    container 'python:3.11'

    input:
    path(spades_qc_files, stageAs: 'spades_qc/s?/*')
    path(pasa_qc_files, stageAs: 'pasa_qc/p?/*')
    path(bakta_files, stageAs: 'bakta/b?/*')
    path(refs_files, stageAs: 'refs/r?/*')
    val(taxonomy_data)  // List of maps: [sample_id: X, genus: Y, species: Z, family: W]

    output:
    path("pipeline_summary.tsv"), emit: summary

    script:
    // Convert taxonomy data to JSON
    def taxonomy_json = new groovy.json.JsonBuilder(taxonomy_data).toString()
    """
    #!/usr/bin/env python3
    import json
    import csv
    import os
    import glob

    # Parse taxonomy data
    taxonomy_list = json.loads('${taxonomy_json}')
    taxonomy_map = {t['sample_id']: t for t in taxonomy_list if t.get('sample_id')}

    # Collect all sample IDs from available files
    sample_ids = set()

    # Get sample IDs from SPAdes QC files (recursive glob for indexed subdirs)
    for f in glob.glob('spades_qc/**/*.json', recursive=True):
        try:
            with open(f, 'r') as fp:
                data = json.load(fp)
                if 'sample_id' in data:
                    sample_ids.add(data['sample_id'])
        except:
            pass

    # Also check PASA QC files
    for f in glob.glob('pasa_qc/**/*.json', recursive=True):
        try:
            with open(f, 'r') as fp:
                data = json.load(fp)
                if 'sample_id' in data:
                    sample_ids.add(data['sample_id'])
        except:
            pass

    # Add from taxonomy
    for sid in taxonomy_map.keys():
        sample_ids.add(sid)

    print(f"Found {len(sample_ids)} samples to summarize")

    # TSV header with SPAdes and PASA comparison columns
    header = [
        "sample_id",
        "family",
        "genus",
        "species",
        "spades_completeness",
        "spades_contamination",
        "spades_genome_size",
        "pasa_completeness",
        "pasa_contamination",
        "pasa_genome_size",
        "completeness_change",
        "contamination_change",
        "total_genes",
        "cds_count",
        "best_reference",
        "best_reference_ani"
    ]

    # Helper function to find JSON file for a sample (recursive for indexed subdirs)
    def find_sample_json(directory, sample_id):
        for f in glob.glob(f'{directory}/**/*.json', recursive=True):
            try:
                with open(f, 'r') as fp:
                    data = json.load(fp)
                    if data.get('sample_id') == sample_id:
                        return data
            except:
                continue
        return None

    # Helper to find best PASA QC (highest completeness) for a sample
    # This handles cases where both standard and sensitive QC files exist
    def find_best_pasa_qc(directory, sample_id):
        best_data = None
        best_completeness = -1
        for f in glob.glob(f'{directory}/**/*.json', recursive=True):
            try:
                with open(f, 'r') as fp:
                    data = json.load(fp)
                    if data.get('sample_id') == sample_id:
                        completeness = float(data.get('completeness', 0))
                        if completeness > best_completeness:
                            best_completeness = completeness
                            best_data = data
            except:
                continue
        return best_data

    # Build reference lookup from all CSV files (sample_id -> best reference)
    # Each CSV now has sample_id column for matching
    def build_refs_lookup(directory):
        refs_lookup = {}
        for f in glob.glob(f'{directory}/**/*.csv', recursive=True):
            try:
                with open(f, 'r') as fp:
                    reader = csv.DictReader(fp)
                    for row in reader:
                        sid = row.get('sample_id')
                        if sid and sid not in refs_lookup:
                            # Take first (best) reference for each sample
                            refs_lookup[sid] = {
                                'accession': row.get('accession', 'NA'),
                                'ani': row.get('ani', 'NA')
                            }
            except:
                continue
        return refs_lookup

    # Build reference lookup once
    refs_lookup = build_refs_lookup('refs')

    rows = []
    for sample_id in sorted(sample_ids):
        row = {
            "sample_id": sample_id,
            "family": "NA",
            "genus": "NA",
            "species": "NA",
            "spades_completeness": "NA",
            "spades_contamination": "NA",
            "spades_genome_size": "NA",
            "pasa_completeness": "NA",
            "pasa_contamination": "NA",
            "pasa_genome_size": "NA",
            "completeness_change": "NA",
            "contamination_change": "NA",
            "total_genes": "NA",
            "cds_count": "NA",
            "best_reference": "NA",
            "best_reference_ani": "NA"
        }

        # Get taxonomy
        if sample_id in taxonomy_map:
            tax = taxonomy_map[sample_id]
            row["family"] = tax.get('family', 'NA')
            row["genus"] = tax.get('genus', 'NA')
            row["species"] = tax.get('species', 'NA')

        # Get SPAdes QC
        spades_data = find_sample_json('spades_qc', sample_id)
        if spades_data:
            row["spades_completeness"] = spades_data.get('completeness', 'NA')
            row["spades_contamination"] = spades_data.get('contamination', 'NA')
            row["spades_genome_size"] = spades_data.get('genome_size', 'NA')

        # Get PASA QC (best result if both standard and sensitive exist)
        pasa_data = find_best_pasa_qc('pasa_qc', sample_id)
        if pasa_data:
            row["pasa_completeness"] = pasa_data.get('completeness', 'NA')
            row["pasa_contamination"] = pasa_data.get('contamination', 'NA')
            row["pasa_genome_size"] = pasa_data.get('genome_size', 'NA')

        # Calculate changes (PASA - SPAdes)
        try:
            if row["spades_completeness"] != "NA" and row["pasa_completeness"] != "NA":
                change = float(row["pasa_completeness"]) - float(row["spades_completeness"])
                row["completeness_change"] = f"{change:+.2f}"
        except:
            pass

        try:
            if row["spades_contamination"] != "NA" and row["pasa_contamination"] != "NA":
                change = float(row["pasa_contamination"]) - float(row["spades_contamination"])
                row["contamination_change"] = f"{change:+.2f}"
        except:
            pass

        # Get Bakta gene counts
        bakta_data = find_sample_json('bakta', sample_id)
        if bakta_data and 'annotation_stats' in bakta_data:
            stats = bakta_data['annotation_stats']
            row["total_genes"] = stats.get('total_features', 'NA')
            counts = stats.get('feature_counts', {})
            # CDS count - try different possible keys
            cds = counts.get('cds') or counts.get('CDS') or counts.get('gene')
            if cds:
                row["cds_count"] = cds

        # Get reference info from pre-built lookup
        if sample_id in refs_lookup:
            refs_data = refs_lookup[sample_id]
            row["best_reference"] = refs_data.get('accession', 'NA')
            row["best_reference_ani"] = refs_data.get('ani', 'NA')

        rows.append(row)

    # Write TSV
    with open("pipeline_summary.tsv", 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter='\\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Generated summary for {len(rows)} samples")
    print(f"Columns: {', '.join(header)}")
    """
}
