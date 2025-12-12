/*
 * MERGE_GFF_FASTA - Merge GFF and FASTA files for panta compatibility
 *
 * Takes downloaded GFF files and their corresponding FNA files,
 * merges them into proper GFF3 format with embedded ##FASTA section.
 * Updates genome_metadata.csv to point to the new .gff3 files.
 */
process MERGE_GFF_FASTA {
    tag "${genus}"

    storeDir "${params.genome_database}/${family}/${genus}"

    container 'python:3.11'

    input:
    tuple val(genus), val(family), path("genome_metadata_download.csv")

    output:
    tuple val(genus), val(family), path("genome_metadata_merged.csv"), emit: metadata

    script:
    """
    #!/usr/bin/env python3
    import csv
    from pathlib import Path

    GENUS_DIR = Path("${params.genome_database}/${family}/${genus}")

    # Read input metadata
    with open("genome_metadata_download.csv") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames

    print(f"Processing {len(rows)} genomes")
    merged_count = 0

    for row in rows:
        accession = row['accession']
        gff_path = row.get('gff_path', '')

        if not gff_path or not gff_path.endswith('.gff'):
            continue

        gff_file = GENUS_DIR / gff_path
        fna_path = gff_path.replace('.gff', '.fna')
        fna_file = GENUS_DIR / fna_path
        gff3_path = gff_path.replace('.gff', '.gff3')
        gff3_file = GENUS_DIR / gff3_path

        if not gff_file.exists() or not fna_file.exists():
            continue

        # Merge GFF + FNA into GFF3
        content = gff_file.read_text().rstrip()
        if content.endswith('###'):
            content = content[:-3].rstrip()

        gff3_file.write_text(content + '\\n##FASTA\\n' + fna_file.read_text())
        row['gff_path'] = gff3_path
        merged_count += 1

    # Write updated metadata
    with open("genome_metadata_merged.csv", "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Merged {merged_count} GFF+FASTA pairs into .gff3 files")
    """
}
