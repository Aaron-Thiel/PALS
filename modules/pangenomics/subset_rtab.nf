/*
 * Subset Rtab by Genus/Species
 *
 * Takes the full pangenome Rtab file and creates per-genus and per-species
 * subsets based on taxonomy. Only creates subsets for genera/species that
 * have at least one internal (PAL*) sample.
 *
 * Uses the Rtab (binary 0/1 matrix) which is the correct input format
 * for PANGENOME_PLOTS.
 */

process SUBSET_RTAB {
    tag "$cohort_id"
    label 'process_medium'

    publishDir "${params.outdir}/pangenomics/subsets", mode: 'copy', pattern: 'manifest.tsv'

    container 'python:3.11'

    input:
    tuple val(cohort_id), path(rtab_file)
    path(pipeline_summary)
    path(taxonomy_csv)

    output:
    path("subsets/*.Rtab"), emit: rtab_files
    path("manifest.tsv"), emit: manifest

    script:
    def min_samples = params.subset_min_samples_pangenomics ?: 3
    """
    #!/usr/bin/env python3

    import csv
    import os

    print("=" * 60)
    print("Subset Rtab by Genus/Species")
    print("=" * 60)
    print(f"Rtab file: ${rtab_file}")
    print(f"Pipeline summary: ${pipeline_summary}")
    print(f"Taxonomy CSV: ${taxonomy_csv}")
    print(f"Min samples per group: ${min_samples}")
    print("")

    min_samples = int(${min_samples})

    # ================================================================
    # Step 1: Read Rtab header to get sample names
    # ================================================================
    with open("${rtab_file}") as f:
        header = f.readline().strip().split("\\t")

    gene_col = header[0]  # "Gene"
    sample_names = header[1:]
    print(f"Total samples in Rtab: {len(sample_names)}")

    # ================================================================
    # Step 2: Build taxonomy mapping from pipeline_summary.tsv (internal)
    # ================================================================
    internal_taxonomy = {}  # sample_id -> (genus, species)
    internal_genera = set()
    internal_species = set()

    with open("${pipeline_summary}") as f:
        reader = csv.DictReader(f, delimiter="\\t")
        for row in reader:
            sid = row['sample_id']
            genus = row.get('genus', '')
            species = row.get('species', '')
            if genus and genus != 'NA' and genus != '':
                internal_taxonomy[sid] = (genus, species if species and species != 'NA' else '')
                internal_genera.add(genus)
                if species and species != 'NA':
                    internal_species.add(species)

    print(f"Internal samples with taxonomy: {len(internal_taxonomy)}")
    print(f"Internal genera: {sorted(internal_genera)}")
    print(f"Internal species: {len(internal_species)}")

    # ================================================================
    # Step 3: Build taxonomy mapping from taxonomy.csv (external)
    # ================================================================
    external_taxonomy = {}  # genome_id -> (genus, species)

    with open("${taxonomy_csv}") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gid = row['genome_id']
            genus = row.get('genus', '')
            species = row.get('species', '')
            if genus and genus != '':
                # Normalize species: replace spaces with underscores
                if species:
                    species = species.replace(' ', '_')
                external_taxonomy[gid] = (genus, species if species else '')

    print(f"External genomes with taxonomy: {len(external_taxonomy)}")

    # ================================================================
    # Step 4: Map each sample in Rtab to genus/species
    # ================================================================
    sample_genus = {}
    sample_species = {}

    for sample in sample_names:
        if sample in internal_taxonomy:
            genus, species = internal_taxonomy[sample]
            sample_genus[sample] = genus
            if species:
                sample_species[sample] = species
        elif sample in external_taxonomy:
            genus, species = external_taxonomy[sample]
            sample_genus[sample] = genus
            if species:
                sample_species[sample] = species

    print(f"Samples mapped to genus: {len(sample_genus)}")
    print(f"Samples mapped to species: {len(sample_species)}")

    # ================================================================
    # Step 5: Build groups (only for genera/species in internal samples)
    # ================================================================
    groups = {}  # group_id -> list of (sample_name, column_index)

    # Build column index map
    col_index = {name: i for i, name in enumerate(sample_names)}

    # Genus-level groups
    for genus in sorted(internal_genera):
        group_id = f"genus_{genus}"
        members = [(s, col_index[s]) for s in sample_names
                   if sample_genus.get(s) == genus]
        n_internal = sum(1 for s, _ in members if s in internal_taxonomy)
        if len(members) >= min_samples:
            groups[group_id] = {
                'members': members,
                'level': 'genus',
                'name': genus,
                'n_samples': len(members),
                'n_internal': n_internal
            }
            print(f"  {group_id}: {len(members)} samples ({n_internal} internal)")
        else:
            print(f"  {group_id}: SKIPPED ({len(members)} < {min_samples} min)")

    # Species-level groups
    for species in sorted(internal_species):
        group_id = f"species_{species}"
        members = [(s, col_index[s]) for s in sample_names
                   if sample_species.get(s) == species]
        n_internal = sum(1 for s, _ in members if s in internal_taxonomy)
        if len(members) >= min_samples:
            groups[group_id] = {
                'members': members,
                'level': 'species',
                'name': species,
                'n_samples': len(members),
                'n_internal': n_internal
            }
            print(f"  {group_id}: {len(members)} samples ({n_internal} internal)")
        else:
            print(f"  {group_id}: SKIPPED ({len(members)} < {min_samples} min)")

    print(f"\\nTotal groups to create: {len(groups)}")

    if not groups:
        print("WARNING: No qualifying groups found. Creating empty manifest.")
        os.makedirs("subsets", exist_ok=True)
        # Create a dummy Rtab so output glob doesn't fail
        with open("subsets/empty.Rtab", 'w') as f:
            f.write("Gene\\n")
        with open("manifest.tsv", 'w') as f:
            f.write("group_id\\tlevel\\tname\\tn_samples\\tn_internal\\n")
        exit(0)

    # ================================================================
    # Step 6: Read full Rtab and write subsets
    # ================================================================
    os.makedirs("subsets", exist_ok=True)

    print("\\nReading full Rtab and writing subsets...")

    # Pre-compute column indices for each group
    group_cols = {}
    for group_id, info in groups.items():
        group_cols[group_id] = [idx for _, idx in info['members']]

    # Open all output files
    out_files = {}
    for group_id, info in groups.items():
        out_path = f"subsets/{group_id}.Rtab"
        out_files[group_id] = open(out_path, 'w')
        # Write header
        member_names = [name for name, _ in info['members']]
        out_files[group_id].write(gene_col + "\\t" + "\\t".join(member_names) + "\\n")

    # Stream through Rtab rows
    row_count = 0
    with open("${rtab_file}") as f:
        next(f)  # skip header
        for line in f:
            row_count += 1
            parts = line.strip().split("\\t")
            gene_name = parts[0]
            values = parts[1:]  # 0/1 values

            for group_id, cols in group_cols.items():
                subset_values = [values[c] for c in cols]
                # Only write if at least one sample has the gene
                if any(v != '0' for v in subset_values):
                    out_files[group_id].write(gene_name + "\\t" + "\\t".join(subset_values) + "\\n")

            if row_count % 50000 == 0:
                print(f"  Processed {row_count} genes...")

    # Close all files
    for f in out_files.values():
        f.close()

    print(f"Processed {row_count} total genes")

    # ================================================================
    # Step 7: Write manifest
    # ================================================================
    with open("manifest.tsv", 'w') as f:
        f.write("group_id\\tlevel\\tname\\tn_samples\\tn_internal\\n")
        for group_id, info in sorted(groups.items()):
            f.write(f"{group_id}\\t{info['level']}\\t{info['name']}\\t{info['n_samples']}\\t{info['n_internal']}\\n")

    print("\\nSubset Rtab creation completed!")
    for group_id, info in sorted(groups.items()):
        print(f"  {group_id}: {info['n_samples']} samples")
    """
}
