/*
 * Subset Core Genes by Genus/Species
 *
 * Takes the full PANTA output and extracts core genes for each genus/species
 * group that has at least one internal (PAL*) sample. Reads the CSV file
 * (not Rtab) because it contains per-sample sequence IDs needed for
 * extracting nucleotide sequences.
 *
 * Outputs one core_genes directory and genus_map.tsv per qualifying group,
 * plus a manifest.tsv for downstream channel construction.
 */

process SUBSET_CORE_GENES {
    tag "$cohort_id"
    label 'process_high'

    publishDir "${params.outdir}/phylogenetics/subsets", mode: 'copy', pattern: 'manifest.tsv'

    container 'python:3.9'

    input:
    tuple val(cohort_id), path(panta_dir)
    path(pipeline_summary)
    path(taxonomy_csv)

    output:
    path("groups/*/core_genes"), emit: core_genes_dirs, optional: true
    path("groups/*/genus_map.tsv"), emit: genus_maps, optional: true
    path("manifest.tsv"), emit: manifest

    script:
    def min_samples = params.subset_min_samples_phylo ?: 4
    def core_threshold = params.phylo_core_threshold ?: 0.95
    """
    #!/usr/bin/env python3

    import csv
    import os
    import math
    from pathlib import Path
    from collections import defaultdict

    print("=" * 60)
    print("Subset Core Genes by Genus/Species")
    print("=" * 60)
    print(f"PANTA dir: ${panta_dir}")
    print(f"Min samples per group: ${min_samples}")
    print(f"Core threshold: ${core_threshold}")
    print("")

    min_samples = int(${min_samples})
    core_threshold = float(${core_threshold})

    # ================================================================
    # Step 1: Build taxonomy mappings
    # ================================================================
    internal_taxonomy = {}
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

    external_taxonomy = {}
    with open("${taxonomy_csv}") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gid = row['genome_id']
            genus = row.get('genus', '')
            species = row.get('species', '')
            if genus and genus != '':
                if species:
                    species = species.replace(' ', '_')
                external_taxonomy[gid] = (genus, species if species else '')

    print(f"Internal samples: {len(internal_taxonomy)}, genera: {len(internal_genera)}, species: {len(internal_species)}")
    print(f"External genomes: {len(external_taxonomy)}")

    # ================================================================
    # Step 2: Read gene_presence_absence.csv header to identify samples
    # ================================================================
    gpa_path = Path("${panta_dir}") / "gene_presence_absence.csv"
    if not gpa_path.exists():
        print(f"ERROR: {gpa_path} not found")
        # Write empty manifest
        with open("manifest.tsv", 'w') as f:
            f.write("group_id\\tsample_count\\tn_core_genes\\n")
        exit(0)

    # Parse CSV header
    with open(gpa_path) as f:
        reader = csv.reader(f)
        header = next(reader)
        header = [h.strip('"').strip() for h in header]

    # Find sample columns (skip metadata columns)
    metadata_cols = {
        'Gene', 'Annotation', 'Non-unique Gene name', 'No. isolates',
        'No. sequences', 'No. sequences per isolate', 'Avg sequences per isolate',
        'Genome Fragment', 'Order within Fragment', 'Accessory Fragment',
        'Accessory Order with Fragment', 'QC', 'Min group size nuc',
        'Max group size nuc', 'Avg group size nuc', 'gene_id', 'cluster_id'
    }

    sample_start = 0
    for i, col in enumerate(header):
        if col not in metadata_cols:
            sample_start = i
            break

    all_sample_names = header[sample_start:]
    print(f"Total samples in CSV: {len(all_sample_names)}")

    # Map samples to genus/species
    sample_genus = {}
    sample_species = {}
    for sample in all_sample_names:
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

    # ================================================================
    # Step 3: Build groups
    # ================================================================
    groups = {}

    for genus in sorted(internal_genera):
        group_id = f"genus_{genus}"
        members = [s for s in all_sample_names if sample_genus.get(s) == genus]
        if len(members) >= min_samples:
            groups[group_id] = {
                'members': members,
                'level': 'genus',
                'name': genus
            }
            print(f"  {group_id}: {len(members)} samples")
        else:
            print(f"  {group_id}: SKIPPED ({len(members)} < {min_samples})")

    for species in sorted(internal_species):
        group_id = f"species_{species}"
        members = [s for s in all_sample_names if sample_species.get(s) == species]
        if len(members) >= min_samples:
            groups[group_id] = {
                'members': members,
                'level': 'species',
                'name': species
            }
            print(f"  {group_id}: {len(members)} samples")
        else:
            print(f"  {group_id}: SKIPPED ({len(members)} < {min_samples})")

    print(f"\\nTotal groups: {len(groups)}")

    if not groups:
        print("No qualifying groups. Writing empty manifest.")
        with open("manifest.tsv", 'w') as f:
            f.write("group_id\\tsample_count\\tn_core_genes\\n")
        exit(0)

    # ================================================================
    # Step 4: Identify core genes per group
    # ================================================================
    print("\\nIdentifying core genes per group...")

    # Pre-compute member indices
    sample_col_map = {name: i for i, name in enumerate(all_sample_names)}
    group_member_indices = {}
    for group_id, info in groups.items():
        group_member_indices[group_id] = [sample_col_map[m] for m in info['members']]

    # Track core genes per group: group_id -> list of (gene_id, {sample: seq_id})
    group_core_genes = {gid: [] for gid in groups}

    # Compute min_presence per group
    group_min_presence = {}
    for group_id, info in groups.items():
        n = len(info['members'])
        min_pres = max(2, math.ceil(n * core_threshold))
        if n <= 5:
            min_pres = n
        group_min_presence[group_id] = min_pres
        print(f"  {group_id}: need >= {min_pres}/{n} samples for core")

    # Stream through CSV
    print("\\nScanning gene_presence_absence.csv...")
    row_count = 0
    with open(gpa_path) as f:
        reader = csv.reader(f)
        next(reader)  # skip header

        for row in reader:
            if len(row) <= sample_start:
                continue
            row_count += 1
            gene_id = row[0].strip('"')
            sample_values = row[sample_start:]

            for group_id, indices in group_member_indices.items():
                # Count presence and collect sequence IDs
                present_samples = {}
                for idx in indices:
                    if idx < len(sample_values):
                        val = sample_values[idx].strip('"').strip()
                        if val and val != '-' and val != '':
                            present_samples[all_sample_names[idx]] = val

                if len(present_samples) >= group_min_presence[group_id]:
                    group_core_genes[group_id].append((gene_id, present_samples))

            if row_count % 50000 == 0:
                print(f"  Processed {row_count} genes...")

    print(f"Processed {row_count} total genes")

    for group_id in groups:
        print(f"  {group_id}: {len(group_core_genes[group_id])} core genes")

    # ================================================================
    # Step 5: Load sequences and write per-group outputs
    # ================================================================
    print("\\nLoading sequences from PANTA samples directory...")

    samples_dir = Path("${panta_dir}") / "samples"

    # Determine which samples we actually need sequences for
    needed_samples = set()
    for group_id, core_genes in group_core_genes.items():
        if not core_genes:
            continue
        for _, present_samples in core_genes:
            needed_samples.update(present_samples.keys())

    print(f"Need sequences for {len(needed_samples)} samples")

    # Load sequences for needed samples
    all_sequences = {}
    loaded_count = 0
    for sample_name in sorted(needed_samples):
        sample_fasta = samples_dir / sample_name / f"{sample_name}.fna"
        if not sample_fasta.exists():
            continue

        with open(sample_fasta) as f:
            current_id = None
            current_seq = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        all_sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id:
                all_sequences[current_id] = ''.join(current_seq)

        loaded_count += 1
        if loaded_count % 100 == 0:
            print(f"  Loaded {loaded_count} samples...")

    print(f"Loaded sequences from {loaded_count} samples ({len(all_sequences)} total sequences)")

    # ================================================================
    # Step 6: Write per-group core gene files and genus maps
    # ================================================================
    print("\\nWriting per-group outputs...")

    manifest_rows = []

    for group_id, info in sorted(groups.items()):
        core_genes = group_core_genes[group_id]

        if not core_genes:
            print(f"  {group_id}: no core genes, skipping")
            continue

        group_dir = Path(f"groups/{group_id}")
        core_dir = group_dir / "core_genes"
        core_dir.mkdir(parents=True, exist_ok=True)

        # Write core gene multi-FASTA files
        extracted = 0
        for gene_id, present_samples in core_genes:
            safe_name = gene_id.replace('/', '_').replace('\\\\', '_').replace(':', '_').replace(' ', '_')
            outpath = core_dir / f"{safe_name}.fasta"

            seqs_written = 0
            with open(outpath, 'w') as f:
                for sample_name, seq_id in sorted(present_samples.items()):
                    if seq_id in all_sequences:
                        f.write(f">{sample_name}\\n")
                        seq = all_sequences[seq_id]
                        for i in range(0, len(seq), 80):
                            f.write(seq[i:i+80] + "\\n")
                        seqs_written += 1

            if seqs_written >= 2:
                extracted += 1
            else:
                os.remove(outpath)

        print(f"  {group_id}: {extracted} core gene files written")

        if extracted == 0:
            # Remove empty directory
            import shutil
            shutil.rmtree(group_dir)
            continue

        # Write genus map for tree visualization
        with open(group_dir / "genus_map.tsv", 'w') as f:
            for sample in info['members']:
                genus = sample_genus.get(sample, 'Unknown')
                f.write(f"{sample}\\t{genus}\\n")

        manifest_rows.append({
            'group_id': group_id,
            'sample_count': len(info['members']),
            'n_core_genes': extracted
        })

    # ================================================================
    # Step 7: Write manifest
    # ================================================================
    with open("manifest.tsv", 'w') as f:
        f.write("group_id\\tsample_count\\tn_core_genes\\n")
        for row in manifest_rows:
            f.write(f"{row['group_id']}\\t{row['sample_count']}\\t{row['n_core_genes']}\\n")

    print(f"\\nCompleted! {len(manifest_rows)} groups with core genes written.")
    """
}
