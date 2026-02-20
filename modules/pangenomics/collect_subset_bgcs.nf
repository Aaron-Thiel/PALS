/*
 * Collect antiSMASH BGC files by family/genus/species subset
 *
 * Groups antiSMASH region .gbk files by family, genus and species using taxonomy
 * from pipeline_summary.tsv. Creates per-group directories structured for
 * BiG-SCAPE input (subdirectories per sample with region .gbk files).
 *
 * Only creates groups for family/genera/species with at least min_samples internal samples.
 */

process COLLECT_SUBSET_BGCS {
    tag "collect_bgcs"
    label 'process_low'

    container 'python:3.11'

    input:
    path(antismash_dirs)
    path(pipeline_summary)

    output:
    path("groups/*"), emit: group_dirs
    path("bgc_manifest.tsv"), emit: manifest

    script:
    def min_samples = params.bigscape_min_samples ?: 2
    """
    #!/usr/bin/env python3

    import csv
    import os
    import shutil
    import glob

    print("=" * 60)
    print("Collect Subset BGCs for BiG-SCAPE")
    print("=" * 60)
    print(f"Pipeline summary: ${pipeline_summary}")
    print(f"Min samples per group: ${min_samples}")
    print("")

    min_samples = int(${min_samples})

    # ================================================================
    # Step 1: Read pipeline_summary.tsv for sample -> family/genus/species mapping
    # ================================================================
    sample_taxonomy = {}  # sample_id -> (family, genus, species)

    with open("${pipeline_summary}") as f:
        reader = csv.DictReader(f, delimiter="\\t")
        for row in reader:
            sid = row['sample_id']
            family = row.get('family', '')
            genus = row.get('genus', '')
            species = row.get('species', '')
            if family and family != 'NA' and family != '':
                sample_taxonomy[sid] = (
                    family,
                    genus if genus and genus != 'NA' else '',
                    species if species and species != 'NA' else '',
                )

    print(f"Samples with taxonomy: {len(sample_taxonomy)}")

    # ================================================================
    # Step 2: Find all antiSMASH directories and their region .gbk files
    # ================================================================
    sample_bgcs = {}  # sample_id -> list of region .gbk file paths

    # antiSMASH dirs are staged flat in the work directory
    for entry in os.listdir('.'):
        if not os.path.isdir(entry):
            continue
        # The entry is the sample_id directory from ANTISMASH output
        sample_id = entry
        if sample_id == 'groups':
            continue

        # Find region .gbk files (BiG-SCAPE looks for files with "region" in name)
        region_files = glob.glob(os.path.join(entry, '*.region*.gbk'))
        if region_files:
            sample_bgcs[sample_id] = region_files
            print(f"  {sample_id}: {len(region_files)} BGC regions")

    print(f"\\nSamples with BGC regions: {len(sample_bgcs)}")

    # ================================================================
    # Step 3: Build family/genus/species groups
    # ================================================================
    groups = {}  # group_id -> list of sample_ids

    # Family-level groups
    family_samples = {}
    for sid, (family, genus, species) in sample_taxonomy.items():
        if sid in sample_bgcs:
            family_samples.setdefault(family, []).append(sid)

    for family, samples in sorted(family_samples.items()):
        if len(samples) >= min_samples:
            group_id = f"family_{family}"
            groups[group_id] = samples
            n_bgcs = sum(len(sample_bgcs[s]) for s in samples)
            print(f"  {group_id}: {len(samples)} samples, {n_bgcs} BGCs")
        else:
            print(f"  family_{family}: SKIPPED ({len(samples)} < {min_samples} min)")

    # Genus-level groups
    genus_samples = {}
    for sid, (family, genus, species) in sample_taxonomy.items():
        if sid in sample_bgcs and genus:
            genus_samples.setdefault(genus, []).append(sid)

    for genus, samples in sorted(genus_samples.items()):
        if len(samples) >= min_samples:
            group_id = f"genus_{genus}"
            groups[group_id] = samples
            n_bgcs = sum(len(sample_bgcs[s]) for s in samples)
            print(f"  {group_id}: {len(samples)} samples, {n_bgcs} BGCs")
        else:
            print(f"  genus_{genus}: SKIPPED ({len(samples)} < {min_samples} min)")

    # Species-level groups
    species_samples = {}
    for sid, (family, genus, species) in sample_taxonomy.items():
        if sid in sample_bgcs and species:
            species_samples.setdefault(species, []).append(sid)

    for species, samples in sorted(species_samples.items()):
        if len(samples) >= min_samples:
            group_id = f"species_{species}"
            groups[group_id] = samples
            n_bgcs = sum(len(sample_bgcs[s]) for s in samples)
            print(f"  {group_id}: {len(samples)} samples, {n_bgcs} BGCs")
        else:
            print(f"  species_{species}: SKIPPED ({len(samples)} < {min_samples} min)")

    print(f"\\nTotal groups: {len(groups)}")

    if not groups:
        print("WARNING: No qualifying groups. Creating empty manifest.")
        os.makedirs("groups/empty", exist_ok=True)
        with open("bgc_manifest.tsv", 'w') as f:
            f.write("group_id\\tlevel\\tname\\tn_samples\\tn_bgcs\\n")
        exit(0)

    # ================================================================
    # Step 4: Create per-group directories with region .gbk files
    # ================================================================
    print("\\nCreating group directories...")

    for group_id, samples in sorted(groups.items()):
        group_dir = os.path.join("groups", group_id)
        total_bgcs = 0

        for sample_id in samples:
            sample_dir = os.path.join(group_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)

            for gbk_path in sample_bgcs[sample_id]:
                # Prefix GBK filename with sample_id to avoid collisions
                # (e.g. contig_1.region001.gbk → PAL001_contig_1.region001.gbk)
                basename = os.path.basename(gbk_path)
                dst = os.path.join(sample_dir, f"{sample_id}_{basename}")
                shutil.copy2(gbk_path, dst)
                total_bgcs += 1

        print(f"  {group_id}: {len(samples)} samples, {total_bgcs} BGC files")

    # ================================================================
    # Step 5: Write manifest
    # ================================================================
    with open("bgc_manifest.tsv", 'w') as f:
        f.write("group_id\\tlevel\\tname\\tn_samples\\tn_bgcs\\n")
        for group_id, samples in sorted(groups.items()):
            parts = group_id.split("_", 1)
            level = parts[0]
            name = parts[1]
            n_bgcs = sum(len(sample_bgcs[s]) for s in samples)
            f.write(f"{group_id}\\t{level}\\t{name}\\t{len(samples)}\\t{n_bgcs}\\n")

    print("\\nBGC collection completed!")
    """
}
