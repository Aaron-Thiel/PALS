/*
 * Build Subtree Groups
 *
 * Reads the main phylogenetic tree + taxonomy data to identify genus/species
 * groups that qualify for subtree visualization. For each qualifying group,
 * outputs a (group_id, comma-separated tip list) pair.
 *
 * A group qualifies when:
 *   - It contains at least one internal (PAL*) sample
 *   - It has >= subtree_min_tips total tips present in the main tree
 */

process BUILD_SUBTREE_GROUPS {
    tag "build_groups"
    label 'process_low'

    container 'python:3.9'

    input:
    path(main_tree)
    path(taxonomy_csv)
    path(pipeline_summary)

    output:
    path("group_tips_*.csv"), emit: tip_files, optional: true

    script:
    def min_tips = params.subtree_min_tips ?: 3
    """
#!/usr/bin/env python3
import csv
import re
import os

print("=" * 60)
print("Build Subtree Groups")
print("=" * 60)

min_tips = int(${min_tips})
print(f"Min tips per group: {min_tips}")

# ================================================================
# Read tree tip labels from newick
# ================================================================
with open("${main_tree}") as f:
    newick = f.read().strip()
tips_in_tree = set(re.findall(r'([\\w.]+)(?=:)', newick))
print(f"Tips in main tree: {len(tips_in_tree)}")

# ================================================================
# Build taxonomy maps
# ================================================================
sample_genus = {}
sample_species = {}
internal_genera = set()
internal_species = set()

# Internal taxonomy from pipeline_summary
if os.path.exists("${pipeline_summary}"):
    with open("${pipeline_summary}") as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            sid = row['sample_id']
            genus = row.get('genus', '')
            species = row.get('species', '')
            if genus and genus != 'NA':
                sample_genus[sid] = genus
                internal_genera.add(genus)
            if species and species != 'NA':
                sp = species.replace(' ', '_')
                sample_species[sid] = sp
                internal_species.add(sp)

# External taxonomy from CSV
with open("${taxonomy_csv}") as f:
    reader = csv.DictReader(f)
    for row in reader:
        gid = row['genome_id']
        genus = row.get('genus', '')
        species = row.get('species', '')
        if genus:
            sample_genus[gid] = genus
        if species:
            sample_species[gid] = species.replace(' ', '_')

print(f"Internal genera: {sorted(internal_genera)}")
print(f"Internal species: {sorted(internal_species)}")
print(f"Total genus mappings: {len(sample_genus)}")
print(f"Total species mappings: {len(sample_species)}")

# ================================================================
# Build groups (only for genera/species with internal samples)
# ================================================================
groups = {}

for genus in sorted(internal_genera):
    group_id = f"genus_{genus}"
    members = sorted([t for t in tips_in_tree if sample_genus.get(t) == genus])
    if len(members) >= min_tips:
        groups[group_id] = members
        print(f"  {group_id}: {len(members)} tips")
    else:
        print(f"  {group_id}: SKIPPED ({len(members)} < {min_tips})")

for species in sorted(internal_species):
    group_id = f"species_{species}"
    members = sorted([t for t in tips_in_tree if sample_species.get(t) == species])
    if len(members) >= min_tips:
        groups[group_id] = members
        print(f"  {group_id}: {len(members)} tips")
    else:
        print(f"  {group_id}: SKIPPED ({len(members)} < {min_tips})")

print(f"\\nTotal qualifying groups: {len(groups)}")

# ================================================================
# Write one file per group (group_id<TAB>comma_separated_tips)
# ================================================================
for group_id, members in groups.items():
    with open(f"group_tips_{group_id}.csv", 'w') as f:
        f.write(f"{group_id}\\t{','.join(members)}\\n")

print("Done.")
"""
}
