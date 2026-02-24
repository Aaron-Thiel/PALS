/*
 * Select Complete Reference Genome
 *
 * For each sample, joins skani ANI results with genome metadata to find
 * the nearest complete genome (assembly_level = "Complete Genome" or "Chromosome")
 * at species level. Emits the reference path only when a complete genome exists.
 */

process SELECT_COMPLETE_REFERENCE {
    tag "$sample_id"
    label 'process_low'

    container null

    input:
    tuple val(sample_id), path(ani_results), path(genome_metadata), val(genome_dir)

    output:
    tuple val(sample_id), path("selected_complete_ref.tsv"), emit: selection
    tuple val(sample_id), env(REF_FNA), env(REF_GFF), emit: ref_paths, optional: true

    script:
    """
    python3 << 'PYEOF'
import csv
import os

ani_file = "${ani_results}"
metadata_file = "${genome_metadata}"
genome_dir = "${genome_dir}"
sample_id = "${sample_id}"

# ----------------------------------------------------------------
# Step 1: Load genome metadata (accession -> assembly_level + paths)
# ----------------------------------------------------------------
metadata = {}
with open(metadata_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        acc = row['accession']
        metadata[acc] = {
            'assembly_level': row.get('assembly_level', ''),
            'completeness': float(row.get('completeness', 0) or 0),
            'contamination': float(row.get('contamination', 100) or 100),
            'file_path': row.get('file_path', ''),
            'gff_path': row.get('gff_path', ''),
            'species': row.get('species', ''),
        }

print(f"Loaded {len(metadata)} entries from genome_metadata.csv")

# ----------------------------------------------------------------
# Step 2: Parse skani ANI results and join with metadata
# ----------------------------------------------------------------
candidates = []
with open(ani_file) as f:
    reader = csv.DictReader(f, delimiter='\\t')
    for row in reader:
        ref_file = row.get('Ref_file', '')
        # Extract accession from path: .../GCA_000392485.2.fna
        basename = os.path.basename(ref_file)
        accession = basename.replace('.fna', '').replace('.fasta', '')

        if accession not in metadata:
            continue

        meta = metadata[accession]
        assembly_level = meta['assembly_level']

        if assembly_level not in ('Complete Genome', 'Chromosome'):
            continue

        ani = float(row.get('ANI', 0))
        af_query = float(row.get('Align_fraction_query', 0))

        # Resolve absolute paths from relative file_path in metadata
        # metadata file_path is like: ./Species_name/GCA_xxx.fna
        # Use genome_dir (original database path) instead of staged metadata location
        fna_path = os.path.normpath(os.path.join(genome_dir, meta['file_path']))
        gff_path = os.path.normpath(os.path.join(genome_dir, meta['gff_path']))

        candidates.append({
            'accession': accession,
            'species': meta['species'],
            'ani': ani,
            'af_query': af_query,
            'assembly_level': assembly_level,
            'completeness': meta['completeness'],
            'contamination': meta['contamination'],
            'fna_path': fna_path,
            'gff_path': gff_path,
        })

print(f"Complete genome candidates: {len(candidates)}")

# ----------------------------------------------------------------
# Step 3: Select best complete genome (highest ANI)
# ----------------------------------------------------------------
if candidates:
    candidates.sort(key=lambda x: (-x['ani'], -x['af_query']))
    best = candidates[0]

    print(f"Selected: {best['accession']} ({best['species']})")
    print(f"  ANI: {best['ani']}%, AF: {best['af_query']}%")
    print(f"  Assembly level: {best['assembly_level']}")
    print(f"  FNA: {best['fna_path']}")
    print(f"  GFF: {best['gff_path']}")

    with open("selected_complete_ref.tsv", 'w') as f:
        f.write("sample_id\\taccession\\tspecies\\tani\\taf_query\\tassembly_level\\tcompleteness\\tcontamination\\tfna_path\\tgff_path\\n")
        f.write(f"{sample_id}\\t{best['accession']}\\t{best['species']}\\t{best['ani']}\\t{best['af_query']}\\t"
                f"{best['assembly_level']}\\t{best['completeness']}\\t{best['contamination']}\\t"
                f"{best['fna_path']}\\t{best['gff_path']}\\n")

    # Write env file for shell to source
    with open("ref_paths.env", 'w') as f:
        f.write(f"export REF_FNA={best['fna_path']}\\n")
        f.write(f"export REF_GFF={best['gff_path']}\\n")

else:
    print(f"No complete genome found for {sample_id} - skipping visualization")

    with open("selected_complete_ref.tsv", 'w') as f:
        f.write("sample_id\\taccession\\tspecies\\tani\\taf_query\\tassembly_level\\tcompleteness\\tcontamination\\tfna_path\\tgff_path\\n")
        f.write(f"{sample_id}\\tNONE\\t\\t\\t\\t\\t\\t\\t\\t\\n")
PYEOF

    # Source env vars for Nextflow env() capture (only if complete genome found)
    if [ -f ref_paths.env ]; then
        source ref_paths.env
    fi
    """
}
