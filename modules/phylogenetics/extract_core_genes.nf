/*
 * Extract Core Genes from PANTA Output
 * Identifies genes present in >=X% of genomes
 */

process EXTRACT_CORE_GENES {
    tag "$sample_id"
    label 'process_low'
    
    publishDir "${params.outdir}/phylogenetics/core_genes", mode: 'copy'
    
    container 'python:3.9'
    
    input:
    tuple val(sample_id), path(panta_dir), val(sample_count)
    
    output:
    tuple val(sample_id), path("core_genes"), emit: core_gene_seqs
    tuple val(sample_id), path("core_genes.txt"), emit: core_gene_list
    tuple val(sample_id), val(sample_count), emit: sample_count
    
    script:
    def core_threshold = params.phylo_core_threshold ?: 0.95
    """
    #!/usr/bin/env python3
    
    import os
    import csv
    import json
    from pathlib import Path
    
    print("==============================================")
    print("Core Gene Extraction")
    print("==============================================")
    print(f"Sample ID: ${sample_id}")
    print(f"Sample count: ${sample_count}")
    print(f"Core threshold: ${core_threshold}")
    print("")
    
    panta_dir = "${panta_dir}"
    sample_count = int(${sample_count})
    core_threshold = float(${core_threshold})

    # For small sample sizes, be stricter with threshold
    # With 3 samples: 95% -> 3 samples (require all)
    # With 10 samples: 95% -> 10 samples
    # With 100 samples: 95% -> 95 samples
    import math
    min_presence = max(2, math.ceil(sample_count * core_threshold))

    # For very small samples (<=5), require presence in ALL samples for true core genes
    if sample_count <= 5:
        min_presence = sample_count
        print(f"Small sample size ({sample_count}) - requiring genes present in ALL samples")

    print(f"Core genes must be present in >= {min_presence} of {sample_count} samples ({core_threshold:.0%})")
    
    # Find gene presence/absence file
    gpa_path = None
    for name in ['gene_presence_absence.csv', 'gene_presence_absence.Rtab', 'clustered_proteins.csv']:
        p = Path(panta_dir) / name
        if p.exists():
            gpa_path = p
            break
    
    if not gpa_path:
        print(f"ERROR: No gene presence/absence file found in {panta_dir}")
        print("Available files:")
        for f in Path(panta_dir).iterdir():
            print(f"  - {f.name}")
        exit(1)
    
    print(f"Reading: {gpa_path.name}")
    
    # Parse and identify core genes
    core_genes = []
    total_genes = 0
    sample_names = []

    with open(gpa_path) as f:
        reader = csv.reader(f)
        header = next(reader)

        # Clean header values (remove quotes)
        header = [h.strip('"').strip() for h in header]

        # Panta CSV format: Gene, Annotation, No. isolates, No. sequences, ..., Sample1, Sample2, ...
        # Find the "No. isolates" column which has the presence count
        isolates_col = None
        for i, col in enumerate(header):
            if col == 'No. isolates':
                isolates_col = i
                break

        # Identify sample columns (after metadata columns)
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

        sample_names = header[sample_start:]
        print(f"Found {len(sample_names)} samples: {', '.join(sample_names)}")
        print(f"Using 'No. isolates' column (index {isolates_col}) for presence count")

        for row in reader:
            if len(row) <= sample_start:
                continue

            total_genes += 1
            gene_id = row[0].strip('"')
            annotation = row[1].strip('"') if len(row) > 1 else ''

            # Use "No. isolates" column if available, otherwise count sample columns
            if isolates_col is not None and isolates_col < len(row):
                try:
                    presence = int(row[isolates_col].strip('"'))
                except ValueError:
                    # Fallback to counting sample columns
                    presence = sum(1 for v in row[sample_start:] if v and v.strip() and v.strip('"') and v.strip('"') != '-')
            else:
                # Count non-empty sample columns
                presence = sum(1 for v in row[sample_start:] if v and v.strip() and v.strip('"') and v.strip('"') != '-')

            if presence >= min_presence:
                core_genes.append({
                    'gene_id': gene_id,
                    'presence': presence,
                    'annotation': annotation
                })
    
    print(f"Identified {len(core_genes)} core genes from {total_genes} total clusters")
    
    # Re-read CSV to get per-sample sequence IDs for core genes
    core_gene_samples = {}  # gene_id -> {sample_name: seq_id}

    with open(gpa_path) as f:
        reader = csv.reader(f)
        header = next(reader)
        header = [h.strip('"').strip() for h in header]

        for row in reader:
            if len(row) <= sample_start:
                continue
            gene_id = row[0].strip('"')
            if gene_id not in [g['gene_id'] for g in core_genes]:
                continue

            # Get sequence IDs for each sample
            sample_seqs = {}
            for i, sample_name in enumerate(sample_names):
                col_idx = sample_start + i
                if col_idx < len(row):
                    seq_id = row[col_idx].strip('"').strip()
                    if seq_id and seq_id != '-':
                        sample_seqs[sample_name] = seq_id

            if sample_seqs:
                core_gene_samples[gene_id] = sample_seqs

    print(f"Found sequence IDs for {len(core_gene_samples)} core genes")

    # Write core gene list
    with open("core_genes.txt", 'w') as f:
        f.write(f"# Core genes (present in >= {min_presence} of {sample_count} samples)\\n")
        f.write(f"# Threshold: {core_threshold:.0%}\\n")
        f.write(f"# Total core genes: {len(core_genes)}\\n")
        f.write(f"# Total gene clusters: {total_genes}\\n")
        f.write("#\\n")
        f.write("gene_id\\tpresence\\tannotation\\n")
        for g in sorted(core_genes, key=lambda x: -x['presence']):
            f.write(f"{g['gene_id']}\\t{g['presence']}\\t{g['annotation']}\\n")

    # Extract sequences for core genes from per-sample fasta files
    os.makedirs("core_genes", exist_ok=True)

    if not core_gene_samples:
        print("WARNING: No core genes identified")
        with open("core_genes/empty.fasta", 'w') as f:
            f.write(">no_core_genes\\nNNNN\\n")
        exit(0)

    # Load all per-sample sequences into memory
    # Panta stores per-sample sequences in samples/{sample_name}/{sample_name}.fna
    print("Loading per-sample sequences...")
    all_sequences = {}  # seq_id -> sequence

    samples_dir = Path(panta_dir) / "samples"
    if not samples_dir.exists():
        print(f"ERROR: samples directory not found at {samples_dir}")
        exit(1)

    for sample_name in sample_names:
        sample_fasta = samples_dir / sample_name / f"{sample_name}.fna"
        if not sample_fasta.exists():
            print(f"  WARNING: No fasta for {sample_name}")
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

        print(f"  Loaded {sample_name}: {len([k for k in all_sequences if k.startswith(sample_name)])} sequences")

    print(f"Total sequences loaded: {len(all_sequences)}")

    # Extract core gene sequences - one multi-fasta per gene with all samples
    extracted = 0
    for gene_id, sample_seqs in core_gene_samples.items():
        safe_name = gene_id.replace('/', '_').replace('\\\\', '_').replace(':', '_').replace(' ', '_')
        outpath = f"core_genes/{safe_name}.fasta"

        seqs_written = 0
        with open(outpath, 'w') as f:
            for sample_name, seq_id in sample_seqs.items():
                if seq_id in all_sequences:
                    # Use sample name as sequence header for alignment
                    f.write(f">{sample_name}\\n")
                    seq = all_sequences[seq_id]
                    for i in range(0, len(seq), 80):
                        f.write(seq[i:i+80] + "\\n")
                    seqs_written += 1

        if seqs_written >= 2:  # Only keep genes with at least 2 samples (for alignment)
            extracted += 1
        else:
            os.remove(outpath)  # Remove single-sequence files

    print(f"Extracted {extracted} core gene files (each with multiple samples)")
    print("Core gene extraction completed")
    """
}