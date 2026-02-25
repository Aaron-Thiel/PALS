/*
 * Build Codon Usage References
 *
 * Two-step process:
 * 1. EXTRACT_CDS - Extract coding sequences using gffread
 * 2. BUILD_*_REFERENCE - Calculate codon usage statistics
 */

/*
 * Extract CDS from FNA + GFF3 files using gffread
 * Input: tuple(name, [fna_files], [gff_files])
 * Output: tuple(name, merged_cds.ffn)
 */
process EXTRACT_CDS {
    tag "$name"
    label 'process_low'

    container 'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4'

    input:
    tuple val(name), path(fna_files), path(gff_files)

    output:
    tuple val(name), path("${name}.merged_cds.ffn"), emit: cds

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Create output file
    touch "${name}.merged_cds.ffn"

    # Process each FNA/GFF3 pair
    fna_array=(${fna_files})
    gff_array=(${gff_files})

    for i in "\${!fna_array[@]}"; do
        fna="\${fna_array[\$i]}"
        gff="\${gff_array[\$i]}"
        accession=\$(basename "\$fna" .fna)

        if [[ -f "\$fna" && -f "\$gff" ]]; then
            # Extract CDS
            if gffread -x "temp_\${accession}.ffn" -g "\$fna" "\$gff" 2>/dev/null; then
                if [[ -s "temp_\${accession}.ffn" ]]; then
                    # Prefix headers with accession and append
                    sed "s/^>/>\${accession}|/" "temp_\${accession}.ffn" >> "${name}.merged_cds.ffn"
                fi
                rm -f "temp_\${accession}.ffn"
            fi
        fi
    done

    # Ensure output exists
    if [[ ! -s "${name}.merged_cds.ffn" ]]; then
        echo "# No CDS extracted for ${name}" > "${name}.merged_cds.ffn"
    fi

    echo "Extracted CDS for ${name}"
    """
}

/*
 * Build species-level codon reference from extracted CDS
 * Input: tuple(species, merged_cds.ffn)
 * Output: tuple(species, codon_reference.json)
 *
 * Calculates per-genome RSCU and frequencies to enable confidence intervals.
 */
process BUILD_SPECIES_REFERENCE {
    tag "$species"
    label 'process_low'

    // Store reference alongside source genomes
    storeDir "${params.reference_genomes_dir}/${species.split('_')[0]}/${species}"

    container 'quay.io/biocontainers/pandas:2.2.1'
    conda 'conda-forge::pandas conda-forge::numpy'

    input:
    tuple val(species), path(cds_file)

    output:
    tuple val(species), path("codon_reference.json"), emit: reference

    script:
    def genus = species.split('_')[0]
    """
    #!/usr/bin/env python3
    import json
    import numpy as np
    from collections import defaultdict

    CODON_TABLE = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    AA_TO_CODONS = defaultdict(list)
    for codon, aa in CODON_TABLE.items():
        AA_TO_CODONS[aa].append(codon)

    SYNONYMOUS_CODONS = [c for c in CODON_TABLE if CODON_TABLE[c] not in ['M', 'W', '*']]

    def calculate_rscu(codon_counts):
        rscu = {}
        for aa, codons in AA_TO_CODONS.items():
            if aa in ['M', 'W', '*']:
                continue
            total = sum(codon_counts.get(c, 0) for c in codons)
            expected = total / len(codons) if len(codons) > 0 else 0
            for codon in codons:
                observed = codon_counts.get(codon, 0)
                rscu[codon] = observed / expected if expected > 0 else 0.0
        return rscu

    def calculate_freqs(codon_counts):
        freqs = {}
        for aa, codons in AA_TO_CODONS.items():
            total = sum(codon_counts.get(c, 0) for c in codons)
            for codon in codons:
                count = codon_counts.get(codon, 0)
                freqs[codon] = count / total if total > 0 else 0.0
        return freqs

    # Count codons per genome for std calculation
    genome_codon_counts = defaultdict(lambda: defaultdict(int))
    total_codon_counts = defaultdict(int)
    total_genes = 0
    total_codons = 0

    current_seq = []
    current_genome = None

    with open("${cds_file}") as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            if line.startswith('>'):
                if current_seq and current_genome:
                    seq = ''.join(current_seq).upper().replace('U', 'T')
                    for i in range(0, len(seq) - 2, 3):
                        codon = seq[i:i+3]
                        if all(b in 'TCAG' for b in codon):
                            genome_codon_counts[current_genome][codon] += 1
                            total_codon_counts[codon] += 1
                            total_codons += 1
                    total_genes += 1
                header = line[1:]
                current_genome = header.split('|')[0] if '|' in header else 'unknown'
                current_seq = []
            else:
                current_seq.append(line)

        # Process last sequence
        if current_seq and current_genome:
            seq = ''.join(current_seq).upper().replace('U', 'T')
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i:i+3]
                if all(b in 'TCAG' for b in codon):
                    genome_codon_counts[current_genome][codon] += 1
                    total_codon_counts[codon] += 1
                    total_codons += 1
            total_genes += 1

    genomes = list(genome_codon_counts.keys())
    n_genomes = len(genomes)

    # Calculate pooled RSCU and frequencies (for mean)
    pooled_rscu = calculate_rscu(dict(total_codon_counts))
    pooled_freqs = calculate_freqs(dict(total_codon_counts))

    # Calculate per-genome RSCU and frequencies, then compute std
    rscu_std = {}
    freq_std = {}

    if n_genomes >= 2:
        # Calculate RSCU and freq for each genome
        genome_rscus = [calculate_rscu(dict(genome_codon_counts[g])) for g in genomes]
        genome_freqs = [calculate_freqs(dict(genome_codon_counts[g])) for g in genomes]

        # Calculate std for each codon
        for codon in SYNONYMOUS_CODONS:
            rscu_values = [gr.get(codon, 0) for gr in genome_rscus]
            freq_values = [gf.get(codon, 0) for gf in genome_freqs]
            rscu_std[codon] = float(np.std(rscu_values)) if rscu_values else 0.0
            freq_std[codon] = float(np.std(freq_values)) if freq_values else 0.0

    reference = {
        'level': 'species',
        'name': '${species}',
        'genus': '${genus}',
        'n_genomes': n_genomes,
        'n_genes': total_genes,
        'n_codons': total_codons,
        'codon_counts': dict(total_codon_counts),
        'rscu': pooled_rscu,
        'rscu_std': rscu_std,
        'freq_std': freq_std
    }

    with open('codon_reference.json', 'w') as f:
        json.dump(reference, f, indent=2)

    print(f"Built species reference: {n_genomes} genomes, {total_genes} genes, {total_codons} codons")
    """
}

/*
 * Build genus-level codon reference from extracted CDS
 * Input: tuple(genus, merged_cds.ffn)
 * Output: tuple(genus, codon_reference.json)
 *
 * Calculates per-genome RSCU and frequencies to enable confidence intervals.
 */
process BUILD_GENUS_REFERENCE {
    tag "$genus"
    label 'process_low'

    // Store reference at genus level
    storeDir "${params.reference_genomes_dir}/${genus}"

    container 'quay.io/biocontainers/pandas:2.2.1'
    conda 'conda-forge::pandas conda-forge::numpy'

    input:
    tuple val(genus), path(cds_file)

    output:
    tuple val(genus), path("codon_reference.json"), emit: reference

    script:
    """
    #!/usr/bin/env python3
    import json
    import numpy as np
    from collections import defaultdict

    CODON_TABLE = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    AA_TO_CODONS = defaultdict(list)
    for codon, aa in CODON_TABLE.items():
        AA_TO_CODONS[aa].append(codon)

    SYNONYMOUS_CODONS = [c for c in CODON_TABLE if CODON_TABLE[c] not in ['M', 'W', '*']]

    def calculate_rscu(codon_counts):
        rscu = {}
        for aa, codons in AA_TO_CODONS.items():
            if aa in ['M', 'W', '*']:
                continue
            total = sum(codon_counts.get(c, 0) for c in codons)
            expected = total / len(codons) if len(codons) > 0 else 0
            for codon in codons:
                observed = codon_counts.get(codon, 0)
                rscu[codon] = observed / expected if expected > 0 else 0.0
        return rscu

    def calculate_freqs(codon_counts):
        freqs = {}
        for aa, codons in AA_TO_CODONS.items():
            total = sum(codon_counts.get(c, 0) for c in codons)
            for codon in codons:
                count = codon_counts.get(codon, 0)
                freqs[codon] = count / total if total > 0 else 0.0
        return freqs

    # Count codons per genome for std calculation
    genome_codon_counts = defaultdict(lambda: defaultdict(int))
    total_codon_counts = defaultdict(int)
    total_genes = 0
    total_codons = 0

    current_seq = []
    current_genome = None

    with open("${cds_file}") as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            if line.startswith('>'):
                if current_seq and current_genome:
                    seq = ''.join(current_seq).upper().replace('U', 'T')
                    for i in range(0, len(seq) - 2, 3):
                        codon = seq[i:i+3]
                        if all(b in 'TCAG' for b in codon):
                            genome_codon_counts[current_genome][codon] += 1
                            total_codon_counts[codon] += 1
                            total_codons += 1
                    total_genes += 1
                header = line[1:]
                current_genome = header.split('|')[0] if '|' in header else 'unknown'
                current_seq = []
            else:
                current_seq.append(line)

        # Process last sequence
        if current_seq and current_genome:
            seq = ''.join(current_seq).upper().replace('U', 'T')
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i:i+3]
                if all(b in 'TCAG' for b in codon):
                    genome_codon_counts[current_genome][codon] += 1
                    total_codon_counts[codon] += 1
                    total_codons += 1
            total_genes += 1

    genomes = list(genome_codon_counts.keys())
    n_genomes = len(genomes)

    # Calculate pooled RSCU and frequencies (for mean)
    pooled_rscu = calculate_rscu(dict(total_codon_counts))
    pooled_freqs = calculate_freqs(dict(total_codon_counts))

    # Calculate per-genome RSCU and frequencies, then compute std
    rscu_std = {}
    freq_std = {}

    if n_genomes >= 2:
        # Calculate RSCU and freq for each genome
        genome_rscus = [calculate_rscu(dict(genome_codon_counts[g])) for g in genomes]
        genome_freqs = [calculate_freqs(dict(genome_codon_counts[g])) for g in genomes]

        # Calculate std for each codon
        for codon in SYNONYMOUS_CODONS:
            rscu_values = [gr.get(codon, 0) for gr in genome_rscus]
            freq_values = [gf.get(codon, 0) for gf in genome_freqs]
            rscu_std[codon] = float(np.std(rscu_values)) if rscu_values else 0.0
            freq_std[codon] = float(np.std(freq_values)) if freq_values else 0.0

    reference = {
        'level': 'genus',
        'name': '${genus}',
        'n_genomes': n_genomes,
        'n_genes': total_genes,
        'n_codons': total_codons,
        'codon_counts': dict(total_codon_counts),
        'rscu': pooled_rscu,
        'rscu_std': rscu_std,
        'freq_std': freq_std
    }

    with open('codon_reference.json', 'w') as f:
        json.dump(reference, f, indent=2)

    print(f"Built genus reference: {n_genomes} genomes, {total_genes} genes, {total_codons} codons")
    """
}
