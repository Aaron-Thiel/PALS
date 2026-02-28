# PALS

**Pangenome-Assisted Lactobacillaceae Scaffolding**

A Nextflow DSL2 pipeline for bacterial genome assembly, pangenome-assisted scaffolding, and downstream comparative genomics analysis. Designed for Lactobacillaceae but adaptable to other bacterial families.

## Overview

PALS consists of two main pipelines that can be run independently or together:

| Pipeline | Entry point | Config | Purpose |
|----------|------------|--------|---------|
| **Assembly** | `assembly.nf` | `assembly.config` | Reads → assembled, annotated, scaffolded genomes |
| **Analysis** | `analysis.nf` | `analysis.config` | Annotations → pangenomics, phylogenetics, BGC analysis |

An orchestrator (`main.nf`) is also provided to run both sequentially.

## Assembly Pipeline

Takes paired-end Illumina reads and produces annotated, scaffolded genomes.

**Steps:**

1. **Preprocessing** *(optional)* — FastQC, fastp trimming, Kraken2 taxonomic filtering
2. **SPAdes** — de novo genome assembly
3. **CheckM2 QC** — completeness and contamination assessment
4. **HQ Filter** — only high-quality assemblies (>90% complete, <5% contamination) proceed
5. **Sourmash Classification** *(optional)* — taxonomic classification, Lactobacillaceae filtering, automated reference genome download via NCBI
6. **Bakta** — bacterial genome annotation
7. **PANTA** — per-sample pangenome analysis (sample vs. references)
8. **PASA** — pangenome-assisted scaffolding using synteny information
9. **PASA Filter** — intelligent scaffold selection based on QC comparison
10. **Summary** — per-sample summary TSV with key metrics

**Input:** A CSV samplesheet with columns `sample_id,read1,read2,reference_dir`.
When `sourmash_classify=true` (default), `reference_dir` can be left empty as references are downloaded automatically.

**Output:** `results/assembly/` with per-sample subdirectories containing scaffolds, annotations, QC reports, and a `pipeline_summary.tsv`.

## Analysis Pipeline

Takes Bakta annotations (from the assembly pipeline or external sources) and runs modular downstream analyses. Each block can be enabled/disabled independently.

**Core (always runs):**
- **BAKTA** — re-annotates PASA scaffolds (or reuses existing SPAdes-stage Bakta with `--spades`)
- **Merge external** — optionally merges pre-annotated external genomes
- **PANTA cohort** — pangenome analysis across all samples (auto-detects cached results)

**Downstream blocks** (toggle via `*_enable` flags in `analysis.config`):

| Block | Flag | Description |
|-------|------|-------------|
| Codon QC | `codon_qc_enable` | Detects assembly anomalies via codon usage patterns |
| Pangenomics | `pangenomics_enable` | EggNOG functional annotation, pangenome plots, genus/species subsets |
| Phylogenetics | `phylogenetics_enable` | Core gene alignment, IQ-TREE phylogeny, subtree pruning |
| BGC | `bgc_enable` | antiSMASH BGC detection, BiG-SCAPE similarity networks |
| New Genes | `new_genes_enable` | Scaffold vs. contig gene comparison, aggregated overview |

When `subset_analysis_enable=true`, each block additionally produces genus- and species-level analyses.

**Input:** Assembly pipeline results directory (default: `results/assembly`), optionally external GFF3/FNA annotations.

**Output:** `results/analysis/` with subdirectories per block (pangenomics, phylogenetics, antismash, genomeviz, etc.).

## Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) >= 21.04.0
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://docs.conda.io/) (used by some visualization processes)

### 1. Clone the repository

```bash
git clone <repository-url>
cd nextflow
```

### 2. Download databases

The pipeline requires several external databases. An interactive download script is provided:

```bash
chmod +x info/download_databases.sh
./info/download_databases.sh              # default: /BGC-data/databases
./info/download_databases.sh /my/path     # custom location
```

This downloads:
- **Bakta** (~30 GB full / ~1.4 GB light) — genome annotation
- **EggNOG** (~45 GB) — functional annotation
- **Pfam-A 35.0** (~300 MB) — BiG-SCAPE BGC networks

antiSMASH databases are bundled in the Docker image and require no separate download.

See [info/DATABASE_SETUP.md](info/DATABASE_SETUP.md) for details on each database and the expected directory layout.

### 3. Prepare input

Create a samplesheet CSV:

```csv
sample_id,read1,read2,reference_dir
SAMPLE1,/path/to/SAMPLE1_R1.fastq.gz,/path/to/SAMPLE1_R2.fastq.gz,
SAMPLE2,/path/to/SAMPLE2_R1.fastq.gz,/path/to/SAMPLE2_R2.fastq.gz,
```

### 4. Run the pipelines

**Assembly:**
```bash
nextflow run assembly.nf -c assembly.config --input input/samplesheet.csv
```

**Analysis** (after assembly completes):
```bash
nextflow run analysis.nf -c analysis.config
```

**Or via the orchestrator:**
```bash
nextflow run main.nf --pipeline assembly --input input/samplesheet.csv
nextflow run main.nf --pipeline analysis
nextflow run main.nf --pipeline all --input input/samplesheet.csv
```

### 5. Resume failed runs

Nextflow caches completed tasks. Resume after a failure with:

```bash
nextflow run assembly.nf -c assembly.config --input input/samplesheet.csv -resume
```

## Configuration

Each pipeline has its own config file with documented parameters:

- **`assembly.config`** — assembly parameters (SPAdes mode, QC thresholds, PASA mode, database paths)
- **`analysis.config`** — analysis parameters (block toggles, PANTA thresholds, phylogenetic model, BiG-SCAPE cutoffs, database paths)
- **`nextflow.config`** — orchestrator settings

Key parameters can also be overridden on the command line:

```bash
nextflow run analysis.nf -c analysis.config \
  --pangenomics_enable true \
  --phylogenetics_enable true \
  --bgc_enable false
```

## Project Structure

```
nextflow/
├── main.nf                  # Pipeline orchestrator
├── assembly.nf              # Assembly pipeline
├── assembly.config          # Assembly configuration
├── analysis.nf              # Analysis pipeline
├── analysis.config          # Analysis configuration
├── nextflow.config          # Shared/orchestrator configuration
├── modules/                 # Individual process definitions
│   ├── bakta.nf
│   ├── spades.nf
│   ├── panta.nf
│   ├── pasa.nf
│   ├── qc.nf
│   ├── hq_filter.nf
│   ├── classification.nf
│   ├── summary.nf
│   ├── pangenomics/         # Pangenomics modules
│   ├── phylogenetics/       # Phylogenetics modules
│   ├── bgc/                 # BGC modules
│   ├── codon/               # Codon bias modules
│   └── new_genes/           # New gene analysis modules
├── workflows/               # Sub-workflow definitions
│   ├── pangenomics.nf
│   ├── phylogenetics.nf
│   ├── bgc.nf
│   ├── codon_bias.nf
│   └── new_genes.nf
├── scripts/                 # Utility scripts
├── input/                   # Samplesheets and input data
├── info/                    # Database setup docs and scripts
│   ├── DATABASE_SETUP.md
│   └── download_databases.sh
└── results/                 # Pipeline output
    ├── assembly/
    └── analysis/
```

## Execution Profiles

The assembly pipeline supports multiple execution profiles:

```bash
# Local execution (default)
nextflow run assembly.nf -c assembly.config --input samplesheet.csv

# SLURM cluster
nextflow run assembly.nf -c assembly.config --input samplesheet.csv -profile slurm

# Singularity instead of Docker
nextflow run assembly.nf -c assembly.config --input samplesheet.csv -profile singularity
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome. See the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines.

## Author

Aaron Thiel
