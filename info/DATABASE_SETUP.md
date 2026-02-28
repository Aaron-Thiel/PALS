# Database Setup

The PALS pipeline requires several external databases. Some are used by both
pipelines, others only by one. antiSMASH databases are bundled inside its
Docker image and need no separate download.

## Quick Start

```bash
chmod +x info/download_databases.sh
./info/download_databases.sh              # default: databases/
./info/download_databases.sh /my/path     # custom location
```

The script prompts **Yes / No / All** before each download.

## Database Overview

| # | Database | Used by | Pipeline | Approx. Size | Path |
|---|----------|---------|----------|-------------|------|
| 1 | **Bakta** | `BAKTA` (genome annotation) | Both | 30 GB (full) / 1.4 GB (light) | `<db>/bakta/db/` |
| 2 | **CheckM2** | `CHECKM2` (completeness/contamination) | Assembly | 3.5 GB | `<db>/checkm2/` |
| 3 | **BUSCO** | `BUSCO` (completeness assessment) | Assembly | ~300 MB | `<db>/busco/lineages/` |
| 4 | **Kraken2** | `KRAKEN2` (read classification) | Assembly | 8 GB (standard-8) | `<db>/kraken2/` |
| 5 | **Sourmash GTDB** | `SOURMASH` (genus classification) | Assembly | ~2 GB | `<db>/sourmash/` |
| 6 | **EggNOG** | `EGGNOG` (functional annotation) | Analysis | 45 GB | `<db>/eggnog/` |
| 7 | **Pfam-A 35.0** | `BIGSCAPE` (BGC similarity networks) | Analysis | 300 MB | `<db>/antismash/pfam/35.0/Pfam-A.hmm` |
| 8 | **antiSMASH** | `ANTISMASH` (BGC detection) | Analysis | built into Docker image | — |

## Database Details

### 1. Bakta

Bacterial genome annotation database. Choose **full** for comprehensive annotation
or **light** for faster downloads with slightly reduced sensitivity.

- Pipeline: Assembly + Analysis
- Docker image: `staphb/bakta:latest`
- Mounted at `/databases/bakta/db` inside the container
- Download tool: `bakta_db download` or direct from [Zenodo](https://zenodo.org/records/10522951)

### 2. CheckM2

Genome quality assessment using machine learning. Evaluates genome completeness
and contamination. Required by both the main QC step and the reference genome
classification workflow.

- Pipeline: Assembly
- Docker image: `staphb/checkm2:latest`
- Mounted at `/databases/checkm2` inside the container
- Key file: `uniref100.KO.1.dmnd` (DIAMOND database)
- Download tool: `checkm2 database --download --path <path>`
- Source: [CheckM2 GitHub](https://github.com/chklovski/CheckM2)

### 3. BUSCO

Benchmarking Universal Single-Copy Orthologs for genome completeness assessment.
The pipeline runs BUSCO in **offline mode**, so the lineage must be pre-downloaded.

- Pipeline: Assembly
- Docker image: `staphb/busco:latest`
- Mounted at `/databases/busco` inside the container
- Download: `busco --download_path <path> --download <lineage>`
- Source: [BUSCO datasets](https://busco-data.ezlab.org/v5/data/lineages/)

**Lineage selection:**

| Lineage | Scope | Markers | Use when |
|---------|-------|---------|----------|
| `lactobacillaceae_odb12` | Lactobacillaceae family | ~402 | Default for this pipeline |
| `lactobacillales_odb10` | Lactobacillales order | ~402 | Broader lactic acid bacteria |
| `bacilli_odb10` | Bacilli class | ~450 | Mixed Firmicutes samples |
| `bacteria_odb10` | All bacteria | ~124 | Any bacterial genome |

The pipeline hardcodes `lactobacillaceae_odb12` in `modules/qc/busco.nf`.
To use a different lineage, update the `-l` flag in that module.

### 4. Kraken2

Taxonomic classification of sequencing reads. Used in the optional preprocessing
step to filter reads by taxonomy before assembly. Only needed when
`run_preprocessing = true`. The pipeline uses Kraken2 output with KrakenTools
to filter for Lactobacillaceae reads and remove contamination.

- Pipeline: Assembly (preprocessing)
- Docker image: `staphb/kraken2:latest`
- Mounted at `/databases/kraken2` inside the container
- Download tool: `kraken2-build --standard --db <path>` or download pre-built from
  [Kraken2 databases](https://benlangmead.github.io/aws-indexes/k2)

**Database selection:**

| Database | Size | Contents | Use when |
|----------|------|----------|----------|
| Standard-8 | ~8 GB | RefSeq archaea, bacteria, viral, human (capped) | Default — good balance of size and sensitivity |
| Standard-16 | ~16 GB | Same as above, higher resolution | More memory available |
| Standard | ~70 GB | Full RefSeq, no capping | Maximum sensitivity, requires lots of RAM |
| PlusPF-8 | ~8 GB | Standard + protozoa + fungi (capped) | Samples with potential eukaryotic contamination |

The download script defaults to **Standard-8**. For Lactobacillaceae-focused
work this is sufficient, as the standard databases include comprehensive
bacterial RefSeq coverage. There is no Lactobacillaceae-specific Kraken2
database — all variants use the same bacterial references.

### 5. Sourmash GTDB

Sourmash sketches of GTDB representative genomes for genus-level taxonomic
classification. Used when `sourmash_classify = true` (default) to automatically
identify the genus and download reference genomes.

- Pipeline: Assembly (classification)
- Docker image: `nanozoo/sourmash:4.8.14--4df0447`
- Mounted at `/databases/sourmash` inside the container
- Required files:
  - `gtdb-reps-rs226-k31.dna.zip` — GTDB RS226 representatives (k=31)
  - `gtdb-rs226.lineages.csv` — taxonomy lineage mapping
- Source: [Sourmash prepared databases](https://sourmash.readthedocs.io/en/latest/databases.html)

### 6. EggNOG

Orthology-based functional annotation. Used by eggnog-mapper in DIAMOND mode
with `Lactobacillaceae` as the target taxon scope.

- Pipeline: Analysis (pangenomics block)
- Docker image: `nanozoo/eggnog-mapper:2.1.13--c16a7d2`
- Download tool: `download_eggnog_data.py --data_dir <path> -y`
- Can also be downloaded from inside the Docker container

### 7. Pfam-A 35.0

Protein family HMM profiles required by BiG-SCAPE 2 (`-p` flag, mandatory).
Only needed if the BGC block is enabled (`bgc_enable = true`).

- Pipeline: Analysis (BGC block)
- Docker image: `ghcr.io/medema-group/big-scape:2.0.0-beta.6`
- Source: [EBI Pfam FTP](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/)
- Running `hmmpress` on it is optional but speeds up the first BiG-SCAPE run

### 8. antiSMASH (likely no download needed)

The `nanozoo/antismash:8.0.0--b6973cb` Docker image is expected to ship with
all required antiSMASH databases pre-installed (similar to the official
`antismash/standalone` image). The pipeline does not pass a `--databases` flag
to antiSMASH, relying on the container's built-in databases.

If antiSMASH fails with database-related errors, you may need to download the
databases separately. See the
[antiSMASH installation docs](https://docs.antismash.secondarymetabolites.org/install/)
for instructions, or consider switching to the official `antismash/standalone`
image which explicitly bundles all databases.

## Expected Directory Layout

All databases under a single root directory (default: `databases/` relative to
the project, or `/BGC-data/databases` for the shared volume mount):

```
databases/
├── bakta/
│   └── db/                                    # Bakta annotation DB
│       ├── version.json
│       └── ...
├── checkm2/                                   # CheckM2 quality DB
│   └── uniref100.KO.1.dmnd
├── busco/                                     # BUSCO lineages
│   └── lineages/
│       └── lactobacillaceae_odb12/
├── kraken2/                                   # Kraken2 taxonomy DB
│   ├── hash.k2d
│   ├── opts.k2d
│   └── taxo.k2d
├── sourmash/                                  # Sourmash GTDB sketches
│   ├── gtdb-reps-rs226-k31.dna.zip
│   └── gtdb-rs226.lineages.csv
├── eggnog/                                    # EggNOG mapper DB
│   ├── eggnog.db
│   ├── eggnog_proteins.dmnd
│   └── ...
└── antismash/
    └── pfam/
        └── 35.0/
            └── Pfam-A.hmm                     # Pfam HMMs for BiG-SCAPE
```

## Docker Volume Mounts

The databases are made available inside Docker containers via two volume mounts
configured in `assembly.config` and `analysis.config`:

```groovy
docker {
    runOptions = "... -v <project>/databases:/databases -v /BGC-data/databases:/BGC-data/databases ..."
}
```

- **`/databases`** — primary mount for most tools (Bakta, CheckM2, BUSCO, Kraken2, Sourmash)
- **`/BGC-data/databases`** — additional mount for EggNOG and Pfam

If you change the database location, update both the Docker `-v` mount and the
relevant `params.*` in the config files.

## Pipeline Config Reference

Database paths are configured in the config files:

**assembly.config:**
```groovy
params {
    database_base   = '/BGC-data/databases'
    genome_database = "${params.database_base}/genomes"
    kraken2_db      = "${params.database_base}/kraken2"
}
```

**analysis.config:**
```groovy
params {
    database_base   = '/BGC-data/databases'
    eggnog_db       = "${params.database_base}/eggnog"
    antismash_db    = "${params.database_base}/antismash"
    pfam_db         = '/BGC-data/databases/antismash/pfam/35.0/Pfam-A.hmm'
}
```

## Disk Space

| Scenario | Approx. Total |
|----------|---------------|
| All databases (Bakta full) | ~90 GB |
| All databases (Bakta light) | ~60 GB |
| Assembly pipeline only | ~45 GB (full Bakta) / ~15 GB (light Bakta) |
| Analysis pipeline only | ~75 GB |

Ensure sufficient free space before running the download script.

## Disclaimer

The download URLs, version numbers, and size estimates in this document and in
`download_databases.sh` are provided as-is with no warranty of correctness or
completeness. External databases are maintained by their respective projects and
may change, move, or become unavailable at any time. Always verify that the
downloaded databases match the versions expected by the tools in use.
Downloading and using these databases is at your own risk. Refer to each
database's official documentation and license terms before use.
