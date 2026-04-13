# Docker Images

All Docker images used by the PALS pipeline. Most are pulled automatically on first run.

## Assembly Pipeline

| Image | Process | Purpose |
|-------|---------|---------|
| `staphb/fastqc:latest` | FASTQC | Read quality assessment |
| `staphb/fastp:latest` | FASTP | Read trimming and adapter removal |
| `staphb/kraken2:latest` | KRAKEN2 | Taxonomic classification of reads |
| `staphb/krakentools:latest` | KRAKENTOOLS | Read filtering by taxon |
| `staphb/spades:latest` | SPADES | De novo genome assembly |
| `staphb/checkm2:latest` | CHECKM2 | Genome completeness and contamination |
| `staphb/quast:latest` | QUAST | Assembly statistics |
| `staphb/busco:latest` | BUSCO | Genome completeness (BUSCO markers) |
| `staphb/bakta:latest` | BAKTA | Bacterial genome annotation |
| `nanozoo/sourmash:4.8.14--4df0447` | SOURMASH | Genus-level taxonomic classification |
| `staphb/skani:latest` | SKANI | Average nucleotide identity (ANI) |
| `aaronthiel/panta-nextflow:latest` | PANTA | Per-sample pangenome analysis |
| `aaronthiel/pasa:latest` | PASA | Pangenome-assisted scaffolding |

## Analysis Pipeline

| Image | Process | Purpose |
|-------|---------|---------|
| `staphb/bakta:latest` | BAKTA | Genome annotation (scaffold re-annotation) |
| `aaronthiel/panta-nextflow:latest` | PANTA_COHORT | Cohort-level pangenome analysis |
| `nanozoo/eggnog-mapper:2.1.13--c16a7d2` | EGGNOG | Functional annotation (COG, GO, KEGG) |
| `aaronthiel/pangenome-plots:latest` | PANGENOME_PLOTS | Pangenome visualization |
| `staphb/mafft:latest` | ALIGN_CORE_GENES | Multiple sequence alignment |
| `python:3.9` | CONCATENATE_ALIGNMENT | Alignment concatenation |
| `staphb/iqtree3:latest` | IQTREE3 | Maximum likelihood phylogenetic trees |
| `aaronthiel/tree-viz:latest` | TREE_VISUALIZATION | Phylogenetic tree plots |
| `nanozoo/antismash:8.0.0--b6973cb` | ANTISMASH | Biosynthetic gene cluster detection |
| `ghcr.io/medema-group/big-scape:2.0.0-beta.6` | BIGSCAPE | BGC similarity networks |
| `aaronthiel/genomeviz:latest` | GENE_COMPARISON | Scaffold vs. contig gene comparison |

## Conda-Only Processes

Some visualization processes use Conda environments instead of Docker:

| Process | Dependencies |
|---------|-------------|
| NEW_GENES_OVERVIEW | `pandas`, `matplotlib` |

## Pre-pulling Images

To avoid download delays during pipeline execution, you can pre-pull all images:

```bash
docker pull staphb/spades:latest
docker pull staphb/bakta:latest
docker pull staphb/checkm2:latest
docker pull staphb/quast:latest
docker pull staphb/busco:latest
docker pull staphb/fastqc:latest
docker pull staphb/fastp:latest
docker pull staphb/kraken2:latest
docker pull staphb/krakentools:latest
docker pull staphb/mafft:latest
docker pull staphb/iqtree3:latest
docker pull staphb/skani:latest
docker pull nanozoo/sourmash:4.8.14--4df0447
docker pull nanozoo/eggnog-mapper:2.1.13--c16a7d2
docker pull nanozoo/antismash:8.0.0--b6973cb
docker pull ghcr.io/medema-group/big-scape:2.0.0-beta.6
docker pull aaronthiel/panta-nextflow:latest
docker pull aaronthiel/pasa:latest
docker pull aaronthiel/pangenome-plots:latest
docker pull aaronthiel/tree-viz:latest
docker pull aaronthiel/genomeviz:latest
docker pull python:3.9
```
