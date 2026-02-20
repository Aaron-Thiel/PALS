# Lost Genes COG Analysis — Implementation Plan

## Concept
"Lost genes" = protein-coding genes present in SPAdes contigs (pre-PASA) but absent from PASA scaffolds. The inverse of the existing "new genes" analysis.

## Implementation: Standalone Python Script

**File:** `scripts/lost_genes_analysis.py`

A standalone script (no Nextflow re-run needed) that reads existing results and produces output in `results/analysis/lost_genes/`.

### Input Data (all existing, no re-run)

| Data | Path |
|------|------|
| Contig proteins | `results/assembly/PAL*/bakta/*.faa` |
| Scaffold proteins | `results/analysis/bakta/PAL*/*.faa` |
| Taxonomy | `results/assembly/pipeline_summary.tsv` |
| SPAdes PANTA CSV | `results/analysis_spades/pangenomics/panta/panta_results/gene_presence_absence.csv` |
| SPAdes EggNOG | `results/analysis_spades/pangenomics/eggnog/pangenomics_cohort.emapper.annotations` |
| SPAdes Rtab | `results/analysis_spades/pangenomics/panta/panta_results/gene_presence_absence.Rtab` |

### Algorithm

1. **Read taxonomy** from `pipeline_summary.tsv` (sample → family/genus/species)

2. **Per-sample protein comparison** (96 samples):
   - Read scaffold `.faa` → build set of protein sequences
   - Read contig `.faa` → any protein whose sequence is NOT in the scaffold set = **lost gene**
   - Record: sample_id, locus_tag, product (from FASTA header), contig name

3. **Load SPAdes EggNOG** → `group_to_cog` dict (same logic as new_genes_overview)

4. **Map lost genes to PANTA groups** via SPAdes gene_presence_absence.csv:
   - Build target key set: `{sample_id}-{contig}-{locus_tag}` for all lost genes
   - Stream PANTA CSV rows, check PAL* columns against target keys
   - Extract group name → map to COG via EggNOG

5. **Compute COG distributions** (global, per-species, per-genus)

6. **Pangenome baseline** from SPAdes Rtab (for comparison plot)

7. **Generate outputs** (same structure as new_genes_overview):

```
results/analysis/lost_genes/
├── lost_genes_overview.tsv          # per-sample summary
├── lost_genes_all_samples.tsv       # detail table (gene_id, product, contig, panta_group, cog_category)
├── cog_summary.tsv                  # global COG distribution + pangenome comparison
├── cog_lost_genes_distribution.png/pdf  # bar chart of COG categories
├── cog_lost_vs_pangenome.png/pdf        # comparison: lost genes vs full pangenome
├── species/
│   ├── cog_summary_{species}.tsv
│   ├── cog_{species}.png/pdf
│   └── ...
└── genus/
    ├── cog_summary_{genus}.tsv
    ├── cog_{genus}.png/pdf
    └── ...
```

### Plots (matching new_genes_overview style)
- **Plot 1**: COG distribution of lost genes (single color, no classification stacking — lost genes don't have the Other Source/Contextual/Combined classification)
- **Plot 2**: COG comparison: lost genes vs full SPAdes pangenome (side-by-side bars, proportions)
- **Plot 3**: Individual per-species COG plots
- **Plot 4**: Individual per-genus COG plots

### Notes
- Protein matching by **exact sequence comparison** (reliable since both are Bakta annotations of the same underlying DNA; minor re-annotations at boundaries are genuinely different gene calls)
- Uses the **SPAdes** PANTA/EggNOG (not the scaffold PANTA/EggNOG), since lost genes are contig-level genes
- Script is self-contained with only stdlib + matplotlib + numpy dependencies (same conda env as new_genes_overview)
- All 96 PAL* samples have both contig and scaffold .faa files