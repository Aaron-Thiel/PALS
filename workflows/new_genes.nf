/*
 * New Genes Block
 *
 * Two independent steps:
 *   1. GENE_COMPARISON (mode=comparison): scaffold vs contig comparison
 *      for each sample → produces gene_comparison_report.csv.
 *      Feeds into NEW_GENES_OVERVIEW and LOST_GENES (aggregated analyses).
 *
 *   2. GENOME_ALIGNMENT_VIZ (mode=scaffold): scaffold vs nearest complete
 *      reference genome → produces clean alignment PNGs for publication.
 *      Only runs for samples where a complete genome exists at species level.
 *      SELECT_COMPLETE_REFERENCE picks the best complete genome from skani ANI
 *      results joined with genome_metadata.csv.
 */

include { GENOME_VISUALIZATION as GENE_COMPARISON } from '../modules/new_genes/genomeviz.nf'
include { GENOME_VISUALIZATION as GENOME_ALIGNMENT_VIZ } from '../modules/new_genes/genomeviz.nf'
include { SELECT_COMPLETE_REFERENCE } from '../modules/new_genes/select_complete_reference.nf'
include { NEW_GENES_OVERVIEW } from '../modules/new_genes/new_genes_overview.nf'
include { LOST_GENES } from '../modules/new_genes/lost_genes.nf'

workflow NEW_GENES_BLOCK {

    take:
    ch_genomeviz_input        // tuple(sample_id, scaffold_fna, scaffold_gff, scaffold_faa,
                              //       ref_fna, ref_gff, contig_fna, contig_gff, contig_faa)
    ch_skani_metadata         // tuple(sample_id, ani_results, genome_metadata, genome_dir)
    ch_pipeline_summary       // path (value channel)
    ch_gene_presence_absence  // path: gene_presence_absence.csv
    ch_eggnog_annotations     // path: eggnog .emapper.annotations
    ch_rtab                   // path: gene_presence_absence.Rtab

    main:
    // =================================================================
    // Step 1: Per-sample gene comparison (scaffold vs contig)
    // =================================================================
    ch_genomeviz_input
        .map { sample_id, scaffold_fna, scaffold_gff, scaffold_faa,
               ref_fna, ref_gff, contig_fna, contig_gff, contig_faa ->
            tuple(sample_id, "comparison", scaffold_fna, scaffold_gff, scaffold_faa,
                  ref_fna, ref_gff, contig_fna, contig_gff, contig_faa)
        }
        .set { ch_comparison_input }

    GENE_COMPARISON(ch_comparison_input)

    // =================================================================
    // Step 2: Find nearest complete genome per sample
    // =================================================================
    SELECT_COMPLETE_REFERENCE(ch_skani_metadata)

    // =================================================================
    // Step 3: Scaffold alignment against complete genome
    // =================================================================
    // Extract scaffold + contig files per sample (for stageAs requirements)
    ch_sample_files = ch_genomeviz_input
        .map { sample_id, scaffold_fna, scaffold_gff, scaffold_faa,
               ref_fna, ref_gff, contig_fna, contig_gff, contig_faa ->
            tuple(sample_id, scaffold_fna, scaffold_gff, scaffold_faa,
                  contig_fna, contig_gff, contig_faa)
        }

    // Join complete reference paths with sample files (only samples with complete genomes)
    SELECT_COMPLETE_REFERENCE.out.ref_paths
        .map { sample_id, ref_fna_str, ref_gff_str ->
            tuple(sample_id, file(ref_fna_str), file(ref_gff_str))
        }
        .join(ch_sample_files)
        .map { sample_id, ref_fna, ref_gff, scaffold_fna, scaffold_gff, scaffold_faa,
               contig_fna, contig_gff, contig_faa ->
            tuple(sample_id, "all", scaffold_fna, scaffold_gff, scaffold_faa,
                  ref_fna, ref_gff, contig_fna, contig_gff, contig_faa)
        }
        .set { ch_alignment_input }

    GENOME_ALIGNMENT_VIZ(ch_alignment_input)

    // =================================================================
    // Aggregated analyses (require subset_analysis_enable)
    // =================================================================
    if (params.subset_analysis_enable) {
        GENE_COMPARISON.out.results
            .map { sample_id, results_dir -> results_dir }
            .collect()
            .set { ch_all_comparison }

        NEW_GENES_OVERVIEW(
            ch_all_comparison,
            ch_pipeline_summary,
            ch_gene_presence_absence,
            ch_eggnog_annotations,
            ch_rtab
        )

        // --- Lost genes analysis (contig genes absent from scaffolds) ---
        ch_genomeviz_input
            .map { sample_id, scaffold_fna, scaffold_gff, scaffold_faa,
                   ref_fna, ref_gff, contig_fna, contig_gff, contig_faa ->
                contig_faa
            }
            .collect()
            .set { ch_contig_faa }

        ch_genomeviz_input
            .map { sample_id, scaffold_fna, scaffold_gff, scaffold_faa,
                   ref_fna, ref_gff, contig_fna, contig_gff, contig_faa ->
                contig_gff
            }
            .collect()
            .set { ch_contig_gff }

        ch_genomeviz_input
            .map { sample_id, scaffold_fna, scaffold_gff, scaffold_faa,
                   ref_fna, ref_gff, contig_fna, contig_gff, contig_faa ->
                scaffold_faa
            }
            .collect()
            .set { ch_scaffold_faa }

        LOST_GENES(
            ch_contig_faa,
            ch_contig_gff,
            ch_scaffold_faa,
            ch_pipeline_summary,
            ch_gene_presence_absence,
            ch_eggnog_annotations,
            ch_rtab
        )
    }

    emit:
    comparison_results = GENE_COMPARISON.out.results
    alignment_results  = GENOME_ALIGNMENT_VIZ.out.results
}
