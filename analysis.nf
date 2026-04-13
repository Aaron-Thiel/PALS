#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    PALS - Analysis Pipeline
========================================================================================
    Modular analysis pipeline for pangenomics, phylogenetics, and BGC analysis.

    Architecture:
      CORE (always runs):
        BAKTA → MERGE_EXTERNAL → PANTA_COHORT (auto-detects cache)

      DOWNSTREAM BLOCKS (independent, optional):
        --codon_qc_enable       Codon bias QC
        --pangenomics_enable    EGGNOG + pangenome plots
        --phylogenetics_enable  Core gene tree
        --bgc_enable            antiSMASH BGC detection
        --new_genes_enable      GenomeViz scaffold vs reference

      SUBSET ANALYSIS (--subset_analysis_enable):
        When enabled, adds genus/species-level analysis within each block:
        pangenomics  → subset Rtab + plots per genus/species
        phylogenetics → subtree pruning per genus/species
        bgc          → BiG-SCAPE networks per genus/species
        new_genes    → aggregated overview per family/genus/species

    Run with: nextflow run analysis.nf -c analysis.config --internal results/assembly
    SPAdes:   nextflow run analysis.nf -c analysis.config --internal results/assembly --spades

    Author: Aaron Thiel
========================================================================================
*/

// =========================================================================
// Import core processes
// =========================================================================

include { BAKTA } from './modules/bakta'
include { PANTA_COHORT } from './modules/panta_cohort.nf'

// =========================================================================
// Import downstream blocks
// =========================================================================

include { CODON_BLOCK } from './workflows/codon_bias'
include { PANGENOMICS_BLOCK } from './workflows/pangenomics'
include { PHYLOGENETICS_BLOCK } from './workflows/phylogenetics'
include { BGC_BLOCK } from './workflows/bgc'
include { NEW_GENES_BLOCK } from './workflows/new_genes'

// =========================================================================
// Main Workflow
// =========================================================================

workflow {

    log.info """
    ===================================
    PALS - Analysis Pipeline
    ===================================
    Input:  ${params.internal}${params.external ? ' + ' + params.external : ''}
    Mode:   ${params.spades ? 'SPAdes (pre-existing Bakta)' : 'Normal (PASA scaffolds → Bakta)'}
    Output: ${params.outdir}

    Blocks:
      Codon QC:       ${params.codon_qc_enable}
      Codon Correlation: ${params.codon_correlation_enable}
      Pangenomics:    ${params.pangenomics_enable}
      Phylogenetics:  ${params.phylogenetics_enable}
      BGC:            ${params.bgc_enable}
      New Genes:      ${params.new_genes_enable && !params.spades}
      Subset analysis: ${params.subset_analysis_enable}
      HQ Plots:       ${params.hq_plots_enable}${params.hq_plots_enable ? ' (min ' + params.scaffold_min_completeness + '%)' : ''}
    ===================================
    """.stripIndent()

    // =====================================================================
    // Shared input channels
    // =====================================================================

    ch_internal_taxonomy = channel
        .fromPath("${params.internal}/pipeline_summary.tsv", checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.sample_id, row.family ?: 'Unknown', row.genus ?: 'Unknown', row.species ?: 'Unknown')
        }

    ch_pipeline_summary = channel
        .fromPath("${params.internal}/pipeline_summary.tsv", checkIfExists: true)

    ch_taxonomy_csv = channel
        .fromPath(params.taxonomy_csv, checkIfExists: true)

    // =====================================================================
    // CORE: Resolve Bakta annotations
    // =====================================================================

    if (params.spades) {
        ch_internal_gff = channel.fromPath("${params.internal}/*/bakta/*.gff3", checkIfExists: true)
            .map { gff -> tuple(gff.baseName, gff) }
        ch_internal_fna = channel.fromPath("${params.internal}/*/bakta/*.fna", checkIfExists: true)
            .map { fna -> tuple(fna.baseName, fna) }
        ch_internal_ffn = channel.fromPath("${params.internal}/*/bakta/*.ffn", checkIfExists: true)
            .map { ffn -> tuple(ffn.baseName, ffn) }
        ch_internal_faa = channel.fromPath("${params.internal}/*/bakta/*.faa", checkIfExists: true)
            .map { faa -> tuple(faa.baseName, faa) }
    } else {
        ch_scaffolds = channel
            .fromPath("${params.internal}/*/pasa_filter/selected_scaffold.fasta", checkIfExists: true)
            .map { fasta -> tuple(fasta.parent.parent.name, fasta) }

        BAKTA(ch_scaffolds)

        ch_internal_gff = BAKTA.out.gff
        ch_internal_fna = BAKTA.out.fna
        ch_internal_ffn = BAKTA.out.nucleotides
        ch_internal_faa = BAKTA.out.proteins
    }

    // =====================================================================
    // CORE: Merge external annotations
    // =====================================================================

    if (params.external) {
        ch_external_gff = channel
            .fromPath(["${params.external}/*.gff3", "${params.external}/**/*.gff3"], checkIfExists: false)
            .map { gff -> tuple(gff.baseName, gff) }
        ch_external_fna = channel
            .fromPath(["${params.external}/*.fna", "${params.external}/**/*.fna"], checkIfExists: false)
            .map { fna -> tuple(fna.baseName, fna) }
        ch_all_gff = ch_internal_gff.mix(ch_external_gff)
        ch_all_fna = ch_internal_fna.mix(ch_external_fna)
    } else {
        ch_all_gff = ch_internal_gff
        ch_all_fna = ch_internal_fna
    }

    // =====================================================================
    // CORE: PANTA pangenome analysis (auto-detect cache)
    // =====================================================================

    def cached_rtab = file("${params.outdir}/pangenomics/panta/panta_results/gene_presence_absence.Rtab")

    if (cached_rtab.exists()) {
        log.info "Found cached PANTA results at ${params.outdir}/pangenomics/panta/"
        def panta_path = file("${params.outdir}/pangenomics/panta/panta_results")
        def sample_count = cached_rtab.readLines()[0].split('\t').length - 1

        ch_panta_dir    = channel.of(tuple("pangenomics_cohort", panta_path))
        ch_rtab         = channel.of(tuple("pangenomics_cohort", cached_rtab))
        ch_sample_count = channel.of(tuple("pangenomics_cohort", sample_count))
    } else {
        ch_all_gff
            .map { _sample_id, gff -> gff }
            .collect()
            .map { gff_list -> tuple("pangenomics_cohort", gff_list, gff_list.size()) }
            .set { ch_panta_input }

        PANTA_COHORT(ch_panta_input)

        ch_panta_dir    = PANTA_COHORT.out.panta_dir
        ch_rtab         = PANTA_COHORT.out.rtab
        ch_sample_count = PANTA_COHORT.out.sample_count
    }

    // =====================================================================
    // DOWNSTREAM BLOCKS
    // =====================================================================

    if (params.codon_qc_enable) {
        ch_gpa_csv = ch_panta_dir
            .map { _cohort_id, panta_dir -> file("${panta_dir}/gene_presence_absence.csv") }

        CODON_BLOCK(
            ch_internal_ffn, ch_internal_taxonomy,
            ch_rtab,
            ch_pipeline_summary.first(), ch_taxonomy_csv.first(),
            ch_gpa_csv.first()
        )
    }

    if (params.pangenomics_enable) {
        PANGENOMICS_BLOCK(
            ch_panta_dir, ch_rtab,
            ch_pipeline_summary.first(), ch_taxonomy_csv.first()
        )
    }

    if (params.phylogenetics_enable) {
        PHYLOGENETICS_BLOCK(
            ch_panta_dir, ch_sample_count, ch_internal_taxonomy,
            ch_taxonomy_csv.first(), ch_pipeline_summary.first()
        )
    }

    if (params.bgc_enable) {
        BGC_BLOCK(ch_internal_fna, ch_internal_gff, ch_pipeline_summary.first())
    }

    if (params.new_genes_enable && !params.spades) {
        // Resolve reference paths from pipeline_summary
        ch_genomeviz_taxonomy = channel
            .fromPath("${params.internal}/pipeline_summary.tsv", checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.sample_id, row.genus ?: '', row.species ?: '', row.best_reference ?: '') }
            .filter { _id, genus, _sp, ref ->
                genus && genus != 'NA' && genus != '' && ref && ref != 'NA' && ref != ''
            }

        // Contig annotations (assembly bakta, pre-PASA)
        ch_contig_fna = channel.fromPath("${params.internal}/*/bakta/*.fna", checkIfExists: true)
            .map { fna -> tuple(fna.baseName, fna) }
        ch_contig_gff = channel.fromPath("${params.internal}/*/bakta/*.gff3", checkIfExists: true)
            .map { gff -> tuple(gff.baseName, gff) }
        ch_contig_faa = channel.fromPath("${params.internal}/*/bakta/*.faa", checkIfExists: true)
            .map { faa -> tuple(faa.baseName, faa) }

        // Build 9-tuple for GENE_COMPARISON (scaffold vs contig)
        ch_genomeviz_taxonomy
            .map { sample_id, genus, species, best_ref ->
                def ref_fna = file("${params.reference_genomes_dir}/${genus}/${species}/${best_ref}.fna")
                def ref_gff = file("${params.reference_genomes_dir}/${genus}/${species}/${best_ref}.gff3")
                tuple(sample_id, ref_fna, ref_gff)
            }
            .filter { _id, fna, gff -> fna.exists() && gff.exists() }
            .join(ch_internal_fna).join(ch_internal_gff).join(ch_internal_faa)
            .join(ch_contig_fna).join(ch_contig_gff).join(ch_contig_faa)
            .map { sample_id, ref_fna, ref_gff, scaffold_fna, scaffold_gff, scaffold_faa,
                   contig_fna, contig_gff, contig_faa ->
                tuple(sample_id, scaffold_fna, scaffold_gff, scaffold_faa,
                      ref_fna, ref_gff, contig_fna, contig_gff, contig_faa)
            }
            .set { ch_genomeviz_input }

        // Build skani + metadata channel for SELECT_COMPLETE_REFERENCE
        // genome_dir is the original database path (not staged) for resolving relative paths in metadata
        ch_genomeviz_taxonomy
            .map { sample_id, genus, _species, _best_ref ->
                def ani_tsv = file("${params.internal}/${sample_id}/classification/skani/ani_results.tsv")
                def metadata_csv = file("${params.reference_genomes_dir}/${genus}/genome_metadata.csv")
                def genome_dir = "${params.reference_genomes_dir}/${genus}"
                (ani_tsv.exists() && metadata_csv.exists()) ? tuple(sample_id, ani_tsv, metadata_csv, genome_dir) : null
            }
            .filter { item -> item != null }
            .set { ch_skani_metadata }

        // Derive PANTA gene_presence_absence.csv from panta_dir
        ch_gene_presence_absence_csv = ch_panta_dir
            .map { _cohort_id, panta_dir -> file("${panta_dir}/gene_presence_absence.csv") }

        // EggNOG annotations (auto-detect from cache, fallback to NO_EGGNOG)
        def cached_eggnog = file("${params.outdir}/pangenomics/eggnog/pangenomics_cohort.emapper.annotations")
        if (cached_eggnog.exists()) {
            ch_eggnog_for_newgenes = channel.of(cached_eggnog)
        } else {
            log.warn "  [new_genes] EggNOG annotations not found at ${cached_eggnog}. COG analysis will be skipped."
            ch_eggnog_for_newgenes = channel.of(file("NO_EGGNOG"))
        }

        // Rtab file for pangenome-wide COG baseline
        ch_rtab_file = ch_rtab.map { _cohort_id, rtab -> rtab }

        NEW_GENES_BLOCK(
            ch_genomeviz_input,
            ch_skani_metadata,
            ch_pipeline_summary.first(),
            ch_gene_presence_absence_csv.first(),
            ch_eggnog_for_newgenes.first(),
            ch_rtab_file.first()
        )
    }
}
