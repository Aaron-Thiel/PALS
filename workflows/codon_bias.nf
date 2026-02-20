/*
 * Codon Bias Analysis Block
 *
 * Phase 1: Reference building — builds species/genus codon references from external genomes
 * Phase 2: Internal QC — per-sample codon bias QC against references
 * Phase 4: Correlation — codon bias vs pangenome categories (family-level, internal samples)
 * Phase 5: Subset correlation — per genus/species (gated on subset_analysis_enable)
 *
 * References are cached via storeDir and auto-detected on re-runs.
 */

// --- Existing processes (moved from modules/codon.nf) ---
include { EXTRACT_CDS } from '../modules/codon/build_reference'
include { EXTRACT_CDS as EXTRACT_CDS_GENUS } from '../modules/codon/build_reference'
include { BUILD_SPECIES_REFERENCE } from '../modules/codon/build_reference'
include { BUILD_GENUS_REFERENCE } from '../modules/codon/build_reference'
include { CODON_QC } from '../modules/codon/qc_analysis'
include { AGGREGATE_CODON_QC } from '../modules/codon/aggregate'

// --- New processes ---
include { EXTERNAL_CODON_METRICS } from '../modules/codon/external_metrics'
include { CODON_CORRELATION } from '../modules/codon/correlation'
include { CODON_CORRELATION as SUBSET_CODON_CORRELATION } from '../modules/codon/correlation'

workflow CODON_BLOCK {

    take:
    ch_internal_ffn       // tuple(sample_id, ffn)
    ch_internal_taxonomy  // tuple(sample_id, family, genus, species)
    ch_rtab               // tuple(cohort_id, rtab)
    ch_pipeline_summary   // path (value channel)
    ch_taxonomy_csv       // path (value channel)
    ch_gpa_csv            // path (value channel) - gene_presence_absence.csv

    main:

    // =========================================================================
    // Build meta map: join FFN with taxonomy
    // =========================================================================

    ch_codon_samples = ch_internal_ffn
        .join(ch_internal_taxonomy)
        .map { sample_id, ffn, family, genus, species ->
            def f = (family && family != 'NA' && family != '') ? family : 'Unknown'
            def g = (genus && genus != 'NA' && genus != '') ? genus : 'Unknown'
            def s = (species && species != 'NA' && species != '') ? species : 'Unknown'
            [[id: sample_id, family: f, genus: g, species: s], ffn]
        }

    // =========================================================================
    // Phase 1: Reference building (from modules/codon.nf Steps 1-4)
    // =========================================================================

    // Quality thresholds
    def min_completeness = params.codon_ref_min_completeness ?: 90.0
    def max_contamination = params.codon_ref_max_contamination ?: 5.0
    def min_species_refs = params.codon_min_species_refs ?: 10

    // Load genome metadata from all genera
    def family_dir = file(params.reference_genomes_dir)

    // Build lookup of existing references and HQ genomes
    def species_refs = [:]
    def genus_refs = [:]
    def species_genomes = [:]
    def genus_genomes = [:]

    if (family_dir.exists()) {
        family_dir.eachDir { genus_dir ->
            def genus = genus_dir.name

            def genus_ref = file("${genus_dir}/codon_reference.json")
            if (genus_ref.exists()) {
                genus_refs[genus] = genus_ref
            }

            def metadata_file = file("${genus_dir}/genome_metadata.csv")
            def hq_genomes_by_species = [:]

            if (metadata_file.exists()) {
                metadata_file.eachLine { line, idx ->
                    if (idx == 1) return

                    def fields = line.split(',')
                    if (fields.size() >= 10) {
                        def accession = fields[0]
                        def completeness
                        def contamination
                        try {
                            completeness = fields[4] as Double
                            contamination = fields[5] as Double
                        } catch (NumberFormatException e) {
                            return
                        }

                        def fna_path = fields[8]
                        def species_folder = fna_path.replaceAll('^\\./', '').split('/')[0]

                        if (completeness >= min_completeness && contamination <= max_contamination) {
                            def fna_file = file("${genus_dir}/${species_folder}/${accession}.fna")
                            def gff_file = file("${genus_dir}/${species_folder}/${accession}.gff3")

                            if (fna_file.exists() && gff_file.exists()) {
                                if (!hq_genomes_by_species.containsKey(species_folder)) {
                                    hq_genomes_by_species[species_folder] = []
                                }
                                hq_genomes_by_species[species_folder] << [fna_file, gff_file]

                                if (!genus_genomes.containsKey(genus)) {
                                    genus_genomes[genus] = []
                                }
                                genus_genomes[genus] << [fna_file, gff_file]
                            }
                        }
                    }
                }
            }

            genus_dir.eachDir { species_dir ->
                def species = species_dir.name
                def species_ref = file("${species_dir}/codon_reference.json")
                if (species_ref.exists()) {
                    species_refs[species] = species_ref
                }
                if (hq_genomes_by_species.containsKey(species)) {
                    species_genomes[species] = hq_genomes_by_species[species]
                }
            }
        }
    }

    log.info "  [codon] Found ${species_refs.size()} species-level and ${genus_refs.size()} genus-level codon references"
    log.info "  [codon] HQ genomes available: ${species_genomes.size()} species, ${genus_genomes.size()} genera"

    // Determine which references need building
    ch_species_needed = ch_codon_samples
        .map { meta, _ffn -> meta.species }
        .filter { it && it != 'Unknown' && it != 'NA' }
        .unique()

    ch_genera_needed = ch_codon_samples
        .map { meta, _ffn -> meta.genus }
        .filter { it && it != 'Unknown' && it != 'NA' }
        .unique()

    ch_species_to_build = ch_species_needed
        .filter { species -> !species_refs.containsKey(species) && species_genomes.containsKey(species) }
        .map { species ->
            def genomes = species_genomes[species]
            [species, genomes.collect { it[0] }, genomes.collect { it[1] }]
        }

    ch_genera_to_build = ch_genera_needed
        .filter { genus -> !genus_refs.containsKey(genus) && genus_genomes.containsKey(genus) }
        .map { genus ->
            def genomes = genus_genomes[genus]
            [genus, genomes.collect { it[0] }, genomes.collect { it[1] }]
        }

    // Extract CDS and build missing references
    EXTRACT_CDS(ch_species_to_build)
    BUILD_SPECIES_REFERENCE(EXTRACT_CDS.out.cds)

    EXTRACT_CDS_GENUS(ch_genera_to_build)
    BUILD_GENUS_REFERENCE(EXTRACT_CDS_GENUS.out.cds)

    // Collect newly built references
    ch_new_species_refs = BUILD_SPECIES_REFERENCE.out.reference
        .toList()
        .map { list -> def m = [:]; list.each { item -> m[item[0]] = item[1] }; m }

    ch_new_genus_refs = BUILD_GENUS_REFERENCE.out.reference
        .toList()
        .map { list -> def m = [:]; list.each { item -> m[item[0]] = item[1] }; m }

    // =========================================================================
    // Phase 2: Internal QC (match samples with refs, run CODON_QC, aggregate)
    // =========================================================================

    ch_qc_input = ch_codon_samples
        .combine(ch_new_species_refs)
        .combine(ch_new_genus_refs)
        .map { meta, ffn, new_species_map, new_genus_map ->
            def species = meta.species
            def genus = meta.genus
            def species_genome_count = species_genomes[species]?.size() ?: 0
            def all_species_refs = species_refs + new_species_map
            def all_genus_refs = genus_refs + new_genus_map
            def ref = null
            def ref_level = 'none'

            if (all_species_refs[species] && species_genome_count >= min_species_refs) {
                ref = all_species_refs[species]
                ref_level = 'species'
            } else if (all_genus_refs[genus]) {
                ref = all_genus_refs[genus]
                ref_level = 'genus'
            }

            if (ref) {
                [meta + [ref_level: ref_level], ffn, file(ref.toString())]
            } else {
                log.warn "No codon reference for ${meta.id} (species: ${species}, genus: ${genus})"
                [meta + [ref_level: 'none'], ffn, file('NO_FILE')]
            }
        }

    CODON_QC(ch_qc_input)

    AGGREGATE_CODON_QC(
        CODON_QC.out.summary
            .map { meta, f -> [meta.id, f] }
            .groupTuple()
            .map { _id, files -> files.first() }
            .collect()
    )

    // Log failures
    CODON_QC.out.report
        .filter { meta, report -> report.text.contains('"qc_status": "FAIL"') }
        .subscribe { meta, _report -> log.warn "Codon QC FAILED: ${meta.id}" }

    // =========================================================================
    // Phase 3: Reference genome codon metrics (from external annotations)
    // =========================================================================

    // Load taxonomy CSV at construction time for external genome genus mapping
    def taxonomy_map = [:]
    def taxonomy_file = file(params.taxonomy_csv)
    if (taxonomy_file.exists()) {
        taxonomy_file.eachLine { line, idx ->
            if (idx == 1) return
            def fields = line.split(',')
            if (fields.size() >= 4) {
                taxonomy_map[fields[0].trim()] = [genus: fields[2].trim()]
            }
        }
    }

    // Discover external FFN files grouped by genus
    def ext_dir = params.external_annotations_dir
    def ext_ffn_by_genus = [:]
    if (ext_dir) {
        taxonomy_map.each { genome_id, tax ->
            def ffn = file("${ext_dir}/${genome_id}/${genome_id}.ffn")
            if (ffn.exists()) {
                def genus = tax.genus
                if (!ext_ffn_by_genus.containsKey(genus)) {
                    ext_ffn_by_genus[genus] = []
                }
                ext_ffn_by_genus[genus] << ffn
            }
        }
    }

    log.info "  [codon] External FFN files: ${ext_ffn_by_genus.values().sum { it.size() } ?: 0} genomes across ${ext_ffn_by_genus.size()} genera"

    ch_ext_by_genus = Channel.fromList(
        ext_ffn_by_genus.collect { genus, ffns -> tuple(genus, ffns) }
    ).filter { genus, ffns -> ffns.size() > 0 }

    ch_ext_input = ch_ext_by_genus
        .combine(ch_new_genus_refs)
        .map { genus, ffns, new_genus_map ->
            def all_genus_refs_combined = genus_refs + new_genus_map
            def ref = all_genus_refs_combined[genus]
            ref ? tuple(genus, ffns, file(ref.toString())) : null
        }
        .filter { it != null }

    EXTERNAL_CODON_METRICS(ch_ext_input)

    // =========================================================================
    // Phase 4: Family-level correlation (all internal samples)
    // =========================================================================

    if (params.codon_correlation_enable) {

        // Collect internal gene_codon_analysis.tsv files
        ch_internal_codon_files = CODON_QC.out.gene_analysis
            .map { meta, f -> f }
            .collect()

        // Collect reference codon metrics TSVs
        ch_reference_codon_files = EXTERNAL_CODON_METRICS.out.metrics
            .map { genus, f -> f }
            .collect()
            .ifEmpty([file('NO_REFERENCE_DATA')])

        // Family-level correlation (all internal samples, no taxonomy filter)
        ch_family_scope = Channel.of(
            tuple("family", "internal", "correlation")
        )

        def rare_thresh = params.codon_rare_threshold ?: 0.15
        def core_thresh = params.codon_core_threshold ?: 0.95

        CODON_CORRELATION(
            ch_family_scope,
            ch_rtab.map { _id, rtab -> rtab }.first(),
            ch_gpa_csv,
            ch_internal_codon_files,
            ch_reference_codon_files,
            ch_pipeline_summary,
            ch_taxonomy_csv,
            rare_thresh,
            core_thresh
        )

        // =================================================================
        // Phase 5: Subset correlation (per genus/species with internal samples)
        // =================================================================

        if (params.subset_analysis_enable) {

            // Build list of qualifying groups from internal taxonomy
            def min_samples = params.subset_min_samples_pangenomics ?: 3
            def internal_genera = [:]   // genus -> count of internal samples
            def internal_species = [:]  // species -> count of internal samples

            // Collect internal taxonomy at construction time
            def ps_file = file("${params.internal}/pipeline_summary.tsv")
            if (ps_file.exists()) {
                ps_file.eachLine { line, idx ->
                    if (idx == 1) return
                    def fields = line.split('\t')
                    if (fields.size() >= 4) {
                        def genus = fields[2]?.trim()
                        def species = fields[3]?.trim()
                        if (genus && genus != 'NA' && genus != '') {
                            internal_genera[genus] = (internal_genera[genus] ?: 0) + 1
                            if (species && species != 'NA' && species != '') {
                                def sp_key = "${genus}_${species}"
                                internal_species[sp_key] = (internal_species[sp_key] ?: 0) + 1
                            }
                        }
                    }
                }
            }

            // Build qualifying group list (only count internal samples)
            def subset_groups = []
            internal_genera.each { genus, count ->
                if (count >= min_samples) {
                    subset_groups << ["genus_${genus}", genus, null]
                }
            }
            internal_species.each { sp_key, count ->
                def parts = sp_key.split('_', 2)
                def genus = parts[0]
                def species = parts.size() > 1 ? parts[1] : ''
                if (count >= min_samples) {
                    subset_groups << ["species_${sp_key}", genus, species]
                }
            }

            log.info "  [codon] Subset correlation: ${subset_groups.size()} groups (${internal_genera.size()} genera, ${internal_species.size()} species with internal samples)"

            if (subset_groups) {
                // One correlation per group (internal samples only)
                ch_subset_scopes = Channel.fromList(
                    subset_groups.collect { group ->
                        def group_id = group[0]
                        tuple(group_id, "internal", "subsets/${group_id}")
                    }
                )

                SUBSET_CODON_CORRELATION(
                    ch_subset_scopes,
                    ch_rtab.map { _id, rtab -> rtab }.first(),
                    ch_gpa_csv,
                    ch_internal_codon_files,
                    ch_reference_codon_files,
                    ch_pipeline_summary,
                    ch_taxonomy_csv,
                    rare_thresh,
                    core_thresh
                )
            }
        }
    }

    emit:
    summary       = CODON_QC.out.summary
    report        = CODON_QC.out.report
    batch_summary = AGGREGATE_CODON_QC.out.batch_summary
    batch_report  = AGGREGATE_CODON_QC.out.batch_report
}
