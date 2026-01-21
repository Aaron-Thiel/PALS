/*
 * Codon Bias QC Workflow
 *
 * Analyzes codon usage patterns to detect assembly anomalies.
 * Compares internal samples against species-level references (preferred)
 * or genus-level references (fallback when species has <10 reference genomes).
 *
 * References are built from high-quality database genomes:
 *   - Completeness >= 90%
 *   - Contamination <= 5%
 *
 * Directory structure:
 *   {genomes_dir}/{family}/{genus}/{species}/codon_reference.json  (species-level)
 *   {genomes_dir}/{family}/{genus}/codon_reference.json            (genus-level)
 */

include { EXTRACT_CDS } from './codon/build_reference'
include { EXTRACT_CDS as EXTRACT_CDS_GENUS } from './codon/build_reference'
include { BUILD_SPECIES_REFERENCE } from './codon/build_reference'
include { BUILD_GENUS_REFERENCE } from './codon/build_reference'
include { CODON_QC } from './codon/qc_analysis'
include { AGGREGATE_CODON_QC } from './codon/aggregate'

workflow CODON_BIAS_QC {
    take:
    ch_samples  // tuple(meta, ffn) where meta has [id, genus, species, family]

    main:

    // =========================================================================
    // Step 1: Scan for existing references and determine what needs building
    // =========================================================================

    // Quality thresholds
    def min_completeness = params.codon_ref_min_completeness ?: 90.0
    def max_contamination = params.codon_ref_max_contamination ?: 5.0
    def min_species_refs = params.codon_min_species_refs ?: 10  // Fall back to genus if fewer refs

    // Load genome metadata from all genera
    def genomes_base = file(params.reference_genomes_dir).parent  // /BGC-data/databases/genomes
    def family_dir = file(params.reference_genomes_dir)  // /BGC-data/databases/genomes/Lactobacillaceae

    // Build lookup of existing references and HQ genomes
    def species_refs = [:]      // species -> reference file path
    def genus_refs = [:]        // genus -> reference file path
    def species_genomes = [:]   // species -> list of [fna, gff3] for HQ genomes
    def genus_genomes = [:]     // genus -> list of [fna, gff3] for HQ genomes

    if (family_dir.exists()) {
        family_dir.eachDir { genus_dir ->
            def genus = genus_dir.name

            // Check for genus-level reference
            def genus_ref = file("${genus_dir}/codon_reference.json")
            if (genus_ref.exists()) {
                genus_refs[genus] = genus_ref
            }

            // Load metadata for this genus
            def metadata_file = file("${genus_dir}/genome_metadata.csv")
            def hq_genomes_by_species = [:]  // species_folder -> list of [fna, gff3]

            if (metadata_file.exists()) {
                metadata_file.eachLine { line, idx ->
                    if (idx == 1) return  // Skip header (eachLine index starts at 1)

                    def fields = line.split(',')
                    if (fields.size() >= 10) {
                        def accession = fields[0]
                        def species_name = fields[1]  // Full species name

                        // Parse numeric fields with error handling
                        def completeness
                        def contamination
                        try {
                            completeness = fields[4] as Double
                            contamination = fields[5] as Double
                        } catch (NumberFormatException e) {
                            return  // Skip rows with invalid numbers
                        }

                        def assembly_level = fields[7]
                        def fna_path = fields[8]
                        def gff_path = fields[9]

                        // Extract species folder name from fna_path (e.g., "./Limosilactobacillus_reuteri/GCA_xxx.fna")
                        def species_folder = fna_path.replaceAll('^\\./', '').split('/')[0]

                        // Check quality thresholds
                        if (completeness >= min_completeness && contamination <= max_contamination) {
                            // Build full paths
                            def fna_file = file("${genus_dir}/${species_folder}/${accession}.fna")
                            def gff_file = file("${genus_dir}/${species_folder}/${accession}.gff3")

                            if (fna_file.exists() && gff_file.exists()) {
                                // Add to species-level collection
                                if (!hq_genomes_by_species.containsKey(species_folder)) {
                                    hq_genomes_by_species[species_folder] = []
                                }
                                hq_genomes_by_species[species_folder] << [fna_file, gff_file]

                                // Add to genus-level collection
                                if (!genus_genomes.containsKey(genus)) {
                                    genus_genomes[genus] = []
                                }
                                genus_genomes[genus] << [fna_file, gff_file]
                            }
                        }
                    }
                }
            }

            // Check for species-level references and store genome lists
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

    log.info "Found ${species_refs.size()} species-level and ${genus_refs.size()} genus-level codon references"
    log.info "HQ genomes available: ${species_genomes.size()} species, ${genus_genomes.size()} genera"

    // =========================================================================
    // Step 2: Determine which references need building
    // =========================================================================

    // Get unique species and genera from samples
    ch_species_needed = ch_samples
        .map { meta, _ffn -> meta.species }
        .filter { it && it != 'Unknown' && it != 'NA' }
        .unique()

    ch_genera_needed = ch_samples
        .map { meta, _ffn -> meta.genus }
        .filter { it && it != 'Unknown' && it != 'NA' }
        .unique()

    // Species that need building (not already cached AND have HQ genomes)
    ch_species_to_build = ch_species_needed
        .filter { species ->
            !species_refs.containsKey(species) && species_genomes.containsKey(species)
        }
        .map { species ->
            def genomes = species_genomes[species]
            def fna_files = genomes.collect { it[0] }
            def gff_files = genomes.collect { it[1] }
            [species, fna_files, gff_files]
        }

    // Genera that need building (not already cached AND have HQ genomes)
    ch_genera_to_build = ch_genera_needed
        .filter { genus ->
            !genus_refs.containsKey(genus) && genus_genomes.containsKey(genus)
        }
        .map { genus ->
            def genomes = genus_genomes[genus]
            def fna_files = genomes.collect { it[0] }
            def gff_files = genomes.collect { it[1] }
            [genus, fna_files, gff_files]
        }

    // =========================================================================
    // Step 3: Extract CDS and build missing references
    // =========================================================================

    // Extract CDS for species that need building
    EXTRACT_CDS(ch_species_to_build)
    BUILD_SPECIES_REFERENCE(EXTRACT_CDS.out.cds)

    // Extract CDS for genera that need building
    EXTRACT_CDS_GENUS(ch_genera_to_build)
    BUILD_GENUS_REFERENCE(EXTRACT_CDS_GENUS.out.cds)

    // Collect newly built references into maps
    // These will be combined with pre-existing refs for sample matching
    ch_new_species_refs = BUILD_SPECIES_REFERENCE.out.reference
        .toList()
        .map { list ->
            def m = [:]
            list.each { item -> m[item[0]] = item[1] }
            m
        }

    ch_new_genus_refs = BUILD_GENUS_REFERENCE.out.reference
        .toList()
        .map { list ->
            def m = [:]
            list.each { item -> m[item[0]] = item[1] }
            m
        }

    // =========================================================================
    // Step 4: Match samples with references (species preferred, genus fallback)
    // =========================================================================

    // Combine pre-existing refs (Groovy variables) with newly built refs (channels)
    // Wait for all references to be built before matching samples
    ch_qc_input = ch_samples
        .combine(ch_new_species_refs)
        .combine(ch_new_genus_refs)
        .map { meta, ffn, new_species_map, new_genus_map ->
            def species = meta.species
            def genus = meta.genus

            // Check species genome count
            def species_genome_count = species_genomes[species]?.size() ?: 0

            // Merge pre-existing refs with newly built refs
            def all_species_refs = species_refs + new_species_map
            def all_genus_refs = genus_refs + new_genus_map

            // Try species-level first, but only if we have enough reference genomes
            def ref = null
            def ref_level = 'none'

            if (all_species_refs[species] && species_genome_count >= min_species_refs) {
                ref = all_species_refs[species]
                ref_level = 'species'
            } else if (all_genus_refs[genus]) {
                // Fall back to genus-level (either no species ref or too few genomes)
                ref = all_genus_refs[genus]
                ref_level = 'genus'
                if (all_species_refs[species] && species_genome_count < min_species_refs) {
                    log.info "Using genus-level ref for ${meta.id}: species ${species} has only ${species_genome_count} refs (min: ${min_species_refs})"
                }
            }

            if (ref) {
                // Ensure ref is a proper file path
                [meta + [ref_level: ref_level], ffn, file(ref.toString())]
            } else {
                log.warn "No codon reference for ${meta.id} (species: ${species}, genus: ${genus})"
                [meta + [ref_level: 'none'], ffn, file('NO_FILE')]
            }
        }

    // =========================================================================
    // Step 5: Run codon QC analysis
    // =========================================================================

    CODON_QC(ch_qc_input)

    // =========================================================================
    // Step 6: Aggregate results
    // =========================================================================

    AGGREGATE_CODON_QC(
        CODON_QC.out.summary
            .map { meta, f -> [meta.id, f] }
            .groupTuple()
            .map { _id, files -> files.first() }
            .collect()
    )

    // Log failures
    CODON_QC.out.report
        .filter { meta, report ->
            report.text.contains('"qc_status": "FAIL"')
        }
        .subscribe { meta, _report ->
            log.warn "Codon QC FAILED: ${meta.id}"
        }

    emit:
    summary       = CODON_QC.out.summary
    report        = CODON_QC.out.report
    batch_summary = AGGREGATE_CODON_QC.out.batch_summary
    batch_report  = AGGREGATE_CODON_QC.out.batch_report
}
