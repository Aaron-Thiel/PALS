/*
 * CLASSIFICATION WORKFLOW
 *
 * Progressive database building with proper caching:
 * 1. SOURMASH: Classify samples → genus/family
 * 2. GENUS_DOWNLOAD: Download references (storeDir cached per genus)
 * 3. MERGE_GFF_FASTA: Merge GFF + FASTA into proper GFF3 files
 * 4. CHECKM2_GENUS: Quality assessment (only pending genomes)
 * 5. BAKTA_GENUS: Annotate genomes missing GFF3 files
 * 6. SKANI_BUILD_DB: Build sketch database (storeDir cached per genus)
 * 7. SKANI_SEARCH: Compare each sample (uses cached db)
 * 8. REFERENCE_SELECTOR: Select best references per sample
 */

include { SOURMASH } from './classification/sourmash.nf'
include { GENUS_DOWNLOAD } from './classification/genome_download.nf'
include { MERGE_GFF_FASTA } from './classification/merge_gff_fasta.nf'
include { CHECKM2_GENUS } from './classification/checkm2_genomes.nf'
include { BAKTA_GENUS } from './classification/bakta_genus.nf'
include { SKANI_BUILD_DB; SKANI_SEARCH } from './classification/skani.nf'
include { REFERENCE_SELECTOR } from './classification/reference_selector.nf'

workflow CLASSIFICATION {
    take:
    genomes_ch  // tuple val(sample_id), path(genome_fasta)
    
    main:
    
    log.info """
    ════════════════════════════════════════════════════════════════
    CLASSIFICATION WORKFLOW
    ════════════════════════════════════════════════════════════════
    Genome database: ${params.genome_database}
    Output directory: ${params.outdir}
    Family filter: ${params.filter_family ? params.filter_target_family : 'disabled'}
    ════════════════════════════════════════════════════════════════
    """

    //=========================================================================
    // STEP 1: Taxonomic classification with sourmash
    //=========================================================================
    SOURMASH(genomes_ch)
    
    // Create sample -> taxonomy mapping
    sample_taxonomy = SOURMASH.out.genus
        .join(SOURMASH.out.family)
        .map { sample_id, genus_file, family_file ->
            def genus = genus_file.text.trim()
            def family = family_file.text.trim()
            tuple(sample_id, genus, family)
        }
        .filter { sample_id, genus, family ->
            // Validate classification
            def valid = genus && genus != "Unclassified" && genus != "Unknown" &&
                       family && family != "Unclassified" && family != "Unknown"
            
            if (!valid) {
                log.warn "Sample ${sample_id}: Could not classify (genus=${genus}, family=${family})"
                return false
            }
            
            // Family filter if enabled
            if (params.filter_family && family != params.filter_target_family) {
                log.info "Sample ${sample_id}: Filtered out (family ${family} != ${params.filter_target_family})"
                return false
            }
            
            return true
        }

    //=========================================================================
    // STEP 2: Deduplicate genera for download (one download per genus)
    //=========================================================================
    genera_to_process = sample_taxonomy
        .map { sample_id, genus, family -> tuple(genus, family, sample_id) }
        .groupTuple(by: [0, 1])  // Group by genus+family
        .map { genus, family, sample_ids ->
            // Use first sample as representative for the download
            tuple(sample_ids[0], genus, family)
        }

    //=========================================================================
    // STEP 3: Download reference genomes (storeDir cached)
    //=========================================================================
    GENUS_DOWNLOAD(genera_to_process)

    // Filter out empty genera (no genomes in NCBI)
    valid_downloads = GENUS_DOWNLOAD.out.metadata
        .join(GENUS_DOWNLOAD.out.info, by: [0, 1])
        .filter { genus, _family, metadata, _info ->
            // Check if genus has any genomes by reading the metadata
            def lines = metadata.readLines()
            def hasGenomes = lines.size() > 1  // More than just header
            if (!hasGenomes) {
                log.warn "Genus ${genus}: No genomes in NCBI, skipping downstream analysis"
            }
            return hasGenomes
        }
        .map { genus, family, metadata, _info ->
            tuple(genus, family, metadata)
        }

    //=========================================================================
    // STEP 3.5: Merge GFF + FASTA into proper GFF3 files
    //=========================================================================
    MERGE_GFF_FASTA(valid_downloads)

    //=========================================================================
    // STEP 4: Quality assessment with CheckM2 (only pending genomes)
    //=========================================================================
    CHECKM2_GENUS(MERGE_GFF_FASTA.out.metadata)

    //=========================================================================
    // STEP 5: Annotate genomes missing GFF3 with Bakta
    //=========================================================================
    BAKTA_GENUS(CHECKM2_GENUS.out.metadata)

    //=========================================================================
    // STEP 6: Build skani sketch databases (storeDir cached per genus)
    //=========================================================================
    SKANI_BUILD_DB(BAKTA_GENUS.out.metadata)

    //=========================================================================
    // STEP 7: Run skani search for each sample
    //=========================================================================
    // Create channel with sample info + genus/family
    skani_search_input = sample_taxonomy
        .join(genomes_ch)  // (sample_id, genus, family, fasta)
        .map { sample_id, genus, family, fasta ->
            tuple(sample_id, fasta, genus, family)
        }

    // Wait for skani DB to be ready, then run search
    // Only process samples whose genus has a valid skani DB
    skani_ready = skani_search_input
        .map { sample_id, fasta, genus, family -> 
            tuple(genus, sample_id, fasta, family) 
        }
        .combine(
            SKANI_BUILD_DB.out.database.map { genus, _family, database -> tuple(genus, database) },
            by: 0
        )
        .map { genus, sample_id, fasta, family, database ->
            tuple(sample_id, fasta, genus, family, database)
        }

    SKANI_SEARCH(skani_ready)

    //=========================================================================
    // STEP 8: Reference selection
    //=========================================================================
    // Join ANI results with quality metadata on genus
    reference_input = SKANI_SEARCH.out.ani_results
        .map { sample_id, genus, ani_results -> 
            tuple(genus, sample_id, ani_results) 
        }
        .combine(
            BAKTA_GENUS.out.metadata.map { genus, _family, metadata -> tuple(genus, metadata) },
            by: 0
        )
        .map { genus, sample_id, ani_results, metadata ->
            tuple(sample_id, genus, ani_results, metadata)
        }

    REFERENCE_SELECTOR(reference_input)

    //=========================================================================
    // OUTPUT CHANNELS
    //=========================================================================
    
    // Taxonomy output: sample_id, genus, species, family
    // Include all classified samples, even if their genus had no NCBI genomes
    taxonomy_output = sample_taxonomy
        .join(
            SKANI_SEARCH.out.species.map { sample_id, _genus, species_file ->
                tuple(sample_id, species_file)
            },
            remainder: true
        )
        .map { sample_id, genus, family, species_file ->
            def species = species_file ? species_file.text.trim() : "${genus}_sp"
            tuple(sample_id, genus, species, family)
        }

    emit:
    // Primary outputs
    taxonomy = taxonomy_output
    selected_references = REFERENCE_SELECTOR.out.references
    classification = SOURMASH.out.classification
    
    // Summaries for reporting
    sourmash_summary = SOURMASH.out.summary
    download_info = GENUS_DOWNLOAD.out.info
    checkm2_summary = CHECKM2_GENUS.out.summary
    bakta_summary = BAKTA_GENUS.out.summary
    skani_summary = SKANI_SEARCH.out.summary
    selection_summary = REFERENCE_SELECTOR.out.summary
}
