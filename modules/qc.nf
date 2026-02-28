/*
 * Comprehensive QC mini-workflow
 * Calls individual tool modules and aggregates results with MultiQC
 * This creates a nested workflow approach
 */

// Import individual QC tool modules
include { BUSCO } from './qc/busco.nf'
include { CHECKM2 } from './qc/checkm2.nf'
include { QUAST } from './qc/quast.nf'

// MultiQC aggregation process
process MULTIQC_AGGREGATE {
    tag "$sample_id"
    
    publishDir "${params.outdir}/${sample_id}/qc_${stage_name}", mode: 'copy'
    
    container 'staphb/multiqc:latest'
    
    input:
    tuple val(sample_id), path(qc_results, stageAs: 'qc_results/*')
    val stage_name
    
    output:
    tuple val(sample_id), path("multiqc_report.html"), emit: multiqc_report
    tuple val(sample_id), path("multiqc_report_data"), emit: multiqc_data
    tuple val(sample_id), path("qc_summary.txt"), emit: summary
    tuple val(sample_id), path("*.json"), emit: json_reports, optional: true
    
    script:
    """
    # Create comprehensive summary
    cat > qc_summary.txt << EOF
Comprehensive QC Analysis for ${sample_id}
Stage: ${stage_name}
==========================================

Assemblies analyzed: \$(find qc_results -name "*.fasta" | wc -l)
EOF
    
    # Create JSON summaries for MultiQC
    python3 << 'EOF'
import json
import os
import glob
import re

def parse_busco_summary(file_path, assembly_name):
    data = {"assembly_name": assembly_name, "stage": "${stage_name}"}
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        
        complete_match = re.search(r'C:([0-9.]+)%', content)
        single_match = re.search(r'S:([0-9.]+)%', content)
        duplicated_match = re.search(r'D:([0-9.]+)%', content)
        fragmented_match = re.search(r'F:([0-9.]+)%', content)
        missing_match = re.search(r'M:([0-9.]+)%', content)
        
        if complete_match:
            data['complete_percent'] = float(complete_match.group(1))
        if single_match:
            data['single_copy_percent'] = float(single_match.group(1))
        if duplicated_match:
            data['duplicated_percent'] = float(duplicated_match.group(1))
        if fragmented_match:
            data['fragmented_percent'] = float(fragmented_match.group(1))
        if missing_match:
            data['missing_percent'] = float(missing_match.group(1))
            
    except Exception as e:
        print(f"Error parsing BUSCO {file_path}: {e}")
    
    return data

# Collect all QC data
all_data = []

# Parse BUSCO results - search recursively for short_summary files
busco_summaries = glob.glob("qc_results/**/short_summary*.txt", recursive=True)
for summary_file in busco_summaries:
    # Extract assembly name from path
    assembly_name = os.path.basename(os.path.dirname(summary_file))
    if assembly_name.startswith('run_'):
        # If in run_lineage directory, go up one level
        assembly_name = os.path.basename(os.path.dirname(os.path.dirname(summary_file)))
    data = parse_busco_summary(summary_file, assembly_name)
    data['tool'] = 'busco'
    all_data.append(data)

# Save comprehensive JSON with stage name to avoid file name collisions
with open("qc_comprehensive_summary_${stage_name}.json", 'w') as f:
    json.dump(all_data, f, indent=2)

print("Comprehensive QC JSON summary created for MultiQC")
EOF
    
    # Create MultiQC config
    cat > multiqc_config.yaml << EOF
title: "Comprehensive QC Report - ${sample_id}"
subtitle: "Stage: ${stage_name}"
intro_text: "This report aggregates quality control metrics from multiple tools across all analyzed assemblies."

report_header_info:
  - Contact E-mail: 'user@example.com'
  - Pipeline: 'PALS Comprehensive QC'
  - Stage: '${stage_name}'

module_order:
  - busco
  - checkm2
  - quast
  - custom_content
EOF
    
    # Run MultiQC with explicit module list and better file discovery
    multiqc \\
        --config multiqc_config.yaml \\
        --title "Comprehensive QC Report - ${sample_id} (${stage_name})" \\
        --filename multiqc_report.html \\
        --force \\
        --verbose \\
        --module busco \\
        --module checkm2 \\
        --module quast \\
        qc_results/
    
    echo "Comprehensive QC analysis with MultiQC aggregation completed for ${sample_id} at stage ${stage_name}"
    """
}

// Main QC workflow
workflow QC {
    take:
    assemblies_ch     // tuple val(sample_id), path(assemblies)
    stage_name        // val
    
    main:
    // Prepare assembly channel for individual tools
    assemblies_ch
        .map { sample_id, assemblies ->
            def assembly_list = assemblies instanceof List ? assemblies : [assemblies]
            tuple(sample_id, assembly_list)
        }
        .set { ch_assemblies }
    
    // Run individual QC tools in parallel
    BUSCO(ch_assemblies, stage_name)
    
    // Run CHECKM2 and QUAST on each assembly individually
    ch_assemblies
        .flatMap { sample_id, assemblies ->
            assemblies.collect { assembly ->
                tuple(sample_id, assembly)
            }
        }
        .set { ch_single_assemblies }
    
    CHECKM2(ch_single_assemblies, stage_name)

    // For QUAST, we need to check if there's a reference genome available
    // This will be enhanced later to use MASH results for reference selection
    // Create an empty reference file
    def empty_ref = file("${workflow.workDir}/empty_reference.fasta")
    if (!empty_ref.exists()) {
        empty_ref.text = ""
    }

    QUAST(ch_single_assemblies, channel.value(empty_ref), stage_name)
    
    // Collect all QC results from all tools
    ch_assemblies
        .join(BUSCO.out.results)
        .join(CHECKM2.out.checkm2_results.groupTuple())
        .join(QUAST.out.results.groupTuple())
        .map { sample_id, _assemblies, busco_results, checkm2_results, quast_results ->
            // Create a directory structure that MultiQC expects
            def all_results = []
            all_results.add(busco_results)
            all_results.addAll(checkm2_results)
            all_results.addAll(quast_results)
            tuple(sample_id, all_results.flatten())
        }
        .set { ch_collected_results }
    
    // Run MultiQC aggregation
    MULTIQC_AGGREGATE(ch_collected_results, stage_name)
    
    emit:
    results = MULTIQC_AGGREGATE.out.multiqc_report
    summary = MULTIQC_AGGREGATE.out.summary
    multiqc_report = MULTIQC_AGGREGATE.out.multiqc_report
    multiqc_data = MULTIQC_AGGREGATE.out.multiqc_data
    json_reports = MULTIQC_AGGREGATE.out.json_reports
    checkm2_json = CHECKM2.out.json_reports  // Direct CheckM2 JSON for HQ filtering
    quast_tsv = QUAST.out.tsv_reports         // QUAST transposed report (L90, N50 etc.)
}