/*
 * antiSMASH process for biosynthetic gene cluster (BGC) detection
 * Runs antiSMASH on individual genomes using .fna and .gff3 from Bakta
 * Identifies and annotates secondary metabolite biosynthesis gene clusters
 * Based on nf-core/modules antismash module
 */

process ANTISMASH {
    tag "$sample_id"

    publishDir "${params.outdir}/${sample_id}/pangenomics/antismash", mode: 'copy'

    // Use nanozoo antiSMASH container (includes procps for Nextflow metrics)
    container 'nanozoo/antismash:8.0.0--b6973cb'

    input:
    tuple val(sample_id), path(fna_file), path(gff_file)

    output:
    tuple val(sample_id), path("${sample_id}"), emit: results_dir
    tuple val(sample_id), path("${sample_id}/*.gbk"), emit: genbank, optional: true
    tuple val(sample_id), path("${sample_id}/*.json"), emit: json, optional: true
    tuple val(sample_id), path("${sample_id}/index.html"), emit: html_report, optional: true
    tuple val(sample_id), path("${sample_id}_antismash_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    script:
    def extra_args = params.antismash_extra_args ?: ''
    """
    # Run antiSMASH with Bakta GFF3 annotation (standalone image has built-in databases)
    # Note: procps installed via beforeScript in nextflow.config for task metrics
    antismash \\
        "${fna_file}" \\
        --genefinding-gff3 "${gff_file}" \\
        -c ${task.cpus} \\
        --output-dir ${sample_id} \\
        --output-basename ${sample_id} \\
        --genefinding-tool none \\
        --logfile ${sample_id}/${sample_id}.log \\
        --taxon bacteria \\
        --cb-general \\
        --cb-knownclusters \\
        --cb-subclusters \\
        --asf \\
        --pfam2go \\
        --smcog-trees \\
        --cc-mibig ${extra_args}

    # Create summary table of BGCs found
    echo -e "sample\\tregion\\ttype\\tstart\\tend\\tproduct" > ${sample_id}_antismash_summary.tsv

    # Parse JSON for BGC information
    python3 << 'EOF'
import json
import glob
import os

sample = "${sample_id}"
result_dir = "${sample_id}"

json_files = glob.glob(os.path.join(result_dir, "*.json"))
for jf in json_files:
    if "regions" not in os.path.basename(jf):
        try:
            with open(jf) as f:
                data = json.load(f)
            if "records" in data:
                for record in data["records"]:
                    if "areas" in record:
                        for i, area in enumerate(record["areas"], 1):
                            products = ",".join(area.get("products", ["unknown"]))
                            start = area.get("start", "NA")
                            end = area.get("end", "NA")
                            print(f"{sample}\tregion{i}\tBGC\t{start}\t{end}\t{products}")
        except Exception as e:
            pass
EOF
    >> ${sample_id}_antismash_summary.tsv

    # Version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(echo \$(antismash --version) | sed 's/antiSMASH //;s/-.*//g')
    END_VERSIONS

    echo "antiSMASH analysis completed for ${sample_id}"
    """
}
