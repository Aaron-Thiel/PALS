/*
 * Pangenome plotting process (cohort-level)
 * Generates publication-ready pangenome plots from Panta Rtab file and EggNOG COG annotations
 * Uses pangenome-plots Docker image with all dependencies pre-installed
 */

process PANGENOME_PLOTS {
    tag "$cohort_id"

    publishDir "${params.outdir}/pangenomics/plots", mode: 'copy'

    // Custom pangenome-plots image with numpy, pandas, matplotlib, scipy, and procps
    container 'pangenome-plots:latest'

    input:
    tuple val(cohort_id), path(rtab_file), path(eggnog_annotations)

    output:
    tuple val(cohort_id), path("plots"), emit: plots_dir
    tuple val(cohort_id), path("plots/*.png"), emit: png_plots, optional: true
    tuple val(cohort_id), path("plots/*.pdf"), emit: pdf_plots, optional: true
    tuple val(cohort_id), path("plots/summary.txt"), emit: summary
    tuple val(cohort_id), path("plots/statistics.csv"), emit: statistics

    script:
    """
    echo "=============================================="
    echo "Pangenome Plotting"
    echo "=============================================="
    echo "Rtab file: ${rtab_file}"
    echo "EggNOG annotations: ${eggnog_annotations}"
    echo ""

    # Run pangenome_plots with EggNOG COG annotations
    pangenome_plots ${rtab_file} -o plots --cog ${eggnog_annotations} --style all -n 100

    echo "Pangenome plotting completed for cohort ${cohort_id}"
    """
}
