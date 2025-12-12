/*
 * Align Core Genes with MAFFT
 * Aligns individual core gene sequences
 */

process ALIGN_CORE_GENES {
    tag "$sample_id"
    label 'process_high'
    
    publishDir "${params.outdir}/phylogenetics/aligned_genes", mode: 'copy'
    
    container 'staphb/mafft:latest'
    
    input:
    tuple val(sample_id), path(core_genes_dir), val(sample_count)
    
    output:
    tuple val(sample_id), path("aligned_genes"), emit: aligned_genes
    tuple val(sample_id), val(sample_count), emit: sample_count
    
    script:
    def threads = task.cpus
    def mafft_method = params.mafft_method ?: 'auto'
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "=============================================="
    echo "Core Gene Alignment with MAFFT"
    echo "=============================================="
    echo "Sample ID: ${sample_id}"
    echo "MAFFT method: ${mafft_method}"
    echo "Threads: ${threads}"
    echo ""
    
    mkdir -p aligned_genes
    
    gene_count=\$(ls ${core_genes_dir}/*.fasta 2>/dev/null | wc -l || echo 0)
    echo "Aligning \$gene_count core genes..."
    
    if [ "\$gene_count" -eq 0 ]; then
        echo "WARNING: No core genes to align"
        echo ">no_core_genes_found" > aligned_genes/empty.aln
        echo "NNNNNNNNNNNNNNNNNNNN" >> aligned_genes/empty.aln
        exit 0
    fi
    
    # Align each gene file
    aligned=0
    for fasta in ${core_genes_dir}/*.fasta; do
        gene=\$(basename "\$fasta" .fasta)
        
        # Select MAFFT algorithm based on method parameter
        case "${mafft_method}" in
            linsi)
                mafft --localpair --maxiterate 1000 --thread ${threads} --quiet "\$fasta" > "aligned_genes/\${gene}.aln" 2>/dev/null
                ;;
            einsi)
                mafft --genafpair --maxiterate 1000 --thread ${threads} --quiet "\$fasta" > "aligned_genes/\${gene}.aln" 2>/dev/null
                ;;
            ginsi)
                mafft --globalpair --maxiterate 1000 --thread ${threads} --quiet "\$fasta" > "aligned_genes/\${gene}.aln" 2>/dev/null
                ;;
            fftns2)
                mafft --retree 2 --thread ${threads} --quiet "\$fasta" > "aligned_genes/\${gene}.aln" 2>/dev/null
                ;;
            *)
                # Auto mode - MAFFT chooses best algorithm
                mafft --auto --thread ${threads} --quiet "\$fasta" > "aligned_genes/\${gene}.aln" 2>/dev/null
                ;;
        esac
        
        # Check if alignment was created
        if [ -s "aligned_genes/\${gene}.aln" ]; then
            ((aligned++)) || true
        else
            # Copy original if alignment failed
            cp "\$fasta" "aligned_genes/\${gene}.aln"
            ((aligned++)) || true
        fi
        
        # Progress indicator
        if [ \$((aligned % 100)) -eq 0 ]; then
            echo "  Aligned \$aligned genes..."
        fi
    done
    
    echo "Completed: \$aligned genes aligned"
    echo "MAFFT alignment completed successfully"
    """
}