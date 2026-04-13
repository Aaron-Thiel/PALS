/*
 * IQ-TREE 3 Phylogenetic Inference
 * Builds maximum likelihood phylogenetic tree from core genome alignment
 * 
 * IQ-TREE 3 features used:
 * - Efficient ML tree search
 * - UFBoot2 ultrafast bootstrap
 * - ModelFinder for automatic model selection
 * - Partition models for gene-specific rates
 */

process IQTREE3 {
    tag "$sample_id"
    label 'process_high'
    
    publishDir "${params.outdir}/phylogenetics/tree", mode: 'copy'
    
    container 'staphb/iqtree3:latest'
    
    input:
    tuple val(sample_id), path(alignment), path(partition), val(sample_count)
    
    output:
    tuple val(sample_id), path("phylogeny.treefile"), emit: tree
    tuple val(sample_id), path("phylogeny.contree"), emit: consensus_tree, optional: true
    tuple val(sample_id), path("phylogeny.iqtree"), emit: iqtree_report
    tuple val(sample_id), path("phylogeny.log"), emit: iqtree_log
    tuple val(sample_id), path("iqtree_stats.txt"), emit: iqtree_stats
    tuple val(sample_id), path("phylogeny.splits.nex"), emit: splits, optional: true
    
    script:
    def threads = task.cpus
    def model = params.phylo_model ?: 'GTR+G'
    def use_modelfinder = params.phylo_modelfinder ?: false
    def do_bootstrap = params.phylo_bootstrap ?: true
    def bootstrap_reps = params.phylo_bootstrap_replicates ?: 1000
    def use_partition = params.phylo_use_partition ?: false
    def alrt_reps = params.phylo_alrt ?: 0  // SH-aLRT replicates (0 = disabled)
    
    // Build command options
    def model_opt = use_modelfinder ? '-m MFP' : "-m ${model}"
    def bootstrap_opt = do_bootstrap ? "-B ${bootstrap_reps}" : ""
    def alrt_opt = alrt_reps > 0 ? "-alrt ${alrt_reps}" : ""
    def partition_opt = (use_partition && partition && partition.name != 'null') ? "-p ${partition}" : ""
    """
    #!/bin/bash
    set -eu
    # Note: pipefail removed to avoid SIGPIPE issues with large files in containerized environments

    echo "=============================================="
    echo "IQ-TREE 3 Phylogenetic Analysis"
    echo "=============================================="
    echo "Sample ID: ${sample_id}"
    echo "Samples: ${sample_count}"
    echo "Model: ${model} (ModelFinder: ${use_modelfinder})"
    echo "Bootstrap: ${do_bootstrap} (${bootstrap_reps} UFBoot2 replicates)"
    echo "Partition model: ${use_partition}"
    echo "Threads: ${threads}"
    echo ""

    # Validate input alignment
    if [ ! -s "${alignment}" ]; then
        echo "ERROR: Alignment file is empty or missing"
        exit 1
    fi

    # Count sequences and get alignment length (use wc for robustness with large files)
    seq_count=\$(grep -c "^>" "${alignment}" 2>/dev/null || echo "0")
    
    # Get alignment length from first sequence
    aln_length=\$(awk '
        /^>/ { if (seq) print length(seq); seq=""; next }
        { seq = seq \$0 }
        END { if (seq) print length(seq) }
    ' "${alignment}" | head -1)
    
    echo "Input alignment:"
    echo "  - Sequences: \$seq_count"
    echo "  - Length: \$aln_length bp"
    echo ""
    
    # Check minimum sequence requirement
    if [ "\$seq_count" -lt 4 ]; then
        echo "WARNING: Only \$seq_count sequences found."
        echo "IQ-TREE requires at least 4 sequences for meaningful analysis."
        echo "Creating simple tree..."
        
        # Create minimal tree based on sequence count
        if [ "\$seq_count" -eq 3 ]; then
            s1=\$(grep "^>" "${alignment}" | sed -n '1p' | tr -d '>')
            s2=\$(grep "^>" "${alignment}" | sed -n '2p' | tr -d '>')
            s3=\$(grep "^>" "${alignment}" | sed -n '3p' | tr -d '>')
            echo "((\$s1:0.001,\$s2:0.001):0.001,\$s3:0.001);" > phylogeny.treefile
        elif [ "\$seq_count" -eq 2 ]; then
            s1=\$(grep "^>" "${alignment}" | sed -n '1p' | tr -d '>')
            s2=\$(grep "^>" "${alignment}" | sed -n '2p' | tr -d '>')
            echo "(\$s1:0.001,\$s2:0.001);" > phylogeny.treefile
        else
            s1=\$(grep "^>" "${alignment}" | head -1 | tr -d '>' || echo "sample")
            echo "(\$s1:0.0);" > phylogeny.treefile
        fi
        
        # Create placeholder files
        cp phylogeny.treefile phylogeny.contree 2>/dev/null || true
        echo "Minimal tree created due to insufficient sequences" > phylogeny.iqtree
        echo "Minimal tree created due to insufficient sequences" > phylogeny.log
        
        cat > iqtree_stats.txt << EOL
IQ-TREE 3 Analysis Statistics
=============================
Status: MINIMAL TREE (insufficient sequences)
Sequences: \$seq_count (minimum 4 required for full analysis)
Tree type: Simple distance-based placeholder
EOL
        exit 0
    fi
    
    # Run IQ-TREE 3
    echo "Running IQ-TREE 3..."
    echo ""
    
    # Note: IQ-TREE 3 uses 'iqtree3' command (or just 'iqtree' in container)
    # Check which command is available
    if command -v iqtree3 &> /dev/null; then
        IQTREE_CMD="iqtree3"
    elif command -v iqtree &> /dev/null; then
        IQTREE_CMD="iqtree"
    else
        echo "ERROR: Neither iqtree3 nor iqtree command found"
        exit 1
    fi
    
    echo "Using command: \$IQTREE_CMD"
    \$IQTREE_CMD --version || true
    echo ""
    
    # Build and run IQ-TREE command
    \$IQTREE_CMD \\
        -s "${alignment}" \\
        ${partition_opt} \\
        -pre phylogeny \\
        ${model_opt} \\
        ${bootstrap_opt} \\
        ${alrt_opt} \\
        -T ${threads} \\
        --seed 12345 \\
        -redo
    
    # Validate output
    if [ ! -f "phylogeny.treefile" ]; then
        echo "ERROR: IQ-TREE failed to produce tree file"
        echo "Check phylogeny.log for details"
        exit 1
    fi
    
    echo ""
    echo "Tree construction completed successfully"
    
    # =========================================================================
    # Generate statistics summary
    # =========================================================================
    echo ""
    echo "Generating statistics summary..."
    
    cat > iqtree_stats.txt << EOL
IQ-TREE 3 Phylogenetic Analysis Statistics
==========================================
Run ID: ${sample_id}
Timestamp: \$(date -Iseconds)

Input:
- Sequences: \$seq_count
- Alignment length: \$aln_length bp
- Partition model: ${use_partition}

Parameters:
- Substitution model: ${model}
- ModelFinder: ${use_modelfinder}
- Bootstrap: ${do_bootstrap}
- Bootstrap replicates: ${bootstrap_reps} (UFBoot2)
- SH-aLRT: ${alrt_reps > 0 ? "${alrt_reps} replicates" : "disabled"}
- Threads: ${threads}
- Random seed: 12345

EOL

    # Extract key results from IQ-TREE report
    if [ -f "phylogeny.iqtree" ]; then
        echo "Results from IQ-TREE report:" >> iqtree_stats.txt
        echo "-----------------------------" >> iqtree_stats.txt
        
        # Log-likelihood
        grep -E "^Log-likelihood of the tree:" phylogeny.iqtree >> iqtree_stats.txt 2>/dev/null || true
        
        # Best model (if ModelFinder was used)
        grep -E "^Best-fit model:" phylogeny.iqtree >> iqtree_stats.txt 2>/dev/null || true
        
        # Model used
        grep -E "^Model of substitution:" phylogeny.iqtree >> iqtree_stats.txt 2>/dev/null || true
        
        # Tree length
        grep -E "^Total tree length" phylogeny.iqtree >> iqtree_stats.txt 2>/dev/null || true
        
        # Number of parsimony-informative sites
        grep -E "parsimony-informative" phylogeny.iqtree >> iqtree_stats.txt 2>/dev/null || true
        
        # Runtime
        grep -E "^Total wall-clock time" phylogeny.iqtree >> iqtree_stats.txt 2>/dev/null || true
        grep -E "^Total CPU time" phylogeny.iqtree >> iqtree_stats.txt 2>/dev/null || true
    fi
    
    echo "" >> iqtree_stats.txt
    echo "Output files:" >> iqtree_stats.txt
    echo "-------------" >> iqtree_stats.txt
    ls -lh phylogeny.* 2>/dev/null | awk '{print "  " \$NF ": " \$5}' >> iqtree_stats.txt
    
    echo "" >> iqtree_stats.txt
    echo "Tree file: phylogeny.treefile" >> iqtree_stats.txt
    echo "Consensus tree (with bootstrap): phylogeny.contree" >> iqtree_stats.txt
    
    echo ""
    echo "=============================================="
    echo "IQ-TREE 3 analysis completed"
    echo "=============================================="
    cat iqtree_stats.txt
    """
}
