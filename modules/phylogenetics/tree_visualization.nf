/*
 * Tree Visualization Module (toytree)
 * Generates publication-quality phylogenetic tree visualizations
 *
 * Features:
 * - Taxonomy-based tip coloring by genus
 * - Bootstrap support values on nodes
 * - Distance scale bar showing branch length units
 * - Circular and rectangular layouts
 * - Multiple formats: PNG, SVG, PDF
 */

process TREE_VISUALIZATION {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.outdir}/phylogenetics/tree/visualizations", mode: 'copy'

    container 'aaronthiel/tree-viz:latest'

    input:
    tuple val(sample_id), path(tree_file), path(internal_genus_map)

    output:
    tuple val(sample_id), path("${sample_id}_circular.png"), path("${sample_id}_rectangular.png"), emit: png, optional: true
    tuple val(sample_id), path("${sample_id}_circular.svg"), path("${sample_id}_rectangular.svg"), emit: svg, optional: true
    tuple val(sample_id), path("${sample_id}_circular.pdf"), path("${sample_id}_rectangular.pdf"), emit: pdf, optional: true
    tuple val(sample_id), path("${sample_id}_species_circular.png"), path("${sample_id}_species_rectangular.png"), emit: species_png, optional: true
    tuple val(sample_id), path("${sample_id}_species_circular.svg"), path("${sample_id}_species_rectangular.svg"), emit: species_svg, optional: true
    tuple val(sample_id), path("${sample_id}_species_circular.pdf"), path("${sample_id}_species_rectangular.pdf"), emit: species_pdf, optional: true
    tuple val(sample_id), path("tree_stats.txt"), emit: stats
    path("taxonomy_mapping.tsv"), emit: taxonomy_map

    script:
    def width = params.tree_width ?: 20
    def height = params.tree_height ?: 20
    """
#!/usr/bin/env python3
import toytree
import toyplot
import toyplot.png
import toyplot.svg
import toyplot.pdf
import pandas as pd
import os

print("=" * 50)
print("Tree Visualization with Genus Coloring")
print("=" * 50)
print(f"Sample ID: ${sample_id}")
print(f"Tree file: ${tree_file}")

# Load tree
tree = toytree.tree("${tree_file}")
print(f"Loaded tree with {tree.ntips} tips")

# ========================================
# Load taxonomy data
# ========================================
print("\\nLoading taxonomy data...")
taxonomy_data = []

# Load external data from taxonomy.csv (GCA_* genomes)
external_tax_file = '/home/azureuser/BGC-link/master_thesis/external_data/taxonomy.csv'
species_data = []  # For species-level tree
if os.path.exists(external_tax_file):
    try:
        tax_df = pd.read_csv(external_tax_file)
        tax_df = tax_df[tax_df['genus'].notna() & (tax_df['genus'] != '')]
        # Keep species column for species-level tree
        if 'species' in tax_df.columns:
            species_df = tax_df[['genome_id', 'genus', 'species']].rename(columns={'genome_id': 'tip'})
            species_df = species_df[species_df['species'].notna() & (species_df['species'] != '')]
            species_data.append(species_df)
        tax_df = tax_df[['genome_id', 'genus']].rename(columns={'genome_id': 'tip'})
        taxonomy_data.append(tax_df)
        print(f"  Loaded {len(tax_df)} external genomes from taxonomy.csv")
    except Exception as e:
        print(f"  Warning: Could not load taxonomy.csv: {e}")

# Load internal data from internal_genus_map.tsv (parsed from pipeline_summary.tsv)
internal_genus_file = "${internal_genus_map}"
if os.path.exists(internal_genus_file):
    try:
        internal_df = pd.read_csv(internal_genus_file, sep='\\t', header=None, names=['tip', 'genus'])
        internal_df = internal_df[internal_df['genus'].notna() & (internal_df['genus'] != '') & (internal_df['genus'] != 'NA')]
        taxonomy_data.append(internal_df)
        print(f"  Loaded {len(internal_df)} internal samples from internal_genus_map.tsv")
    except Exception as e:
        print(f"  Warning: Could not load internal_genus_map.tsv: {e}")

# Combine taxonomy data
if taxonomy_data:
    taxonomy_map = pd.concat(taxonomy_data, ignore_index=True)
else:
    taxonomy_map = pd.DataFrame(columns=['tip', 'genus'])

print(f"Total taxonomy mappings: {len(taxonomy_map)}")

# Save taxonomy mapping
taxonomy_map.to_csv("taxonomy_mapping.tsv", sep="\\t", index=False)

# ========================================
# Create color mapping for genera
# ========================================
unique_genera = sorted(taxonomy_map['genus'].unique())
print(f"Unique genera: {len(unique_genera)}")

# Distinct, saturated color palette (no white or light colors)
color_palette = [
    '#e6194B', '#3cb44b', '#4363d8', '#f58231', '#911eb4',
    '#42d4f4', '#f032e6', '#469990', '#9370db', '#800000',
    '#2e8b57', '#808000', '#000075', '#8b4513', '#1e90ff',
    '#ff6347', '#40e0d0', '#4b0082', '#008b8b', '#20b2aa',
    '#4169e1', '#b8860b', '#cd853f', '#0000cd', '#ba55d3',
    '#3cb371', '#7b68ee', '#c71585', '#191970', '#228b22',
    '#dc143c', '#8b0000', '#006400', '#00008b', '#8b008b'
]

# Create genus to color mapping
genus_colors = {genus: color_palette[i % len(color_palette)]
                for i, genus in enumerate(unique_genera)}

# Map tip labels to colors
tip_labels = tree.get_tip_labels()
tip_genus_map = taxonomy_map.set_index('tip')['genus'].to_dict()

tip_colors = []
tip_genera = []
for tip in tip_labels:
    genus = tip_genus_map.get(tip, "Unknown")
    tip_genera.append(genus)
    if genus == "Unknown":
        tip_colors.append("#808080")  # Gray for unknown
    else:
        tip_colors.append(genus_colors[genus])

unknown_count = tip_genera.count("Unknown")
print(f"Mapped tips: {len(tip_colors)} ({unknown_count} unknown)")

# Count samples per genus
genus_counts = {}
for genus in tip_genera:
    genus_counts[genus] = genus_counts.get(genus, 0) + 1

# Create simple visualizations
def render_tree(layout, filename):
    if tree.ntips < 3:
        print(f"  Skipping {layout} tree: only {tree.ntips} tip(s), need at least 3")
        return
    print(f"  Generating {layout} tree...")

    # Draw tree - circular needs edge_type='c'
    # Scale dimensions based on tree size
    w = ${width}
    h = ${height}

    # For rectangular layout, scale height with number of tips (min 10 pixels per tip)
    if layout == "rectangular":
        min_height = tree.ntips * 10
        h = max(h, min_height)
        print(f"    Scaled height to {h} pixels for {tree.ntips} tips")

    if layout == "circular":
        canvas, axes, mark = tree.draw(
            layout='c',
            edge_type='c',  # Required for circular layout
            width=w,
            height=h,
            tip_labels=tip_labels,
            tip_labels_colors=tip_colors,
            tip_labels_style={"font-size": "8px"},
        )
    else:
        canvas, axes, mark = tree.draw(
            layout='r',
            width=w,
            height=h,
            tip_labels=tip_labels,
            tip_labels_colors=tip_colors,
            tip_labels_style={"font-size": "8px"},
        )

    # Set white background for all export formats
    canvas.style = {"background-color": "white"}

    # Add scale bar for rectangular trees only (showing branch length units)
    if layout == "rectangular":
        # Get tree dimensions for calculating scale bar length
        tree_height = tree.treenode.height  # Total tree length from root to tips

        # Calculate nice tick intervals based on tree height
        import math
        
        # Find a nice round interval (0.1, 0.2, 0.5, 1.0, 2.0, etc.)
        raw_interval = tree_height / 5  # Aim for ~5 tick marks
        magnitude = 10 ** math.floor(math.log10(raw_interval))
        normalized = raw_interval / magnitude
        
        if normalized < 1.5:
            nice_interval = 1 * magnitude
        elif normalized < 3.5:
            nice_interval = 2 * magnitude
        elif normalized < 7.5:
            nice_interval = 5 * magnitude
        else:
            nice_interval = 10 * magnitude
        
        # Round tree_height up to next nice interval for clean end point
        max_scale = math.ceil(tree_height / nice_interval) * nice_interval
        
        # Create a separate axes for the scale bar at bottom-left
        # Position below the tree, aligned with tree x-coordinates
        scale_axes = canvas.cartesian(
            bounds=(50, w - 450, h - 80, h - 50),  # Leave room for legend on right
            show=False
        )
        
        # Set the x-axis range to match tree coordinates
        scale_axes.x.domain.min = 0
        scale_axes.x.domain.max = max_scale
        
        # Draw the main scale bar line spanning full tree width
        scale_axes.plot(
            [0, tree_height],
            [0, 0],
            color="black",
            stroke_width=3
        )
        
        # Add tick marks at regular intervals
        tick_height = 0.4
        num_ticks = int(max_scale / nice_interval) + 1
        
        for i in range(num_ticks):
            tick_x = i * nice_interval
            if tick_x <= tree_height + 0.001:  # Only draw ticks within tree range
                # Draw vertical tick
                scale_axes.plot(
                    [tick_x, tick_x],
                    [-tick_height, tick_height],
                    color="black",
                    stroke_width=2
                )
                
                # Add tick label
                # Format label nicely (remove unnecessary decimals)
                if nice_interval >= 1:
                    label = f"{tick_x:.0f}"
                elif nice_interval >= 0.1:
                    label = f"{tick_x:.1f}"
                else:
                    label = f"{tick_x:.2f}"
                    
                scale_axes.text(
                    tick_x, -tick_height - 0.3,
                    label,
                    style={"font-size": "10px", "text-anchor": "middle", "fill": "black"}
                )
        
        # Add final tick at exact tree height if not already there
        last_tick = (num_ticks - 1) * nice_interval
        if abs(tree_height - last_tick) > nice_interval * 0.1:
            scale_axes.plot(
                [tree_height, tree_height],
                [-tick_height, tick_height],
                color="black",
                stroke_width=2
            )

        # Add scale bar title/label centered below
        canvas.text(
            (50 + w - 450) / 2,  # Center of scale bar area
            h - 15,
            "Evolutionary distance (substitutions/site)",
            style={"font-size": "12px", "text-anchor": "middle", "fill": "black", "font-weight": "bold"}
        )

    # Add legend with colored text showing all genera with counts
    legend_x = w - 400  # Increased width for longer text with counts
    legend_y_start = 80

    canvas.text(
        legend_x, legend_y_start - 20,
        "Genus (sample count)",
        style={"font-size": "12px", "font-weight": "bold", "text-anchor": "start"}
    )

    # Show all genera with their counts
    for i, genus in enumerate(unique_genera):
        y_pos = legend_y_start + i * 18
        count = genus_counts.get(genus, 0)
        canvas.text(
            legend_x, y_pos,
            f"{genus} (n={count})",
            style={"font-size": "10px", "fill": genus_colors[genus], "font-weight": "bold", "text-anchor": "start"}
        )

    # Add unknown count if present
    if unknown_count > 0:
        y_pos = legend_y_start + len(unique_genera) * 18
        canvas.text(
            legend_x, y_pos,
            f"Unknown (n={unknown_count})",
            style={"font-size": "10px", "fill": "#808080", "font-weight": "bold", "text-anchor": "start"}
        )

    # Save outputs
    toyplot.png.render(canvas, f"{filename}.png")
    toyplot.svg.render(canvas, f"{filename}.svg")
    toyplot.pdf.render(canvas, f"{filename}.pdf")
    print(f"    -> {filename}.[png,svg,pdf]")

# Generate both layouts
render_tree("circular", "${sample_id}_circular")
render_tree("rectangular", "${sample_id}_rectangular")

# ========================================
# Generate Species-Level Tree (collapsed by species, colored by genus)
# ========================================
print("\\n" + "=" * 50)
print("Generating Species-Level Tree")
print("=" * 50)

# Build species mapping from collected data
if species_data:
    full_species_map = pd.concat(species_data, ignore_index=True)
    tip_species_map = full_species_map.set_index('tip')['species'].to_dict()
    tip_species_genus_map = full_species_map.set_index('tip')['genus'].to_dict()

    # Get species for each tip
    species_per_tip = {}
    for tip in tip_labels:
        species = tip_species_map.get(tip, None)
        if species:
            species_per_tip[tip] = species

    print(f"Species mapping available for {len(species_per_tip)}/{len(tip_labels)} tips")

    if len(species_per_tip) > 0:
        # Group tips by species
        species_to_tips = {}
        for tip, species in species_per_tip.items():
            if species not in species_to_tips:
                species_to_tips[species] = []
            species_to_tips[species].append(tip)

        print(f"Found {len(species_to_tips)} unique species")

        # For species with multiple genomes, keep one representative
        # Strategy: keep the tip that's most central in the clade (or just first alphabetically for simplicity)
        tips_to_keep = []
        species_labels = {}  # map from kept tip to species name

        for species, tips in species_to_tips.items():
            # Keep first tip alphabetically as representative
            rep_tip = sorted(tips)[0]
            tips_to_keep.append(rep_tip)
            species_labels[rep_tip] = species

        # Only keep tips with species mapping (omit unknown from species tree)
        n_unmapped = len(tip_labels) - len(species_per_tip)
        print(f"Keeping {len(tips_to_keep)} species tips (omitting {n_unmapped} unmapped)")

        if len(tips_to_keep) >= 2:
            # Create pruned tree with only representative tips
            try:
                species_tree = tree.mod.prune(*tips_to_keep)
                print(f"Created species tree with {species_tree.ntips} tips")

                # Get new tip labels and create display labels (species names)
                species_tip_labels = species_tree.get_tip_labels()
                display_labels = []
                display_colors = []
                species_genus_counts = {}

                for tip in species_tip_labels:
                    if tip in species_labels:
                        # Use species name as display label
                        species_name = species_labels[tip]
                        display_labels.append(species_name)
                        # Get genus for coloring
                        genus = tip_species_genus_map.get(tip, tip_genus_map.get(tip, "Unknown"))
                    else:
                        # Keep original label for unmapped tips
                        display_labels.append(tip)
                        genus = tip_genus_map.get(tip, "Unknown")

                    if genus == "Unknown":
                        display_colors.append("#808080")
                    else:
                        display_colors.append(genus_colors.get(genus, "#808080"))

                    # Count species per genus
                    species_genus_counts[genus] = species_genus_counts.get(genus, 0) + 1

                # Render species tree
                def render_species_tree(layout, filename):
                    print(f"  Generating {layout} species tree...")

                    w = ${width}
                    h = ${height}

                    if layout == "rectangular":
                        min_height = species_tree.ntips * 14
                        h = max(h, min_height)

                    if layout == "circular":
                        w = max(w, 600)
                        h = max(h, 600)

                        canvas, axes, mark = species_tree.draw(
                            layout='c',
                            edge_type='c',
                            width=w,
                            height=h,
                            tip_labels=display_labels,
                            tip_labels_colors=display_colors,
                            tip_labels_style={"font-size": "8px"},
                        )
                    else:
                        canvas, axes, mark = species_tree.draw(
                            layout='r',
                            width=w,
                            height=h,
                            tip_labels=display_labels,
                            tip_labels_colors=display_colors,
                            tip_labels_style={"font-size": "9px"},
                        )

                    canvas.style = {"background-color": "white"}

                    # Legend
                    legend_x = w - 380
                    legend_y_start = 80

                    canvas.text(
                        legend_x, legend_y_start - 20,
                        "Genus (species count)",
                        style={"font-size": "12px", "font-weight": "bold", "text-anchor": "start"}
                    )

                    sorted_genera = sorted([g for g in species_genus_counts.keys() if g != "Unknown"])
                    for i, genus in enumerate(sorted_genera):
                        y_pos = legend_y_start + i * 18
                        count = species_genus_counts[genus]
                        canvas.text(
                            legend_x, y_pos,
                            f"{genus} (n={count})",
                            style={"font-size": "10px", "fill": genus_colors.get(genus, "#808080"), "font-weight": "bold", "text-anchor": "start"}
                        )

                    if "Unknown" in species_genus_counts:
                        y_pos = legend_y_start + len(sorted_genera) * 18
                        canvas.text(
                            legend_x, y_pos,
                            f"Unknown (n={species_genus_counts['Unknown']})",
                            style={"font-size": "10px", "fill": "#808080", "font-weight": "bold", "text-anchor": "start"}
                        )

                    canvas.text(
                        w / 2, 30,
                        "Species-Level Phylogeny (collapsed by species)",
                        style={"font-size": "14px", "font-weight": "bold", "text-anchor": "middle"}
                    )

                    toyplot.png.render(canvas, f"{filename}.png")
                    toyplot.svg.render(canvas, f"{filename}.svg")
                    toyplot.pdf.render(canvas, f"{filename}.pdf")
                    print(f"    -> {filename}.[png,svg,pdf]")

                render_species_tree("circular", "${sample_id}_species_circular")
                render_species_tree("rectangular", "${sample_id}_species_rectangular")

            except Exception as e:
                print(f"Warning: Could not create species tree: {e}")
        else:
            print("Not enough tips to create species tree")
    else:
        print("No species mapping available, skipping species-level tree")
else:
    print("No species data loaded, skipping species-level tree")

# Create stats file
genera_list = "\\n".join([f"  - {g}" for g in unique_genera])

# Species tree info
species_tree_info = ""
if species_data and len(species_per_tip) > 0:
    species_tree_info = f\"\"\"
Species-Level Tree:
- Species with mapping: {len(species_per_tip)}
- Unique species: {len(species_to_tips)}
- Tips in species tree: {len(tips_to_keep)}
\"\"\"

with open("tree_stats.txt", "w") as f:
    f.write(f\"\"\"Tree Visualization Statistics
========================================
Sample ID: ${sample_id}
Tree file: ${tree_file}

Genome-Level Tree:
- Tips: {tree.ntips}
- Unique genera: {len(unique_genera)}
- Taxonomy mappings: {len(taxonomy_map)}
- Unknown tips: {unknown_count}
{species_tree_info}
Genera:
{genera_list}
\"\"\")

print("=" * 50)
print("Visualization completed")
print("=" * 50)
"""
}
