/*
 * Subtree Visualization Module
 *
 * Prunes the main phylogenetic tree to a subset of tips (by genus or species)
 * and generates publication-quality visualizations. This avoids re-running
 * alignment + IQ-TREE for each subset — branch lengths come from the full tree.
 *
 * Features:
 * - Prunes main tree to genus/species subsets using taxonomy data
 * - Taxonomy-based tip coloring by genus
 * - Distance scale bar (rectangular layout)
 * - Circular and rectangular layouts
 * - Multiple formats: PNG, SVG, PDF
 */

process SUBTREE_VISUALIZATION {
    tag "$group_id"
    label 'process_low'

    publishDir "${params.outdir}/phylogenetics/subtrees/${group_id}", mode: 'copy'

    container 'aaronthiel/tree-viz:latest'

    input:
    tuple val(group_id), val(tip_list)
    path(main_tree)
    path(taxonomy_csv)
    path(pipeline_summary)

    output:
    tuple val(group_id), path("${group_id}_circular.png"), path("${group_id}_rectangular.png"), emit: png, optional: true
    tuple val(group_id), path("${group_id}_circular.svg"), path("${group_id}_rectangular.svg"), emit: svg, optional: true
    tuple val(group_id), path("${group_id}_circular.pdf"), path("${group_id}_rectangular.pdf"), emit: pdf, optional: true
    tuple val(group_id), path("${group_id}_species_circular.png"), path("${group_id}_species_rectangular.png"), emit: species_png, optional: true
    tuple val(group_id), path("${group_id}_species_circular.svg"), path("${group_id}_species_rectangular.svg"), emit: species_svg, optional: true
    tuple val(group_id), path("${group_id}_species_circular.pdf"), path("${group_id}_species_rectangular.pdf"), emit: species_pdf, optional: true
    tuple val(group_id), path("subtree_stats.txt"), emit: stats, optional: true

    script:
    def width = params.tree_width ?: 3000
    def height = params.tree_height ?: 3000
    def min_tips = params.subtree_min_tips ?: 3
    """
#!/usr/bin/env python3
import toytree
import toyplot
import toyplot.png
import toyplot.svg
import toyplot.pdf
import pandas as pd
import math
import os
import sys

print("=" * 50)
print("Subtree Visualization")
print("=" * 50)
print(f"Group: ${group_id}")
print(f"Main tree: ${main_tree}")

# ========================================
# Load main tree
# ========================================
full_tree = toytree.tree("${main_tree}")
all_tips = set(full_tree.get_tip_labels())
print(f"Main tree tips: {len(all_tips)}")

# ========================================
# Parse tip list for this group
# ========================================
requested_tips = [t.strip() for t in "${tip_list}".split(",") if t.strip()]
print(f"Requested tips: {len(requested_tips)}")

tips_to_keep = [t for t in requested_tips if t in all_tips]
missing = len(requested_tips) - len(tips_to_keep)
if missing > 0:
    print(f"  Warning: {missing} tips not found in main tree")
print(f"Tips present in tree: {len(tips_to_keep)}")

min_tips = int(${min_tips})
if len(tips_to_keep) < min_tips:
    print(f"Only {len(tips_to_keep)} tips (need >= {min_tips}). Skipping visualization.")
    with open("subtree_stats.txt", "w") as f:
        f.write(f"Group: ${group_id}\\n")
        f.write(f"Status: SKIPPED (only {len(tips_to_keep)} tips, need >= {min_tips})\\n")
    sys.exit(0)

# ========================================
# Prune tree to subset
# ========================================
print("\\nPruning tree...")
tree = full_tree.mod.prune(*tips_to_keep)
print(f"Pruned tree: {tree.ntips} tips")

# ========================================
# Load taxonomy data
# ========================================
print("\\nLoading taxonomy data...")
taxonomy_data = []

species_data = []
external_tax_file = "${taxonomy_csv}"
if os.path.exists(external_tax_file):
    try:
        tax_df = pd.read_csv(external_tax_file)
        tax_df = tax_df[tax_df['genus'].notna() & (tax_df['genus'] != '')]
        if 'species' in tax_df.columns:
            species_df = tax_df[['genome_id', 'genus', 'species']].rename(columns={'genome_id': 'tip'})
            species_df = species_df[species_df['species'].notna() & (species_df['species'] != '')]
            species_data.append(species_df)
        tax_df = tax_df[['genome_id', 'genus']].rename(columns={'genome_id': 'tip'})
        taxonomy_data.append(tax_df)
        print(f"  External taxonomy: {len(tax_df)} genomes")
    except Exception as e:
        print(f"  Warning: Could not load taxonomy CSV: {e}")

summary_file = "${pipeline_summary}"
if os.path.exists(summary_file):
    try:
        summary_df = pd.read_csv(summary_file, sep='\\t')
        if 'genus' in summary_df.columns and 'sample_id' in summary_df.columns:
            internal_df = summary_df[['sample_id', 'genus']].rename(columns={'sample_id': 'tip'})
            internal_df = internal_df[internal_df['genus'].notna() & (internal_df['genus'] != '') & (internal_df['genus'] != 'NA')]
            taxonomy_data.append(internal_df)
            print(f"  Internal taxonomy: {len(internal_df)} samples")
            if 'species' in summary_df.columns:
                sp_df = summary_df[['sample_id', 'genus', 'species']].rename(columns={'sample_id': 'tip'})
                sp_df = sp_df[sp_df['species'].notna() & (sp_df['species'] != '') & (sp_df['species'] != 'NA')]
                species_data.append(sp_df)
    except Exception as e:
        print(f"  Warning: Could not load pipeline summary: {e}")

if taxonomy_data:
    taxonomy_map = pd.concat(taxonomy_data, ignore_index=True)
else:
    taxonomy_map = pd.DataFrame(columns=['tip', 'genus'])

print(f"Total taxonomy mappings: {len(taxonomy_map)}")

# ========================================
# Color mapping
# ========================================
tip_labels = tree.get_tip_labels()
tip_genus_map = taxonomy_map.set_index('tip')['genus'].to_dict()

unique_genera = sorted(set(tip_genus_map.get(t, "Unknown") for t in tip_labels) - {"Unknown"})
print(f"Unique genera in subtree: {len(unique_genera)}")

color_palette = [
    '#e6194B', '#3cb44b', '#4363d8', '#f58231', '#911eb4',
    '#42d4f4', '#f032e6', '#469990', '#9370db', '#800000',
    '#2e8b57', '#808000', '#000075', '#8b4513', '#1e90ff',
    '#ff6347', '#40e0d0', '#4b0082', '#008b8b', '#20b2aa',
    '#4169e1', '#b8860b', '#cd853f', '#0000cd', '#ba55d3',
    '#3cb371', '#7b68ee', '#c71585', '#191970', '#228b22',
    '#dc143c', '#8b0000', '#006400', '#00008b', '#8b008b'
]

genus_colors = {genus: color_palette[i % len(color_palette)]
                for i, genus in enumerate(unique_genera)}

tip_colors = []
tip_genera = []
for tip in tip_labels:
    genus = tip_genus_map.get(tip, "Unknown")
    tip_genera.append(genus)
    tip_colors.append(genus_colors.get(genus, "#808080"))

unknown_count = tip_genera.count("Unknown")
print(f"Mapped tips: {len(tip_colors)} ({unknown_count} unknown)")

genus_counts = {}
for genus in tip_genera:
    genus_counts[genus] = genus_counts.get(genus, 0) + 1

# ========================================
# Render function
# ========================================
def render_tree(layout, filename):
    print(f"  Generating {layout} tree...")

    w = ${width}
    h = ${height}

    if layout == "rectangular":
        min_height = tree.ntips * 10
        h = max(h, min_height)
        print(f"    Scaled height to {h} pixels for {tree.ntips} tips")

    if layout == "circular":
        canvas, axes, mark = tree.draw(
            layout='c',
            edge_type='c',
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

    canvas.style = {"background-color": "white"}

    # Scale bar for rectangular layout
    if layout == "rectangular":
        tree_height = tree.treenode.height
        if tree_height > 0:
            raw_interval = tree_height / 5
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

            max_scale = math.ceil(tree_height / nice_interval) * nice_interval

            scale_axes = canvas.cartesian(
                bounds=(50, w - 450, h - 80, h - 50),
                show=False
            )
            scale_axes.x.domain.min = 0
            scale_axes.x.domain.max = max_scale

            scale_axes.plot([0, tree_height], [0, 0], color="black", stroke_width=3)

            tick_height = 0.4
            num_ticks = int(max_scale / nice_interval) + 1

            for i in range(num_ticks):
                tick_x = i * nice_interval
                if tick_x <= tree_height + 0.001:
                    scale_axes.plot(
                        [tick_x, tick_x], [-tick_height, tick_height],
                        color="black", stroke_width=2
                    )
                    if nice_interval >= 1:
                        label = f"{tick_x:.0f}"
                    elif nice_interval >= 0.1:
                        label = f"{tick_x:.1f}"
                    else:
                        label = f"{tick_x:.2f}"
                    scale_axes.text(
                        tick_x, -tick_height - 0.3, label,
                        style={"font-size": "10px", "text-anchor": "middle", "fill": "black"}
                    )

            canvas.text(
                (50 + w - 450) / 2, h - 15,
                "Evolutionary distance (substitutions/site)",
                style={"font-size": "12px", "text-anchor": "middle", "fill": "black", "font-weight": "bold"}
            )

    # Legend
    legend_x = w - 400
    legend_y_start = 80

    canvas.text(
        legend_x, legend_y_start - 20,
        "Genus (sample count)",
        style={"font-size": "12px", "font-weight": "bold", "text-anchor": "start"}
    )

    for i, genus in enumerate(unique_genera):
        y_pos = legend_y_start + i * 18
        count = genus_counts.get(genus, 0)
        canvas.text(
            legend_x, y_pos,
            f"{genus} (n={count})",
            style={"font-size": "10px", "fill": genus_colors[genus], "font-weight": "bold", "text-anchor": "start"}
        )

    if unknown_count > 0:
        y_pos = legend_y_start + len(unique_genera) * 18
        canvas.text(
            legend_x, y_pos,
            f"Unknown (n={unknown_count})",
            style={"font-size": "10px", "fill": "#808080", "font-weight": "bold", "text-anchor": "start"}
        )

    toyplot.png.render(canvas, f"{filename}.png")
    toyplot.svg.render(canvas, f"{filename}.svg")
    toyplot.pdf.render(canvas, f"{filename}.pdf")
    print(f"    -> {filename}.[png,svg,pdf]")

# Generate both layouts
render_tree("circular", "${group_id}_circular")
render_tree("rectangular", "${group_id}_rectangular")

# ========================================
# Species-Level Tree (collapsed by species)
# ========================================
print("\\n" + "=" * 50)
print("Generating Species-Level Tree")
print("=" * 50)

if species_data:
    full_species_map = pd.concat(species_data, ignore_index=True)
    tip_species_map = full_species_map.set_index('tip')['species'].to_dict()
    tip_species_genus_map = full_species_map.set_index('tip')['genus'].to_dict()

    species_per_tip = {}
    for tip in tip_labels:
        species = tip_species_map.get(tip, None)
        if species:
            species_per_tip[tip] = species

    print(f"Species mapping: {len(species_per_tip)}/{len(tip_labels)} tips")

    if len(species_per_tip) > 0:
        species_to_tips = {}
        for tip, species in species_per_tip.items():
            if species not in species_to_tips:
                species_to_tips[species] = []
            species_to_tips[species].append(tip)

        print(f"Unique species: {len(species_to_tips)}")

        tips_to_keep_sp = []
        species_labels = {}
        for species, tips in species_to_tips.items():
            rep_tip = sorted(tips)[0]
            tips_to_keep_sp.append(rep_tip)
            species_labels[rep_tip] = species

        if len(tips_to_keep_sp) >= 3:
            try:
                species_tree = tree.mod.prune(*tips_to_keep_sp)
                print(f"Species tree: {species_tree.ntips} tips")

                species_tip_labels = species_tree.get_tip_labels()
                display_labels = []
                display_colors = []
                species_genus_counts = {}

                for tip in species_tip_labels:
                    if tip in species_labels:
                        species_name = species_labels[tip]
                        display_labels.append(species_name)
                        genus = tip_species_genus_map.get(tip, tip_genus_map.get(tip, "Unknown"))
                    else:
                        display_labels.append(tip)
                        genus = tip_genus_map.get(tip, "Unknown")
                    display_colors.append(genus_colors.get(genus, "#808080"))
                    species_genus_counts[genus] = species_genus_counts.get(genus, 0) + 1

                def render_species_tree(layout, filename):
                    print(f"  Generating {layout} species tree...")
                    w = ${width}
                    h = ${height}

                    if layout == "rectangular":
                        min_height = species_tree.ntips * 14
                        h = max(h, min_height)

                    if layout == "circular":
                        canvas, axes, mark = species_tree.draw(
                            layout='c', edge_type='c',
                            width=w, height=h,
                            tip_labels=display_labels,
                            tip_labels_colors=display_colors,
                            tip_labels_style={"font-size": "8px"},
                        )
                    else:
                        canvas, axes, mark = species_tree.draw(
                            layout='r',
                            width=w, height=h,
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

                    title_text = f"${group_id} - Species-Level Phylogeny (pruned from main tree)"
                    canvas.text(
                        w / 2, 30, title_text,
                        style={"font-size": "14px", "font-weight": "bold", "text-anchor": "middle"}
                    )

                    toyplot.png.render(canvas, f"{filename}.png")
                    toyplot.svg.render(canvas, f"{filename}.svg")
                    toyplot.pdf.render(canvas, f"{filename}.pdf")
                    print(f"    -> {filename}.[png,svg,pdf]")

                render_species_tree("circular", "${group_id}_species_circular")
                render_species_tree("rectangular", "${group_id}_species_rectangular")

            except Exception as e:
                print(f"Warning: Could not create species tree: {e}")
        else:
            print(f"Not enough species ({len(tips_to_keep_sp)}) for species tree")
    else:
        print("No species mapping available")
else:
    print("No species data loaded")

# ========================================
# Stats file
# ========================================
genera_list = "\\n".join([f"  - {g} (n={genus_counts.get(g, 0)})" for g in unique_genera])

with open("subtree_stats.txt", "w") as f:
    f.write(f\"\"\"Subtree Visualization Statistics
========================================
Group: ${group_id}
Source: Pruned from main tree (${main_tree})

Subtree:
- Tips: {tree.ntips}
- Unique genera: {len(unique_genera)}
- Unknown tips: {unknown_count}
- Tree height: {tree.treenode.height:.6f}

Genera:
{genera_list}
\"\"\")

print("\\n" + "=" * 50)
print("Subtree visualization completed")
print("=" * 50)
"""
}
