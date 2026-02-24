/*
 * BiG-SCAPE Network Plot — custom publication-ready BGC similarity network
 *
 * Parses BiG-SCAPE 2 output (network files, record annotations, clustering)
 * and generates a multi-panel figure organized by BGC class with nodes
 * colored by genus/species.
 */

process BIGSCAPE_NETWORK_PLOT {
    tag "$group_id"
    label 'process_low'

    publishDir "${params.outdir}/pangenomics/bigscape/${group_id}", mode: 'copy'

    conda "conda-forge::networkx conda-forge::matplotlib conda-forge::pandas conda-forge::scipy"
    container null

    input:
    tuple val(group_id), path(bigscape_results)
    path(pipeline_summary)

    output:
    tuple val(group_id), path("bgc_network_*.png"), emit: png, optional: true
    tuple val(group_id), path("bgc_network_*.pdf"), emit: pdf, optional: true

    script:
    """
    #!/usr/bin/env python3

    import os
    import sys
    import glob
    import csv
    import re
    import sqlite3
    import networkx as nx
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from collections import defaultdict
    import math

    group_id = "${group_id}"
    results_dir = "${bigscape_results}"
    summary_file = "${pipeline_summary}"

    print("=" * 60)
    print(f"BiG-SCAPE Network Plot -- {group_id}")
    print("=" * 60)

    # ================================================================
    # Step 1: Build sample_id -> genus/species mapping
    # ================================================================
    sample_genus = {}
    sample_species = {}
    with open(summary_file) as f:
        reader = csv.DictReader(f, delimiter="\\t")
        for row in reader:
            sid = row["sample_id"]
            genus = row.get("genus", "")
            species = row.get("species", "")
            if genus and genus != "NA":
                sample_genus[sid] = genus
            if species and species != "NA":
                sample_species[sid] = species

    # ================================================================
    # Step 2: Find BiG-SCAPE output files
    # ================================================================
    annotation_file = None
    network_files = []
    clustering_files = []

    for root, dirs, files in os.walk(results_dir):
        for fname in files:
            fp = os.path.join(root, fname)
            if fname == "record_annotations.tsv":
                annotation_file = fp
            elif "annotation" in fname.lower() and fname.endswith(".tsv") and not annotation_file:
                annotation_file = fp
            elif fname.endswith(".network") and "topolink" not in fname and fname != "full.network":
                network_files.append(fp)
            elif fname.startswith("clustering_") and fname.endswith(".tsv"):
                clustering_files.append(fp)

    # Fallback: use full.network if no filtered network found
    if not network_files:
        for root, dirs, files in os.walk(results_dir):
            for fname in files:
                if fname.endswith(".network"):
                    network_files.append(os.path.join(root, fname))

    print(f"Annotation file: {annotation_file}")
    print(f"Network files: {len(network_files)}")
    print(f"Clustering files: {len(clustering_files)}")

    if not annotation_file or not network_files:
        print("WARNING: Missing BiG-SCAPE output files. No plot generated.")
        sys.exit(0)

    # ================================================================
    # Step 3: Parse record annotations
    # ================================================================
    bgc_info = {}

    with open(annotation_file) as f:
        reader = csv.DictReader(f, delimiter="\\t")
        headers = reader.fieldnames
        print(f"Annotation columns: {headers}")

        # Build case-insensitive column lookup
        col_map = {h.lower(): h for h in headers}

        for row in reader:
            # Case-insensitive field access helper
            def get_ci(key, default=""):
                actual = col_map.get(key.lower())
                return row[actual] if actual and row.get(actual) else default

            bgc_name = get_ci("Record") or get_ci("BGC") or get_ci("bgc_id") or row.get(headers[0], "")

            bgc_class = (
                get_ci("bigscape_class") or
                get_ci("product_class") or
                get_ci("category") or
                get_ci("class") or
                "Unknown"
            )

            product = (
                get_ci("product") or
                get_ci("Product Prediction") or
                get_ci("description") or
                ""
            )

            organism = get_ci("organism")

            bgc_info[bgc_name] = {
                "class": bgc_class,
                "product": product,
                "organism": organism,
                "description": get_ci("description"),
            }

    print(f"BGC annotations loaded: {len(bgc_info)}")

    # ================================================================
    # Step 4: Parse network edges
    # ================================================================
    edges = []

    # Prefer the 'mix' network if available
    chosen_network = network_files[0]
    for nf in network_files:
        if "mix" in os.path.basename(nf).lower():
            chosen_network = nf
            break

    print(f"Using network file: {chosen_network}")

    with open(chosen_network) as f:
        reader = csv.DictReader(f, delimiter="\\t")
        net_headers = reader.fieldnames
        print(f"Network columns: {net_headers}")

        for row in reader:
            bgc1 = row.get(net_headers[0], "")
            bgc2 = row.get(net_headers[1], "")

            dist = 0.5
            for col in ["Raw distance", "raw_distance", "distance", "Distance"]:
                if col in row:
                    try:
                        dist = float(row[col])
                    except (ValueError, TypeError):
                        pass
                    break

            if bgc1 and bgc2:
                edges.append((bgc1, bgc2, dist))

    print(f"Network edges loaded: {len(edges)}")

    if not edges:
        print("WARNING: No network edges found. No plot generated.")
        sys.exit(0)

    # ================================================================
    # Step 5: Parse clustering (GCF assignments)
    # ================================================================
    bgc_gcf = {}
    if clustering_files:
        with open(clustering_files[0]) as f:
            reader = csv.reader(f, delimiter="\\t")
            header = next(reader, None)
            for row in reader:
                if len(row) >= 2:
                    bgc_gcf[row[0]] = row[1]
        print(f"GCF assignments loaded: {len(bgc_gcf)}")

    # ================================================================
    # Step 6: Map BGC names to sample_id and genus using SQLite DB
    # ================================================================

    # Build record_name -> sample_id from BiG-SCAPE database
    record_to_sample = {}
    db_files = glob.glob(os.path.join(results_dir, "**", "*.db"), recursive=True)
    for db_path in db_files:
        try:
            conn = sqlite3.connect(db_path)
            cur = conn.cursor()
            # Join bgc_record with gbk to get path (which contains sample_id)
            cur.execute("SELECT br.id, g.path, br.record_number, g.organism FROM bgc_record br JOIN gbk g ON br.gbk_id = g.id WHERE br.record_type = 'region'")
            for rec_id, gbk_path, rec_num, organism in cur.fetchall():
                # Path format: group_id/SAMPLE_ID/contig.region.gbk or MiBIG path
                parts = gbk_path.replace("\\\\", "/").split("/")
                gbk_basename = parts[-1].replace(".gbk", "")
                record_name = f"{gbk_basename}.gbk_region_{rec_num}"

                if "mibig" in gbk_path.lower() or gbk_path.startswith("BGC"):
                    record_to_sample[record_name] = "__MIBIG__"
                elif len(parts) >= 2:
                    # sample_id is the second-to-last directory
                    sample_id = parts[-2]
                    record_to_sample[record_name] = sample_id
            conn.close()
            print(f"DB mapping loaded: {len(record_to_sample)} records from {db_path}")
        except Exception as e:
            print(f"Warning: Could not read DB {db_path}: {e}")

    # Fallback: regex for PAL\\d+ in name
    def extract_sample_id(bgc_name):
        if bgc_name in record_to_sample:
            return record_to_sample[bgc_name]
        m = re.match(r"(PAL\\d+)", bgc_name)
        if m:
            return m.group(1)
        parts = bgc_name.split("_")
        for i in range(1, len(parts) + 1):
            candidate = "_".join(parts[:i])
            if candidate in sample_genus:
                return candidate
        return None

    bgc_genus = {}
    bgc_sample = {}
    all_bgc_names = set(list(bgc_info.keys()) + [e[0] for e in edges] + [e[1] for e in edges])
    for bgc_name in all_bgc_names:
        sid = extract_sample_id(bgc_name)
        if sid == "__MIBIG__":
            bgc_genus[bgc_name] = "MiBIG"
        elif sid and sid in sample_genus:
            bgc_sample[bgc_name] = sid
            bgc_genus[bgc_name] = sample_genus.get(sid, "Unknown")
        else:
            if bgc_name.startswith("BGC") or "mibig" in bgc_name.lower():
                bgc_genus[bgc_name] = "MiBIG"
            else:
                bgc_genus[bgc_name] = "Unknown"

    # Report mapping stats
    n_mapped = sum(1 for g in bgc_genus.values() if g not in ("Unknown", "MiBIG"))
    n_mibig = sum(1 for g in bgc_genus.values() if g == "MiBIG")
    n_unknown = sum(1 for g in bgc_genus.values() if g == "Unknown")
    print(f"BGC genus mapping: {n_mapped} mapped, {n_mibig} MiBIG, {n_unknown} unknown")

    # ================================================================
    # Step 7: Build network graph grouped by BGC class
    # ================================================================
    G = nx.Graph()

    all_bgcs = set()
    for bgc1, bgc2, dist in edges:
        all_bgcs.add(bgc1)
        all_bgcs.add(bgc2)

    for bgc in all_bgcs:
        info = bgc_info.get(bgc, {})
        G.add_node(bgc,
                   bgc_class=info.get("class", "Unknown"),
                   product=info.get("product", ""),
                   genus=bgc_genus.get(bgc, "Unknown"))

    for bgc1, bgc2, dist in edges:
        G.add_edge(bgc1, bgc2, weight=max(0.01, 1.0 - dist))

    class_nodes = defaultdict(list)
    for node, data in G.nodes(data=True):
        cls = data.get("bgc_class", "Unknown") or "Unknown"
        class_nodes[cls].append(node)

    sorted_classes = sorted(class_nodes.keys(), key=lambda c: len(class_nodes[c]), reverse=True)

    classes_with_content = []
    for cls in sorted_classes:
        nodes = class_nodes[cls]
        subgraph = G.subgraph(nodes)
        if subgraph.number_of_edges() > 0 or len(nodes) >= 2:
            classes_with_content.append(cls)

    if not classes_with_content:
        print("WARNING: No BGC classes with edges. No plot generated.")
        sys.exit(0)

    print(f"\\nBGC classes to plot: {len(classes_with_content)}")
    for cls in classes_with_content:
        print(f"  {cls}: {len(class_nodes[cls])} BGCs")

    # ================================================================
    # Step 8: Genus color palette
    # ================================================================
    genus_colors_map = {
        "Lacticaseibacillus":      "#E74C3C",
        "Companilactobacillus":    "#FF6B35",
        "Lactobacillus":           "#F39C12",
        "Ligilactobacillus":       "#2ECC71",
        "Lactiplantibacillus":     "#27AE60",
        "Lentilactobacillus":      "#1ABC9C",
        "Levilactobacillus":       "#3498DB",
        "Limosilactobacillus":     "#9B59B6",
        "Latilactobacillus":       "#E91E63",
        "Weissella":               "#FFC107",
        "Leuconostoc":             "#00BCD4",
        "Pediococcus":             "#FF9800",
        "Loigolactobacillus":      "#8BC34A",
        "Liquorilactobacillus":    "#795548",
        "Schleiferilactobacillus": "#607D8B",
        "MiBIG":                   "#BDBDBD",
        "Unknown":                 "#9E9E9E",
    }

    all_genera = sorted(set(bgc_genus.values()))
    cmap = plt.cm.get_cmap("tab20", max(20, len(all_genera)))
    ci = 0
    for g in all_genera:
        if g not in genus_colors_map:
            genus_colors_map[g] = matplotlib.colors.rgb2hex(cmap(ci % 20))
            ci += 1

    # ================================================================
    # Step 9: Multi-panel figure
    # ================================================================
    n_panels = len(classes_with_content)
    n_cols = min(4, n_panels)
    n_rows = math.ceil(n_panels / n_cols)

    fig_width = max(12, n_cols * 5)
    fig_height = max(8, n_rows * 5)
    fig, all_axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))

    # Flatten axes to 1D list for easy indexing
    if n_panels == 1:
        ax_list = [all_axes]
    elif n_rows == 1 or n_cols == 1:
        ax_list = list(all_axes.flat) if hasattr(all_axes, 'flat') else [all_axes]
    else:
        ax_list = list(all_axes.flat)

    for idx, cls in enumerate(classes_with_content):
        ax = ax_list[idx]
        nodes = class_nodes[cls]
        subgraph = G.subgraph(nodes).copy()

        if subgraph.number_of_nodes() > 0:
            try:
                k_val = 2.0 / math.sqrt(max(1, subgraph.number_of_nodes()))
                pos = nx.spring_layout(subgraph, k=k_val, iterations=100, seed=42)
            except Exception:
                pos = nx.random_layout(subgraph, seed=42)
        else:
            pos = {}

        # Edges
        nx.draw_networkx_edges(subgraph, pos, ax=ax, alpha=0.15, edge_color="#CCCCCC", width=0.5)

        # Nodes
        node_list = list(subgraph.nodes())
        node_colors = [genus_colors_map.get(G.nodes[n].get("genus", "Unknown"), "#9E9E9E") for n in node_list]
        node_sizes = [30 if bgc_genus.get(n, "") == "MiBIG" else 60 for n in node_list]

        nx.draw_networkx_nodes(subgraph, pos, ax=ax, nodelist=node_list,
                               node_color=node_colors, node_size=node_sizes,
                               alpha=0.85, edgecolors="white", linewidths=0.3)

        # Labels for MiBIG/known products
        labels = {}
        for n in node_list:
            product = G.nodes[n].get("product", "")
            if product and bgc_genus.get(n, "") == "MiBIG":
                label = product.split(",")[0].strip()
                if len(label) > 25:
                    label = label[:22] + "..."
                labels[n] = label

        if labels:
            nx.draw_networkx_labels(subgraph, pos, labels=labels, ax=ax,
                                    font_size=4, font_color="#C0392B", alpha=0.8)

        ax.set_title(cls, fontsize=10, fontweight="bold", pad=4)
        ax.axis("off")

    # Hide unused panels
    for idx in range(n_panels, len(ax_list)):
        ax_list[idx].axis("off")

    # ================================================================
    # Step 10: Legend
    # ================================================================
    present_genera = sorted(set(bgc_genus[n] for n in G.nodes() if n in bgc_genus))
    legend_patches = [mpatches.Patch(color=genus_colors_map.get(g, "#9E9E9E"), label=g)
                      for g in present_genera]

    fig.legend(handles=legend_patches, loc="lower left",
               ncol=min(4, len(legend_patches)),
               fontsize=7, frameon=True, title="Genus", title_fontsize=8,
               bbox_to_anchor=(0.02, 0.01))

    display_name = group_id.replace("genus_", "").replace("species_", "").replace("_", " ")
    fig.suptitle(f"BGC Similarity Network -- {display_name}",
                 fontsize=14, fontweight="bold", y=0.98)

    plt.tight_layout(rect=[0, 0.06, 1, 0.96])

    # ================================================================
    # Step 11: Save
    # ================================================================
    cutoff_match = re.search(r"c?(0\\.\\d+)", os.path.basename(chosen_network))
    cutoff_str = cutoff_match.group(1) if cutoff_match else "0.30"

    png_path = f"bgc_network_{cutoff_str}.png"
    pdf_path = f"bgc_network_{cutoff_str}.pdf"

    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"\\nSaved: {png_path}")
    print(f"Saved: {pdf_path}")
    print(f"Total BGCs plotted: {G.number_of_nodes()}")
    print(f"Total edges: {G.number_of_edges()}")
    print("Network plot completed!")
    """
}
