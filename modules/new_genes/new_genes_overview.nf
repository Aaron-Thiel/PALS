/*
 * New Genes Overview
 *
 * Aggregates per-sample new gene reports from GENOME_VISUALIZATION
 * and groups them by family, genus, and species using taxonomy from
 * pipeline_summary.tsv.
 *
 * When PANTA gene_presence_absence.csv and EggNOG annotations are
 * available, maps new genes to COG categories and generates:
 *   - COG distribution plots (by classification type)
 *   - COG comparison: new genes vs full pangenome
 *   - Per-species and per-genus COG subplot grids
 *   - Enriched detail table with panta_group and cog_category columns
 */

process NEW_GENES_OVERVIEW {
    tag "new_genes_overview"
    label 'process_low'

    publishDir "${params.outdir}/genomeviz/overview", mode: 'copy'

    conda "conda-forge::pandas conda-forge::matplotlib"
    container null

    input:
    path(genomeviz_dirs)
    path(pipeline_summary)
    path(gene_presence_absence_csv)
    path(eggnog_annotations)
    path(rtab_file)

    output:
    path("new_genes_overview.tsv"), emit: summary
    path("new_genes_all_samples.tsv"), emit: details
    path("summary_by_*.tsv"), emit: taxonomy_summaries
    path("cog_*.png"), emit: cog_plots_png, optional: true
    path("cog_*.pdf"), emit: cog_plots_pdf, optional: true
    path("cog_summary.tsv"), emit: cog_summary, optional: true
    path("species"), emit: species_dir, optional: true
    path("genus"), emit: genus_dir, optional: true

    script:
    """
    #!/usr/bin/env python3

    import csv
    import os
    import math
    from collections import defaultdict, Counter

    print("=" * 60)
    print("New Genes Overview")
    print("=" * 60)

    # ================================================================
    # COG constants
    # ================================================================
    COG_ORDER = ['-', 'J', 'A', 'K', 'L', 'B', 'D', 'Y', 'V', 'T', 'M', 'N',
                 'Z', 'W', 'U', 'O', 'C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q', 'R', 'S']

    COG_NAMES = {
        '-': 'Not assigned', 'J': 'Translation', 'A': 'RNA processing',
        'K': 'Transcription', 'L': 'Replication/repair', 'B': 'Chromatin',
        'D': 'Cell cycle', 'Y': 'Nuclear structure', 'V': 'Defense',
        'T': 'Signal transduction', 'M': 'Cell wall/membrane',
        'N': 'Cell motility', 'Z': 'Cytoskeleton', 'W': 'Extracellular',
        'U': 'Intracellular trafficking', 'O': 'Post-translational modification',
        'C': 'Energy production', 'G': 'Carbohydrate metabolism',
        'E': 'Amino acid metabolism', 'F': 'Nucleotide metabolism',
        'H': 'Coenzyme metabolism', 'I': 'Lipid metabolism',
        'P': 'Inorganic ion transport', 'Q': 'Secondary metabolites',
        'R': 'General function', 'S': 'Function unknown',
    }

    CLASSIFICATION_COLORS = {
        'Other Source': '#E74C3C',
        'Contextual Re-prediction': '#3498DB',
        'Combined Contigs': '#2ECC71',
    }

    # ================================================================
    # Step 1: Read taxonomy from pipeline_summary.tsv
    # ================================================================
    sample_taxonomy = {}
    with open("${pipeline_summary}") as f:
        reader = csv.DictReader(f, delimiter="\\t")
        for row in reader:
            sid = row['sample_id']
            family = row.get('family', '')
            genus = row.get('genus', '')
            species = row.get('species', '')
            sample_taxonomy[sid] = {
                'family': family if family and family != 'NA' else 'Unknown',
                'genus': genus if genus and genus != 'NA' else 'Unknown',
                'species': species if species and species != 'NA' else 'Unknown',
            }

    print(f"Samples with taxonomy: {len(sample_taxonomy)}")

    # ================================================================
    # Step 2: Collect new gene reports from each sample
    # ================================================================
    all_genes = []       # list of dicts for detail table
    sample_stats = []    # per-sample summary

    for entry in sorted(os.listdir('.')):
        if not os.path.isdir(entry):
            continue

        report_path = os.path.join(entry, 'assembly_comparison', 'new_genes', 'gene_comparison_report.csv')
        if not os.path.exists(report_path):
            continue

        sample_id = entry
        tax = sample_taxonomy.get(sample_id, {'family': 'Unknown', 'genus': 'Unknown', 'species': 'Unknown'})

        # Parse gene_comparison_report.csv
        classifications = defaultdict(int)
        genes = []
        with open(report_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                classification = row.get('Classification', 'Unknown')
                classifications[classification] += 1
                genes.append(row)

                all_genes.append({
                    'sample_id': sample_id,
                    'family': tax['family'],
                    'genus': tax['genus'],
                    'species': tax['species'],
                    'gene_id': row.get('GeneID', ''),
                    'product': row.get('Product', ''),
                    'scaffold': row.get('Scaffold', ''),
                    'classification': classification,
                    'joined': row.get('Joined', ''),
                    'panta_group': '',
                    'cog_category': '',
                })

        total = len(genes)
        n_other = classifications.get('Other Source', 0)
        n_contextual = classifications.get('Contextual Re-prediction', 0)
        n_novel = total - n_other - n_contextual

        sample_stats.append({
            'sample_id': sample_id,
            'family': tax['family'],
            'genus': tax['genus'],
            'species': tax['species'],
            'total_new_genes': total,
            'other_source': n_other,
            'contextual_reprediction': n_contextual,
            'novel': n_novel,
        })

        print(f"  {sample_id} ({tax['genus']} {tax['species']}): {total} new genes")

    print(f"\\nTotal samples with new genes data: {len(sample_stats)}")
    print(f"Total new genes across all samples: {len(all_genes)}")

    # ================================================================
    # Step 3: Write per-sample overview
    # ================================================================
    with open("new_genes_overview.tsv", 'w') as f:
        f.write("sample_id\\tfamily\\tgenus\\tspecies\\ttotal_new_genes\\tother_source\\tcontextual_reprediction\\tnovel\\n")
        for s in sorted(sample_stats, key=lambda x: x['sample_id']):
            f.write(f"{s['sample_id']}\\t{s['family']}\\t{s['genus']}\\t{s['species']}\\t"
                    f"{s['total_new_genes']}\\t{s['other_source']}\\t{s['contextual_reprediction']}\\t{s['novel']}\\n")

    # ================================================================
    # Step 5: Write taxonomy-level summaries
    # ================================================================
    for level in ['family', 'genus', 'species']:
        groups = defaultdict(lambda: {'n_samples': 0, 'total': 0, 'other_source': 0,
                                       'contextual_reprediction': 0, 'novel': 0})
        for s in sample_stats:
            key = s[level]
            if key == 'Unknown':
                continue
            groups[key]['n_samples'] += 1
            groups[key]['total'] += s['total_new_genes']
            groups[key]['other_source'] += s['other_source']
            groups[key]['contextual_reprediction'] += s['contextual_reprediction']
            groups[key]['novel'] += s['novel']

        with open(f"summary_by_{level}.tsv", 'w') as f:
            f.write(f"{level}\\tn_samples\\ttotal_new_genes\\tmean_new_genes\\tother_source\\tcontextual_reprediction\\tnovel\\n")
            for name in sorted(groups.keys()):
                g = groups[name]
                mean = g['total'] / g['n_samples'] if g['n_samples'] > 0 else 0
                f.write(f"{name}\\t{g['n_samples']}\\t{g['total']}\\t{mean:.1f}\\t"
                        f"{g['other_source']}\\t{g['contextual_reprediction']}\\t{g['novel']}\\n")

        print(f"\\n{level.capitalize()}-level summary: {len(groups)} groups")
        for name in sorted(groups.keys()):
            g = groups[name]
            print(f"  {name}: {g['n_samples']} samples, {g['total']} new genes (mean {g['total']/g['n_samples']:.1f})")

    # ================================================================
    # Step 6: COG Analysis (requires EggNOG + PANTA)
    # ================================================================
    eggnog_file = "${eggnog_annotations}"
    gpa_csv_file = "${gene_presence_absence_csv}"
    rtab_input = "${rtab_file}"

    has_eggnog = os.path.exists(eggnog_file) and os.path.basename(eggnog_file) != "NO_EGGNOG" and os.path.getsize(eggnog_file) > 0
    has_gpa = os.path.exists(gpa_csv_file) and os.path.getsize(gpa_csv_file) > 0
    has_rtab = os.path.exists(rtab_input) and os.path.getsize(rtab_input) > 0

    if not has_eggnog or not has_gpa or not has_rtab or len(all_genes) == 0:
        if not has_eggnog:
            print("\\nSkipping COG analysis: EggNOG annotations not available")
        elif not has_gpa:
            print("\\nSkipping COG analysis: gene_presence_absence.csv not available")
        elif not has_rtab:
            print("\\nSkipping COG analysis: Rtab not available")
        else:
            print("\\nSkipping COG analysis: no new genes found")

        # Write detail table without COG columns
        with open("new_genes_all_samples.tsv", 'w') as f:
            f.write("sample_id\\tfamily\\tgenus\\tspecies\\tgene_id\\tproduct\\tscaffold\\tclassification\\tjoined\\n")
            for g in sorted(all_genes, key=lambda x: (x['sample_id'], x['gene_id'])):
                f.write(f"{g['sample_id']}\\t{g['family']}\\t{g['genus']}\\t{g['species']}\\t"
                        f"{g['gene_id']}\\t{g['product']}\\t{g['scaffold']}\\t{g['classification']}\\t{g['joined']}\\n")

        print("\\nNew genes overview completed!")
        import sys
        sys.exit(0)

    print("\\n" + "=" * 60)
    print("COG Category Analysis")
    print("=" * 60)

    # ----------------------------------------------------------------
    # Step 6a: Load EggNOG annotations (group -> COG category)
    # ----------------------------------------------------------------
    group_to_cog = {}
    with open(eggnog_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\\t")
            if len(parts) >= 7:
                query = parts[0]
                cog_cat = parts[6] if parts[6] and parts[6] != "-" else "-"
                group_to_cog[query] = cog_cat

    print(f"EggNOG annotations loaded: {len(group_to_cog)} gene groups")

    # ----------------------------------------------------------------
    # Step 7: Map new genes to PANTA groups via gene_presence_absence.csv
    # ----------------------------------------------------------------
    # Build target lookup set: {sample_id}-{scaffold}-{gene_id}
    target_keys = set()
    for g in all_genes:
        key = f"{g['sample_id']}-{g['scaffold']}-{g['gene_id']}"
        target_keys.add(key)

    print(f"\\nMapping {len(target_keys)} new genes to PANTA groups...")

    gene_to_group = {}  # (sample_id, gene_id) -> group_name

    with open(gpa_csv_file) as f:
        reader = csv.reader(f)
        header = next(reader)

        # Find internal sample columns (PAL*)
        pal_cols = {}
        for i, h in enumerate(header):
            if h.startswith("PAL"):
                pal_cols[i] = h

        print(f"  Internal sample columns: {len(pal_cols)}")

        # Stream rows
        found = 0
        for row_num, row in enumerate(reader):
            group_name = row[0]

            for col_idx in pal_cols:
                cell = row[col_idx].strip() if col_idx < len(row) else ""
                if cell and cell in target_keys:
                    parts = cell.split("-", 2)
                    if len(parts) == 3:
                        sid, _scaffold, gid = parts
                        gene_to_group[(sid, gid)] = group_name
                        found += 1

            if (row_num + 1) % 100000 == 0:
                print(f"  Scanned {row_num + 1} rows, found {found} matches...")

            if found >= len(target_keys):
                print(f"  All targets found at row {row_num + 1}, stopping early")
                break

    print(f"  Mapped {found}/{len(target_keys)} new genes to PANTA groups")

    # ----------------------------------------------------------------
    # Step 8: Enrich new genes with COG categories
    # ----------------------------------------------------------------
    for g in all_genes:
        group = gene_to_group.get((g['sample_id'], g['gene_id']), '')
        if group:
            g['panta_group'] = group
            g['cog_category'] = group_to_cog.get(group, '-')
        else:
            g['panta_group'] = ''
            g['cog_category'] = ''

    n_mapped = sum(1 for g in all_genes if g['panta_group'])
    n_with_cog = sum(1 for g in all_genes if g['cog_category'] and g['cog_category'] != '-')
    print(f"\\nNew genes mapped to PANTA groups: {n_mapped}/{len(all_genes)}")
    print(f"New genes with COG annotation (non-'-'): {n_with_cog}/{len(all_genes)}")

    # ----------------------------------------------------------------
    # Step 9: Write enriched detail table
    # ----------------------------------------------------------------
    with open("new_genes_all_samples.tsv", 'w') as f:
        f.write("sample_id\\tfamily\\tgenus\\tspecies\\tgene_id\\tproduct\\tscaffold\\t"
                "classification\\tjoined\\tpanta_group\\tcog_category\\n")
        for g in sorted(all_genes, key=lambda x: (x['sample_id'], x['gene_id'])):
            f.write(f"{g['sample_id']}\\t{g['family']}\\t{g['genus']}\\t{g['species']}\\t"
                    f"{g['gene_id']}\\t{g['product']}\\t{g['scaffold']}\\t{g['classification']}\\t"
                    f"{g['joined']}\\t{g['panta_group']}\\t{g['cog_category']}\\n")

    # ----------------------------------------------------------------
    # Step 10: Compute COG distributions
    # ----------------------------------------------------------------

    def normalize_classification(cls):
        if cls.startswith('Combined Contigs'):
            return 'Combined Contigs'
        return cls

    # New genes: per-classification COG counts
    cls_cog_counts = {}  # classification -> Counter
    new_gene_cog_total = Counter()
    for g in all_genes:
        cog = g.get('cog_category', '-') or '-'
        cls = normalize_classification(g['classification'])
        if cls not in cls_cog_counts:
            cls_cog_counts[cls] = Counter()
        for letter in cog:
            if letter in COG_ORDER:
                cls_cog_counts[cls][letter] += 1
                new_gene_cog_total[letter] += 1
            elif cog == '-' or cog == '':
                cls_cog_counts[cls]['-'] += 1
                new_gene_cog_total['-'] += 1
                break

    # Pangenome-wide COG baseline from Rtab gene names
    pangenome_cog_counts = Counter()
    pangenome_gene_count = 0
    with open(rtab_input) as f:
        next(f)  # skip header
        for line in f:
            gene_name = line.split("\\t", 1)[0]
            pangenome_gene_count += 1
            cog = group_to_cog.get(gene_name, '-')
            for letter in cog:
                if letter in COG_ORDER:
                    pangenome_cog_counts[letter] += 1
                elif cog == '-':
                    pangenome_cog_counts['-'] += 1
                    break

    print(f"\\nPangenome genes (from Rtab): {pangenome_gene_count}")

    # Per-species and per-genus COG counts
    species_cls_cog = {}  # species -> classification -> Counter
    genus_cls_cog = {}    # genus -> classification -> Counter
    for g in all_genes:
        cog = g.get('cog_category', '-') or '-'
        cls = normalize_classification(g['classification'])

        for level_name, level_dict in [('species', species_cls_cog), ('genus', genus_cls_cog)]:
            key = g[level_name]
            if key == 'Unknown':
                continue
            if key not in level_dict:
                level_dict[key] = {}
            if cls not in level_dict[key]:
                level_dict[key][cls] = Counter()
            for letter in cog:
                if letter in COG_ORDER:
                    level_dict[key][cls][letter] += 1
                elif cog == '-' or cog == '':
                    level_dict[key][cls]['-'] += 1
                    break

    # ----------------------------------------------------------------
    # Step 11: Write global COG summary table
    # ----------------------------------------------------------------
    total_new = sum(new_gene_cog_total.values()) or 1
    total_pan = sum(pangenome_cog_counts.values()) or 1

    with open("cog_summary.tsv", 'w') as f:
        f.write("cog_category\\tcog_name\\tnew_genes_count\\tnew_genes_pct\\tpangenome_count\\tpangenome_pct\\n")
        for cat in COG_ORDER:
            ng_count = new_gene_cog_total.get(cat, 0)
            ng_pct = 100.0 * ng_count / total_new
            pg_count = pangenome_cog_counts.get(cat, 0)
            pg_pct = 100.0 * pg_count / total_pan
            f.write(f"{cat}\\t{COG_NAMES.get(cat, '')}\\t{ng_count}\\t{ng_pct:.1f}\\t{pg_count}\\t{pg_pct:.1f}\\n")

    print("COG summary table written: cog_summary.tsv")

    # ----------------------------------------------------------------
    # Step 11a: Helper to write per-group COG summary TSV
    # ----------------------------------------------------------------
    def write_cog_summary(filepath, cls_cog_dict, group_name):
        total = Counter()
        for cls_counter in cls_cog_dict.values():
            total += cls_counter
        total_count = sum(total.values()) or 1
        with open(filepath, 'w') as f:
            f.write(f"# {group_name}\\n")
            f.write("cog_category\\tcog_name\\tcount\\tpct\\tother_source\\tcontextual_reprediction\\tcombined_contigs\\n")
            for cat in COG_ORDER:
                count = total.get(cat, 0)
                pct = 100.0 * count / total_count
                os_count = cls_cog_dict.get('Other Source', Counter()).get(cat, 0)
                cr_count = cls_cog_dict.get('Contextual Re-prediction', Counter()).get(cat, 0)
                cc_count = cls_cog_dict.get('Combined Contigs', Counter()).get(cat, 0)
                f.write(f"{cat}\\t{COG_NAMES.get(cat, '')}\\t{count}\\t{pct:.1f}\\t{os_count}\\t{cr_count}\\t{cc_count}\\n")

    # ----------------------------------------------------------------
    # Step 11b: Write per-species and per-genus COG summaries
    # ----------------------------------------------------------------
    os.makedirs("species", exist_ok=True)
    os.makedirs("genus", exist_ok=True)

    for species, cls_cog in species_cls_cog.items():
        safe_name = species.replace(' ', '_').replace('/', '_')
        write_cog_summary(f"species/cog_summary_{safe_name}.tsv", cls_cog, species.replace('_', ' '))
    print(f"Per-species COG summaries written: {len(species_cls_cog)} files")

    for genus, cls_cog in genus_cls_cog.items():
        safe_name = genus.replace(' ', '_').replace('/', '_')
        write_cog_summary(f"genus/cog_summary_{safe_name}.tsv", cls_cog, genus)
    print(f"Per-genus COG summaries written: {len(genus_cls_cog)} files")

    # ----------------------------------------------------------------
    # Step 12: Generate plots
    # ----------------------------------------------------------------
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import numpy as np

    legend_elements = [Patch(facecolor=c, label=l) for l, c in CLASSIFICATION_COLORS.items()]

    # Helper: draw a stacked COG bar chart on an axis
    def plot_stacked_cog(ax, cls_cog, title, fontsize_title=14):
        x = np.arange(len(COG_ORDER))
        width = 0.7
        bottom = np.zeros(len(COG_ORDER))
        for cls_name in ['Other Source', 'Contextual Re-prediction', 'Combined Contigs']:
            if cls_name not in cls_cog:
                continue
            counts = np.array([cls_cog[cls_name].get(cat, 0) for cat in COG_ORDER])
            ax.bar(x, counts, width, bottom=bottom,
                   color=CLASSIFICATION_COLORS.get(cls_name, '#999999'),
                   edgecolor='white', linewidth=0.3)
            bottom += counts
        ax.set_xticks(x)
        ax.set_xticklabels([f'[{c}]' for c in COG_ORDER], fontsize=9, fontweight='bold')
        ax.set_xlabel('COG Category', fontsize=12)
        ax.set_ylabel('Gene Count', fontsize=12)
        ax.set_title(title, fontsize=fontsize_title, fontweight='bold')
        ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    # --- Plot 1: COG distribution by classification type ---
    print("\\nGenerating COG distribution plots...")
    fig, ax = plt.subplots(figsize=(14, 7))
    plot_stacked_cog(ax, cls_cog_counts, 'COG Distribution of New Genes by Classification Type')
    plt.tight_layout()
    fig.savefig('cog_new_genes_by_classification.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig('cog_new_genes_by_classification.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print("  cog_new_genes_by_classification.png")

    # --- Plot 2: COG comparison (new genes vs pangenome) ---
    fig, ax = plt.subplots(figsize=(14, 7))

    pangenome_props = np.array([pangenome_cog_counts.get(cat, 0) / total_pan for cat in COG_ORDER])
    new_props = np.array([new_gene_cog_total.get(cat, 0) / total_new for cat in COG_ORDER])

    x = np.arange(len(COG_ORDER))
    width = 0.35
    ax.bar(x - width/2, pangenome_props, width, label=f'All Pangenome ({pangenome_gene_count:,} genes)',
           color='#95A5A6', edgecolor='white', linewidth=0.5)
    ax.bar(x + width/2, new_props, width, label=f'New Genes ({len(all_genes):,} genes)',
           color='#E74C3C', edgecolor='white', linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels([f'[{c}]' for c in COG_ORDER], fontsize=9, fontweight='bold')
    ax.set_xlabel('COG Category', fontsize=12)
    ax.set_ylabel('Proportion of Genes', fontsize=12)
    ax.set_title('COG Distribution: New Genes vs All Pangenome Genes', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    plt.tight_layout()
    fig.savefig('cog_new_vs_pangenome.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig('cog_new_vs_pangenome.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print("  cog_new_vs_pangenome.png")

    # --- Plot 3: Individual per-species COG plots ---
    species_list = sorted(species_cls_cog.keys())
    if species_list:
        for species in species_list:
            display_name = species.replace('_', ' ')
            fig, ax = plt.subplots(figsize=(14, 7))
            plot_stacked_cog(ax, species_cls_cog[species],
                             f'COG Distribution of New Genes \\u2014 {display_name}')
            plt.tight_layout()
            safe_name = species.replace(' ', '_').replace('/', '_')
            fig.savefig(f'species/cog_{safe_name}.png', dpi=300, bbox_inches='tight', facecolor='white')
            fig.savefig(f'species/cog_{safe_name}.pdf', bbox_inches='tight', facecolor='white')
            plt.close(fig)
        print(f"  Generated {len(species_list)} individual species COG plots in species/")

    # --- Plot 4: Individual per-genus COG plots ---
    genus_list = sorted(genus_cls_cog.keys())
    if genus_list:
        for genus in genus_list:
            fig, ax = plt.subplots(figsize=(14, 7))
            plot_stacked_cog(ax, genus_cls_cog[genus],
                             f'COG Distribution of New Genes \\u2014 {genus}')
            plt.tight_layout()
            safe_name = genus.replace(' ', '_').replace('/', '_')
            fig.savefig(f'genus/cog_{safe_name}.png', dpi=300, bbox_inches='tight', facecolor='white')
            fig.savefig(f'genus/cog_{safe_name}.pdf', bbox_inches='tight', facecolor='white')
            plt.close(fig)
        print(f"  Generated {len(genus_list)} individual genus COG plots in genus/")

    print("\\nNew genes overview completed!")
    """
}
