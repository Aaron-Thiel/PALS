/*
 * Filter Rtab to HQ samples only
 *
 * Removes internal samples whose post-scaffolding completeness
 * (pasa_completeness) falls below a threshold. External samples
 * are always kept.
 */

process FILTER_RTAB_HQ {
    tag "$cohort_id"

    container 'python:3.11'

    input:
    tuple val(cohort_id), path(rtab_file)
    path(pipeline_summary)
    val(min_completeness)

    output:
    tuple val(cohort_id), path("hq_filtered.Rtab"), emit: rtab

    script:
    """
    #!/usr/bin/env python3

    import csv, os, sys

    print("=" * 60)
    print("Filter Rtab to HQ samples")
    print("=" * 60)
    print(f"Rtab file: ${rtab_file}")
    print(f"Min pasa_completeness: ${min_completeness}%")
    print("")

    min_comp = float(${min_completeness})

    # ── Read pipeline_summary to get pasa_completeness per internal sample ──
    exclude = set()
    internal_total = 0

    with open("${pipeline_summary}") as f:
        reader = csv.DictReader(f, delimiter="\\t")
        for row in reader:
            sid = row["sample_id"]
            pasa_comp = row.get("pasa_completeness", "NA")
            if pasa_comp in ("NA", ""):
                # No scaffold → sample won't be in Rtab anyway
                continue
            internal_total += 1
            if float(pasa_comp) < min_comp:
                exclude.add(sid)
                print(f"  EXCLUDE {sid}: pasa_completeness={pasa_comp}%")

    print(f"\\nInternal samples scanned: {internal_total}")
    print(f"Excluded (below {min_comp}%): {len(exclude)}")

    # ── Filter Rtab columns ──
    with open("${rtab_file}") as fin:
        header = fin.readline().strip().split("\\t")

    gene_col = header[0]
    sample_names = header[1:]

    # Keep columns that are NOT in the exclude set
    keep_indices = [i for i, s in enumerate(sample_names) if s not in exclude]
    keep_names = [sample_names[i] for i in keep_indices]

    removed = len(sample_names) - len(keep_names)
    print(f"Samples before filter: {len(sample_names)}")
    print(f"Samples after filter:  {len(keep_names)} (removed {removed})")

    # Write filtered Rtab
    row_count = 0
    kept_genes = 0

    with open("${rtab_file}") as fin, open("hq_filtered.Rtab", "w") as fout:
        next(fin)  # skip header
        fout.write(gene_col + "\\t" + "\\t".join(keep_names) + "\\n")

        for line in fin:
            parts = line.strip().split("\\t")
            gene_name = parts[0]
            values = parts[1:]
            subset = [values[i] for i in keep_indices]

            # Only keep genes present in at least one remaining sample
            if any(v != "0" for v in subset):
                fout.write(gene_name + "\\t" + "\\t".join(subset) + "\\n")
                kept_genes += 1
            row_count += 1

    print(f"Genes before filter: {row_count}")
    print(f"Genes after filter:  {kept_genes}")
    print("\\nHQ filtering complete!")
    """
}
