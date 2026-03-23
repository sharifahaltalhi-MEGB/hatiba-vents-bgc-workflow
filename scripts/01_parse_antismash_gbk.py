import os
import csv
from Bio import SeqIO

root_dir = "."  # run inside antismash8_1136MAGs_HTV (output)
output_csv = "bgc_summary_unique_ids_v8.csv"

def get_genome_name(subdir: str, root: str) -> str:
    rel = os.path.relpath(subdir, root)
    parts = rel.split(os.sep)
    return parts[0] if parts and parts[0] not in [".", ""] else os.path.basename(subdir)

def strip_gbk(filename: str) -> str:
    return filename[:-4] if filename.endswith(".gbk") else filename

rows = []
seen_ids = set()
failed = 0

for subdir, _, files in os.walk(root_dir):
    for file in files:
        if not (file.endswith(".gbk") and ".region" in file):
            continue

        filepath = os.path.join(subdir, file)

        genome_name = get_genome_name(subdir, root_dir)
        bgc_name = strip_gbk(file)  # e.g. S2_k141_2110384.region001
        redsea_bgc_id_base = f"{genome_name}/{bgc_name}"

        try:
            with open(filepath, "r") as handle:
                rec_index = 0
                for record in SeqIO.parse(handle, "genbank"):
                    rec_index += 1

                    cluster_type = "N/A"
                    start = "N/A"
                    end = "N/A"
                    strand = "N/A"

                    chosen_feature = None
                    for feature in record.features:
                        if feature.type in ["protocluster", "cluster"]:
                            chosen_feature = feature
                            break

                    if chosen_feature is not None:
                        cluster_type = chosen_feature.qualifiers.get("product", ["N/A"])[0]
                        start = int(chosen_feature.location.start)
                        end = int(chosen_feature.location.end)
                        strand = chosen_feature.location.strand if chosen_feature.location.strand is not None else "N/A"

                    redsea_bgc_id = redsea_bgc_id_base
                    if rec_index > 1:
                        redsea_bgc_id = f"{redsea_bgc_id_base}__rec{rec_index}"

                    if redsea_bgc_id in seen_ids:
                        bump = 2
                        new_id = f"{redsea_bgc_id}__dup{bump}"
                        while new_id in seen_ids:
                            bump += 1
                            new_id = f"{redsea_bgc_id}__dup{bump}"
                        redsea_bgc_id = new_id

                    seen_ids.add(redsea_bgc_id)

                    rows.append([
                        redsea_bgc_id,
                        genome_name,
                        file,
                        record.id,
                        cluster_type,
                        start,
                        end,
                        strand,
                        os.path.relpath(filepath, root_dir)
                    ])

        except Exception as e:
            failed += 1
            print(f"Failed to parse {filepath}: {e}")

with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "RedSea_BGC_ID",
        "Genome",
        "GBK_File",
        "Contig_ID",
        "Cluster_Type",
        "Start",
        "End",
        "Strand",
        "Relative_Path"
    ])
    writer.writerows(rows)

print(f"Done. Wrote {len(rows)} rows to {output_csv}")
print(f"Unique RedSea_BGC_IDs: {len(seen_ids)}")
print(f"Failed files: {failed}")
