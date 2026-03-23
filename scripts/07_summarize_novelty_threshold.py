# Description:
# Summarizes BiG-SLiCE novelty by extracting the minimum distance
# from each Red Sea BGC to the global reference set.
# BGCs with minimum distance > 900 are classified as novel.

import pandas as pd
import numpy as np

distance_path = "distance_matrix_ids.tsv"
out_path = "novelty_summary_bgc_ids.tsv"
threshold = 900

header = pd.read_csv(distance_path, sep="\t", nrows=0)
cols = header.columns.tolist()

id_col = cols[0]
dist_cols = cols[1:]

chunk_size = 50

min_vals = None
bgc_ids = None

for i in range(0, len(dist_cols), chunk_size):
    usecols = [id_col] + dist_cols[i:i + chunk_size]
    chunk = pd.read_csv(distance_path, sep="\t", usecols=usecols)

    ids = chunk[id_col].to_numpy()
    block = chunk.drop(columns=[id_col]).to_numpy(dtype=np.float32)

    block_min = np.nanmin(block, axis=1)

    if min_vals is None:
        min_vals = block_min
        bgc_ids = ids
    else:
        min_vals = np.minimum(min_vals, block_min)

novel_flags = min_vals > threshold

novelty_summary = pd.DataFrame({
    "RedSea_BGC_ID": bgc_ids,
    "Min_Distance": min_vals,
    "Is_Novel": novel_flags
})

novelty_summary.to_csv(out_path, sep="\t", index=False)

total = len(novelty_summary)
novel = int(novel_flags.sum())

print(f"Total BGCs: {total}")
print(f"Novel BGCs (distance > {threshold}): {novel}")
print(f"Novelty percentage: {100 * novel / total:.2f}%")
