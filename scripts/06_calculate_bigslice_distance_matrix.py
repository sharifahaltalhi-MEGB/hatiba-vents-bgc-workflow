# Description:
# Computes Euclidean distances between Red Sea BGC feature matrix
# and a global reference feature matrix from BiG-SLiCE.
# Outputs a distance matrix used for novelty assessment.

import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances

# Input files (update paths as needed)
your_matrix_path = "RedSea_features_matrix.tsv"
public_matrix_path = "global_features_matrix.tsv"
output_distance_path = "distance_matrix_ids.tsv"

print("Loading feature matrices...")
your_df = pd.read_csv(your_matrix_path, sep="\t", index_col=0)
public_df = pd.read_csv(public_matrix_path, sep="\t", index_col=0)

print("Aligning feature columns...")
common_features = your_df.columns.intersection(public_df.columns)
your_df = your_df[common_features]
public_df = public_df[common_features]

print("Calculating Euclidean distances...")
distances = euclidean_distances(your_df, public_df)

print("Saving distance matrix...")
distance_df = pd.DataFrame(distances, index=your_df.index, columns=public_df.index)
distance_df.to_csv(output_distance_path, sep="\t")

print(f"Distance matrix saved to {output_distance_path}")
