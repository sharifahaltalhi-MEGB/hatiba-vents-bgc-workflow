import sqlite3
import pandas as pd
import numpy as np
from os import path
from sys import argv

def fetch_bgc_features(result_folder):
    print("Loading BGC features...")
    with sqlite3.connect(path.join(result_folder, "result", "data.db")) as con:
        cur = con.cursor()

        # Get BGC IDs and names
        bgc_rows = cur.execute("SELECT id, name FROM bgc ORDER BY id ASC").fetchall()
        bgc_ids = [row[0] for row in bgc_rows]
        bgc_names = [row[1] for row in bgc_rows]  # row index

        # Get HMM features
        hmm_rows = cur.execute("SELECT id, name FROM hmm WHERE db_id=1 ORDER BY id ASC").fetchall()
        hmm_ids = [row[0] for row in hmm_rows]
        hmm_names = [row[1] for row in hmm_rows]

        # Create empty matrix
        bgc_features = pd.DataFrame(
            np.zeros((len(bgc_ids), len(hmm_ids)), dtype=np.uint8),
            index=bgc_names,
            columns=hmm_names
        )

        bgc_id_to_name = dict(bgc_rows)
        hmm_id_to_name = dict(hmm_rows)

        # Fill matrix
        for bgc_id, hmm_id, value in cur.execute("""
            SELECT bgc_id, hmm_id, value
            FROM bgc_features
            WHERE bgc_features.hmm_id IN (
                SELECT id FROM hmm WHERE db_id=1
            )
        """).fetchall():
            row_name = bgc_id_to_name[bgc_id]
            col_name = hmm_id_to_name[hmm_id]
            bgc_features.at[row_name, col_name] = value

    return bgc_features

def main():
    try:
        bigslice_result_folder = argv[1]
        output_csv_path = argv[2]
    except:
        print("usage: python 05_extract_bigslice_features_matrix.py <bigslice_result_folder> <output_tsv_path>")
        return 1

    bigfam_features = fetch_bgc_features(bigslice_result_folder)
    print("Saving to file...")
    bigfam_features.to_csv(output_csv_path, sep="\t")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
