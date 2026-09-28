#!/usr/bin/env python3
import sys
import sqlite3
import pandas as pd

def print_usage_and_exit(exit_code=1):
    usage_text = """
Usage:
    python extract_mini_ligands.py <colab_tsv> <PDB_ligands_db> <out_mini_ligand_tsv> <out_mini_ligand_names_tsv>

Arguments:
    colab_tsv                Path to colab_pack_index.tsv (containing 'dali_id' column)
    PDB_ligands_db           Path to SQLite database containing 'ligand' and 'ligand_names' tables
    out_mini_ligand_tsv      Path to output TSV for matched ligand HETATMs
    out_mini_ligand_names_tsv Path to output TSV for matched ligand annotations

Example:
    python extract_mini_ligands.py colab_index.tsv PDB_ligands.db mini_ligand.tsv mini_ligand_names.tsv
"""
    print(usage_text.strip(), file=sys.stderr if exit_code != 0 else sys.stdout)
    sys.exit(exit_code)

def main():
    # Show usage help if requested or incorrect argument count
    if len(sys.argv) in (2, 3, 4, 5) and sys.argv[1] in ("-h", "--help"):
        print_usage_and_exit(exit_code=0)

    if len(sys.argv) != 5:
        print("Error: Invalid number of arguments.\n", file=sys.stderr)
        print_usage_and_exit(exit_code=1)

    colab_tsv = sys.argv[1]
    PDB_ligands_db = sys.argv[2]
    out_mini_ligand_tsv = sys.argv[3]
    out_mini_ligand_names = sys.argv[4]

    # Parse DALI hits from file
    try:
        df_dali = pd.read_csv(colab_tsv, sep="\t")
        query_pairs = [(cd1[0:4], cd1[4]) for cd1 in df_dali['dali_id']]
    except Exception as e:
        print(f"Error reading input TSV file '{colab_tsv}': {e}", file=sys.stderr)
        sys.exit(1)

    conn = sqlite3.connect(PDB_ligands_db)
    cursor = conn.cursor()

    try:
        # 1. Populate temporary table with query pairs
        cursor.execute("CREATE TEMP TABLE temp_keys (pdb_id TEXT, chain_id TEXT);")
        cursor.executemany("INSERT INTO temp_keys VALUES (?, ?);", query_pairs)
        cursor.execute("CREATE INDEX idx_temp_keys ON temp_keys (pdb_id, chain_id);")

        # 2. Extract matching ligand rows from SQLite
        sql_ligand_query = """
            SELECT l.*
            FROM ligand l
            INNER JOIN temp_keys t
                    ON l.pdb_id = t.pdb_id
                   AND l.chain_id = t.chain_id;
        """
        df_matches = pd.read_sql_query(sql_ligand_query, conn)
        print(f"Loaded {len(df_matches)} matching ligand rows.")

        # Save mini_ligand.tsv
        df_matches.to_csv(out_mini_ligand_tsv, sep="\t", index=False)

        # 3. Fetch matching annotations from SQLite's ligand_names table
        unique_res_names = df_matches["res_name"].dropna().unique().tolist()

        if unique_res_names:
            placeholders = ", ".join(["?"] * len(unique_res_names))
            sql_names_query = f"""
                SELECT DISTINCT res_name, compound
                FROM ligand_names
                WHERE res_name IN ({placeholders});
            """
            df_names_matched = pd.read_sql_query(sql_names_query, conn, params=unique_res_names)

            # Left join to preserve any res_name keys missing in ligand_names table
            keys_df = pd.DataFrame({"res_name": unique_res_names})
            result_names_df = keys_df.merge(df_names_matched, on="res_name", how="left")
        else:
            result_names_df = pd.DataFrame(columns=["res_name", "compound"])

        # Save mini_ligand_names.tsv
        result_names_df.to_csv(out_mini_ligand_names, sep="\t", index=False)
        print(f"Extracted {len(result_names_df)} unique ligand annotations.")

    finally:
        cursor.execute("DROP TABLE IF EXISTS temp_keys;")
        conn.close()

if __name__ == "__main__":
    main()
