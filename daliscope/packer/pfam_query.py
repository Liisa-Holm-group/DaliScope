"""Extract Pfam hits for packed targets from a local annotation database."""
import argparse
import csv
from pathlib import Path
import sqlite3

import pandas as pd

NAMESPACE = {1: "PDB", 2: "AFDB2", 3: "BFVD", 4: "VIRO3D"}


def extract_hits(db_path, dali_list, ns_id, output_tsv):
    namespace = NAMESPACE[int(ns_id)]
    database = Path(db_path).resolve()
    if not database.is_file():
        raise FileNotFoundError(f"Pfam database not found: {database}")
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True, timeout=30) as conn:
        conn.execute("CREATE TEMP TABLE filter_ids (dali_id TEXT PRIMARY KEY)")
        conn.executemany("INSERT OR IGNORE INTO filter_ids VALUES (?)", ((str(d),) for d in dali_list))
        cursor = conn.execute(
            "SELECT h.*, c.CL AS clan FROM hmmer_hits h "
            "LEFT JOIN pfam_entries c ON h.query_accession = c.AC "
            "INNER JOIN filter_ids f ON h.target_name = f.dali_id "
            "WHERE h.database = ?", (namespace,)
        )
        with open(output_tsv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(column[0] for column in cursor.description)
            while rows := cursor.fetchmany(5000):
                writer.writerows(rows)
    print(f"Extracted Pfam hits to {output_tsv}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("namespace", type=int, choices=tuple(NAMESPACE))
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--index", type=Path, default=Path("colab_pack_index.tsv"))
    parser.add_argument("--output", type=Path, default=Path("pfam.tsv"))
    args = parser.parse_args(argv)
    df = pd.read_csv(args.index, sep="\t")
    extract_hits(args.database, df["dali_id"].tolist(), args.namespace, args.output)


if __name__ == "__main__":
    main()
