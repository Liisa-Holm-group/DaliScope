"""Extract Pfam hits for packed targets from a local annotation database."""
import argparse
import csv
from pathlib import Path
import sqlite3

import pandas as pd

NAMESPACE = {1: "PDB", 2: "AFDB2", 3: "BFVD", 4: "VIRO3D", 5: "LOCAL"}

# Pfam 39.0 databases keep each hit together with its clan in hmmer_hits_with_CL.
# Older databases (Pfam 38.2) have hmmer_hits and take the clan from pfam_entries.
# Namespace labels are compared without case: the databases store "Viro3D".
HITS_WITH_CLAN = (
    "SELECT h.*, h.CL AS clan FROM hmmer_hits_with_CL h "
    "INNER JOIN filter_ids f ON h.target_name = f.dali_id "
    "WHERE h.database = ? COLLATE NOCASE"
)
HITS_JOINED = (
    "SELECT h.*, c.CL AS clan FROM hmmer_hits h "
    "LEFT JOIN pfam_entries c ON h.query_accession = c.AC "
    "INNER JOIN filter_ids f ON h.target_name = f.dali_id "
    "WHERE h.database = ? COLLATE NOCASE"
)


def hits_query(conn):
    """Return the hit query that matches the database schema."""
    tables = {row[0].lower() for row in conn.execute("SELECT name FROM sqlite_master WHERE type IN ('table', 'view')")}
    if "hmmer_hits_with_cl" in tables:
        return HITS_WITH_CLAN
    if {"hmmer_hits", "pfam_entries"} <= tables:
        return HITS_JOINED
    raise ValueError("Pfam database has neither hmmer_hits_with_CL nor hmmer_hits with pfam_entries")


def extract_hits(db_path, dali_list, ns_id, output_tsv):
    namespace = NAMESPACE[int(ns_id)]
    database = Path(db_path).resolve()
    if not database.is_file():
        raise FileNotFoundError(f"Pfam database not found: {database}")
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True, timeout=30) as conn:
        query = hits_query(conn)
        conn.execute("CREATE TEMP TABLE filter_ids (dali_id TEXT PRIMARY KEY)")
        conn.executemany("INSERT OR IGNORE INTO filter_ids VALUES (?)", ((str(d),) for d in dali_list))
        cursor = conn.execute(query, (namespace,))
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
    df = pd.read_csv(args.index, sep="\t", dtype={"dali_id": str}, keep_default_na=False)
    extract_hits(args.database, df["dali_id"].tolist(), args.namespace, args.output)


if __name__ == "__main__":
    main()
