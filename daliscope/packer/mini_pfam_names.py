"""Select a pack's Pfam descriptions from the bundled or specified reference."""

import argparse
from pathlib import Path

import pandas as pd

DEFAULT_PFAM_NAMES = Path(__file__).resolve().parent / "data" / "pfam_names.tsv"


def create_mini_pfam_names(pfam_tsv_path: str = "pfam.tsv",
                           names_tsv_path: str | Path = DEFAULT_PFAM_NAMES,
                           output_path: str = "mini_pfam_names.tsv"):
    df_pfam = pd.read_csv(pfam_tsv_path, sep="\t")
    df_names = pd.read_csv(names_tsv_path, sep="\t")
    active_pfams = set(df_pfam["query_accession"].dropna().unique())
    mini_names_df = df_names[df_names["pfam"].isin(active_pfams)]
    mini_names_df.to_csv(output_path, sep="\t", index=False)
    print("Saved {} matching Pfam rows to {}".format(len(mini_names_df), output_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("names_tsv", nargs="?", type=Path, default=DEFAULT_PFAM_NAMES)
    parser.add_argument("output_tsv", nargs="?", type=Path, default=Path("mini_pfam_names.tsv"))
    args = parser.parse_args()
    create_mini_pfam_names(names_tsv_path=args.names_tsv, output_path=args.output_tsv)
