import pandas as pd
import sys

def create_mini_pfam_names(pfam_tsv_path: str = "pfam.tsv", 
                           names_tsv_path: str = "../../pfamdata/pfam_names.tsv", 
                           output_path: str = "mini_pfam_names.tsv"):
    # Load input files (tab-delimited)
    df_pfam = pd.read_csv(pfam_tsv_path, sep='\t')
    df_names = pd.read_csv(names_tsv_path, sep='\t')

    # Extract unique query_accessions present in pfam.tsv
    active_pfams = set(df_pfam['query_accession'].dropna().unique())

    # Filter strictly on matching pfam IDs
    mini_names_df = df_names[df_names['pfam'].isin(active_pfams)]

    # Write to output TSV
    mini_names_df.to_csv(output_path, sep='\t', index=False)
    print("Saved {} matching Pfam rows to {}".format(len(mini_names_df), output_path))

# Execute
if __name__ == "__main__":
    pfam_tsv = 'pfam.tsv'
    names_tsv_path = sys.argv[1]
    output_path = sys.argv[2]
    create_mini_pfam_names(pfam_tsv_path=pfam_tsv, names_tsv_path=names_tsv_path, output_path=output_path)
