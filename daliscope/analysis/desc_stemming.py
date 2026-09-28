import re
import numpy as np
import pandas as pd

from collections import Counter
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


DEFAULT_STOPWORDS = {
    # generic annotation boilerplate
    "protein",
    "putative",
    "probable",
    "possible",
    "uncharacterized",
    "hypothetical",
    "unknown",
    "predicted",
    "candidate",
    "related",
    "like",
    "family",
    "domain",
    "containing",
    "member",
    "fragment",
    "isoform",
    "precursor",
    "chain",
    "and",
    "domain-containing",
    "unassigned",

    # database / annotation boilerplate
    "orf",
    "open",
    "reading",
    "frame",
    "genome",
    "genomic",
    "scaffold",
    "whole",
    "shotgun",
    "sequence",
    "unplaced"
}


def preprocess_label(label, stopwords=None):
    """
    Clean one module_df['label'] entry and return normalized tokens.

    Expected label format is roughly:
        target_id  DESCRIPTION

    Examples of removed boilerplate:
        c9kdA
        RAT:AF-M0RDF1-F1
        swissprot/AFDB:AF-Q9XXXX-F1

    Hyphenated biological terms are preserved, e.g.
        2-oxoglutarate-dependent
        iron-dependent
        lysine-specific

    Returns
    -------
    list[str]
        Cleaned tokens.
    """
    if pd.isna(label):
        return []

    text = str(label).strip().lower()

    # ------------------------------------------------------------
    # Remove target_id (first whitespace-delimited field)
    # ------------------------------------------------------------
    text = re.sub(r"^\S+\s+", "", text)

    # ------------------------------------------------------------
    # Remove database/accession prefixes.
    #
    # Handles examples such as:
    #   rat:af-m0rdf1-f1
    #   arath:af-q9fi36-f1
    #   swissprot/afdb:af-q4wkx0-f1
    #   9euro2/afdb:af-a0a0d2gl62-f1
    # ------------------------------------------------------------
    text = re.sub(
        r"\b[a-z0-9_./-]+:af-[a-z0-9]+-f\d+\b",
        " ",
        text,
        flags=re.IGNORECASE,
    )

    # Some labels have database/accession text without a colon
    text = re.sub(
        r"\b(?:swissprot|afdb|afdb/|uniprot|alphafold)[_/]?(?:af-)?"
        r"[a-z0-9]+-f\d+\b",
        " ",
        text,
        flags=re.IGNORECASE,
    )

    # ------------------------------------------------------------
    # Normalize common annotation punctuation.
    # Keep hyphens because they are biologically informative.
    # ------------------------------------------------------------
    text = text.replace("_", " ")
    text = text.replace("/", " ")
    text = text.replace(";", " ")
    text = text.replace(",", " ")
    text = text.replace("(", " ")
    text = text.replace(")", " ")
    text = text.replace("[", " ")
    text = text.replace("]", " ")
    text = text.replace(":", " ")

    # Collapse whitespace
    text = re.sub(r"\s+", " ", text).strip()

    # ------------------------------------------------------------
    # Tokenization
    #
    # Allows:
    #   letters
    #   numbers
    #   embedded hyphens
    #   e.g. 2-oxoglutarate-dependent
    # ------------------------------------------------------------
    tokens = re.findall(
        r"[a-z0-9]+(?:-[a-z0-9]+)*",
        text,
    )

    if stopwords is None:
        stopwords = DEFAULT_STOPWORDS

    tokens = [
        token
        for token in tokens
        if token not in stopwords
    ]

    return tokens


def add_clean_text(
    modules_df,
    label_col="label",
    output_col="clean_tokens",
    stopwords=None,
):
    """
    Add preprocessed tokens to modules_df.

    Returns a copy of modules_df.
    """
    df = modules_df.copy()

    df[output_col] = df[label_col].map(
        lambda x: preprocess_label(
            x,
            stopwords=stopwords,
        )
    )

    return df

def find_ngram_candidates(
    modules_df,
    module_col="module",
    label_col="label",
    ngram_range=(1, 3),
    min_count=3,
    min_module_count=2,
    min_f1=0.20,
    max_qvalue=0.05,
    stopwords=None,
):
    """
    Find biologically informative n-gram candidates from
    modules_df[['module', 'label']].

    An n-gram is considered in terms of document presence:
    a description either contains it or does not.

    Returns
    -------
    candidates : pd.DataFrame
        Candidate markers ranked within module.
    clean_df : pd.DataFrame
        Input dataframe with cleaned token lists.
    """

    df = add_clean_text(
        modules_df,
        label_col=label_col,
        output_col="_tokens",
        stopwords=stopwords,
    )

    # ------------------------------------------------------------
    # Generate n-gram presence sets per description
    # ------------------------------------------------------------
    presence = []

    for tokens in df["_tokens"]:
        grams = set()

        for n in range(ngram_range[0], ngram_range[1] + 1):
            grams.update(
                " ".join(tokens[i:i+n])
                for i in range(len(tokens) - n + 1)
            )

        presence.append(grams)

    # Global document frequency
    global_counts = Counter()

    for grams in presence:
        global_counts.update(grams)

    valid_ngrams = {
        gram
        for gram, count in global_counts.items()
        if count >= min_count
    }

    modules = sorted(df[module_col].dropna().unique())
    module_values = df[module_col].to_numpy()

    results = []

    # ------------------------------------------------------------
    # Test each n-gram against each module
    # ------------------------------------------------------------
    for gram in valid_ngrams:

        present = np.array(
            [gram in grams for grams in presence],
            dtype=bool,
        )

        total_count = int(present.sum())

        for module in modules:

            in_module = module_values == module

            a = int(np.sum(present & in_module))
            b = int(np.sum(~present & in_module))
            c = int(np.sum(present & ~in_module))
            d = int(np.sum(~present & ~in_module))

            module_size = int(in_module.sum())

            if a < min_module_count:
                continue

            precision = a / module_size
            recall = a / total_count

            f1 = (
                2 * precision * recall / (precision + recall)
                if precision + recall > 0
                else 0.0
            )

            odds_ratio, pvalue = fisher_exact(
                [[a, b], [c, d]],
                alternative="greater",
            )

            results.append({
                "module": module,
                "ngram": gram,
                "n": len(gram.split()),
                "count": a,
                "total_count": total_count,
                "module_size": module_size,
                "precision": precision,
                "recall": recall,
                "F1": f1,
                "odds_ratio": odds_ratio,
                "pvalue": pvalue,
            })

    stats = pd.DataFrame(results)

    if stats.empty:
        return stats, df

    # FDR correction over all module × n-gram tests
    stats["qvalue"] = multipletests(
        stats["pvalue"],
        method="fdr_bh",
    )[1]

    # Apply practical filters
    candidates = stats[
        (stats["F1"] >= min_f1) &
        (stats["qvalue"] <= max_qvalue)
    ].copy()

    # Rank the candidates primarily by F1, then enrichment
    candidates = (
        candidates
        .sort_values(
            ["module", "F1", "odds_ratio", "count"],
            ascending=[True, False, False, False],
        )
        .reset_index(drop=True)
    )

    return candidates, df

from collections import Counter
import pandas as pd


def trim_first_word(text: str) -> str:
    """Drop the first whitespace-delimited token from text."""
    words = text.split()
    return " ".join(words[1:])


def get_stem_prefix(text: str, n_words: int = 2) -> str:
    """Extracts a simple prefix stem taking the first n words, converted to uppercase."""
    words = text.strip().rstrip(";").split()
    if not words:
        return ""
    return " ".join(words[:n_words]).upper()


def summarize_description_stems(df, description_col: str = "description",
                                 n_words: int = 2, min_occurrences: int = 10,
                                 trim_first: bool = False, verbose: bool = True) -> pd.DataFrame:
    """Fixed-length stem aggregation (as before), operating on trimmed text
    (first word dropped) by default."""
    texts = df[description_col].dropna().astype(str)
    if trim_first:
        texts = texts.map(trim_first_word)

    stem_counts = Counter()
    for desc in texts:
        stem = get_stem_prefix(preprocess_label(desc), n_words=n_words)
        if stem:
            stem_counts[stem] += 1

    frequent_stems = {s: c for s, c in stem_counts.items() if c >= min_occurrences}
    sorted_stems = sorted(frequent_stems.items(), key=lambda x: x[1], reverse=True)

    if verbose:
        print(f"{'Stem Prefix':<35} | {'Total Count':<10}")
        print("-" * 50)
        for stem, count in sorted_stems:
            print(f"{stem:<35} | {count:<10}")

    return pd.DataFrame(sorted_stems, columns=["stem", "count"])


def summarize_description_stems_multi(df, description_col: str = "description",
                                       n_words_list=(1, 2, 3, 4), min_occurrences: int = 10,
                                       trim_first: bool = True, verbose: bool = True) -> dict:
    """Option 1: run the fixed-length aggregation at each n_words in
    n_words_list. Returns {n_words: result_df}."""
    results = {}
    for n in n_words_list:
        if verbose:
            print(f"\n=== n_words={n} ===")
        results[n] = summarize_description_stems(
            df, description_col=description_col, n_words=n,
            min_occurrences=min_occurrences, trim_first=trim_first, verbose=verbose,
        )
    return results


def summarize_description_stems_maximal(df, description_col: str = "description",
                                         max_words: int = 8, min_occurrences: int = 10,
                                         trim_first: bool = True, verbose: bool = False) -> pd.DataFrame:
    """
    Option 2: variable-length stems via MAXIMAL frequent prefixes.

    Every prefix length (1..max_words, or the full trimmed description if
    shorter) is counted across all rows. A prefix is "frequent" if its
    count >= min_occurrences. Only MAXIMAL frequent prefixes are reported
    -- ones with no strictly longer frequent extension -- since a shorter
    prefix's count is always >= any longer extension's count (every match
    of the longer prefix also matches the shorter one), so the longest
    still-frequent version is the most specific non-redundant stem.

    Returns a DataFrame ["stem", "count", "n_words"], sorted by count
    descending.
    """
    texts = df[description_col].dropna().astype(str)
    if trim_first:
        texts = texts.map(trim_first_word)

    prefix_counts = Counter()
    for desc in texts:
        words = preprocess_label(desc)[:max_words]
        for k in range(1, len(words) + 1):
            prefix_counts[tuple(words[:k])] += 1

    frequent = {p: c for p, c in prefix_counts.items() if c >= min_occurrences}
    frequent_set = set(frequent.keys())

    maximal = []
    for p, c in frequent.items():
        has_longer_extension = any(
            len(q) == len(p) + 1 and q[:len(p)] == p
            for q in frequent_set
        )
        if not has_longer_extension:
            maximal.append((" ".join(p), c, len(p)))

    maximal.sort(key=lambda x: x[1], reverse=True)

    if verbose:
        print(f"{'Stem Prefix':<35} | {'Total Count':<10} | {'N words':<7}")
        print("-" * 58)
        for stem, count, nw in maximal:
            print(f"{stem:<35} | {count:<10} | {nw:<7}")

    return pd.DataFrame(maximal, columns=["stem", "count", "n_words"])
