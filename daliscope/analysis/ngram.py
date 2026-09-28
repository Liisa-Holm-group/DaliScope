import math
import re
import numpy as np
import pandas as pd
import pydivsufsort
from collections import defaultdict, Counter

import daliscope.analysis.desc_stemming

import numpy as np

def select_representatives(
    original_descriptions: list[str],
    tfidf_scores: dict[str, float],
    top_k: int = 20,
    processed_descriptions: list[str] = None,
    decay_factor: float = 0.0,  # 0.0 = zero-out covered n-grams; 0.5 = halve their weight
) -> list[tuple[str, float]]:
    """
    Greedily selects top representative descriptions using submodular coverage decay.
    Prevents repeated selection of the same functional keyphrase.
    """
    if processed_descriptions is None:
        processed_descriptions = original_descriptions

    # Make a mutable copy of n-gram weights so we don't mutate external state
    current_tfidf = dict(tfidf_scores)

    selected_reps = []
    selected_indices = set()

    for _ in range(min(top_k, len(original_descriptions))):
        # 1. Rescore remaining candidate descriptions with current (damped) TF-IDF weights
        candidates = []
        for idx, (orig, proc) in enumerate(zip(original_descriptions, processed_descriptions)):
            if idx in selected_indices:
                continue

            # Score using current damped TF-IDF
            score = _score_single_description(proc, current_tfidf)
            candidates.append((score, orig, proc, idx))

        if not candidates:
            break

        # 2. Sort deterministically by (score, original_text)
        candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
        best_score, best_orig, best_proc, best_idx = candidates[0]

        if best_score <= 0.0:
            break  # Stop if no remaining informative descriptions

        selected_reps.append((best_orig, best_score))
        selected_indices.add(best_idx)

        # 3. DAMPING STEP: Decay weights of all n-grams covered by the chosen description
        formatted_proc = " " + best_proc.strip()
        for ngram in list(current_tfidf.keys()):
            if current_tfidf[ngram] > 0 and ngram in formatted_proc:
                current_tfidf[ngram] *= decay_factor

    return selected_reps


def _score_single_description(processed_desc: str, tfidf_dict: dict[str, float]) -> float:
    """Helper to compute mean max-TF-IDF for a single processed description."""
    formatted_desc = " " + processed_desc.strip()
    words = formatted_desc.split()
    if not words:
        return 0.0

    char_scores = np.zeros(len(formatted_desc))

    for ngram, val in tfidf_dict.items():
        if val <= 0:
            continue
        start = 0
        while True:
            idx = formatted_desc.find(ngram, start)
            if idx == -1:
                break
            char_scores[idx : idx + len(ngram)] = np.maximum(
                char_scores[idx : idx + len(ngram)], val
            )
            start = idx + 1

    word_scores = []
    curr_pos = 0
    for word in words:
        w_idx = formatted_desc.find(word, curr_pos)
        w_len = len(word)
        if w_idx != -1:
            word_scores.append(char_scores[w_idx : w_idx + w_len].max())
            curr_pos = w_idx + w_len

    return float(np.mean(word_scores)) if word_scores else 0.0

# =====================================================================
# 1. CORE FUNCTIONS (Steps 1 - 8)
# =====================================================================

def extract_maximal_word_ngrams_fast( # NOT USED
    descriptions: list[str],
    min_length: int = 3,
    min_freq: int = 2,
) -> set[str]:

    formatted_descs = [" " + d.strip() for d in descriptions]
    DELIM = "\x00"
    buffer = DELIM.join(formatted_descs) + DELIM
    buffer_bytes = buffer.encode("utf-8")

    sa = pydivsufsort.divsufsort(buffer_bytes)
    lcp = pydivsufsort.kasai(buffer_bytes, sa)

    candidate_ngrams = set()
    i = 0

    while i < len(lcp):
        if lcp[i] < min_length:
            i += 1
            continue

        l_val = lcp[i]
        start_i = i

        while i < len(lcp) and lcp[i] >= l_val:
            i += 1

        freq = (i - start_i) + 1

        if freq < min_freq:
            continue

        match_pos = sa[start_i]
        ngram = buffer[match_pos:match_pos + l_val]

        if DELIM in ngram:
            continue

        # Must begin at a word boundary.
        if not ngram.startswith(" "):
            continue

        candidate_ngrams.add(ngram)

    return candidate_ngrams


def extract_maximal_ngrams(
    descriptions: list[str],
    min_length: int = 3,
    min_freq: int = 2,
) -> set[str]:
    """
    Extract maximal repeated n-grams consisting of complete words.

    Parameters
    ----------
    descriptions
        Input descriptions. Each description is treated as a sequence
        of whitespace-delimited words.

    min_length
        Minimum number of complete words in an n-gram.

    min_freq
        Minimum number of suffix occurrences supporting the n-gram.

    Returns
    -------
    set[str]
        Maximal repeated word n-grams.

    Notes
    -----
    The suffix array is constructed over integer word IDs rather than
    characters. Consequently, the LCP array is measured in complete
    words, not characters.

    A phrase is "maximal" if it cannot be extended by another word
    while still satisfying min_freq.
    """

    # ------------------------------------------------------------------
    # 1. Tokenize descriptions into complete words
    # ------------------------------------------------------------------
    tokenized = [
        d.strip().split()
        for d in descriptions
    ]

    # Remove empty descriptions
    tokenized = [words for words in tokenized if words]

    if not tokenized:
        return set()

    # ------------------------------------------------------------------
    # 2. Encode words as integer IDs
    #
    # 0 is reserved as a description separator.
    # ------------------------------------------------------------------
    word_to_id = {}
    next_id = 1

    encoded = []
    for words in tokenized:
        ids = []

        for word in words:
            if word not in word_to_id:
                word_to_id[word] = next_id
                next_id += 1

            ids.append(word_to_id[word])

        encoded.extend(ids)
        encoded.append(0)  # separator

    sequence = np.asarray(encoded, dtype=np.int64)

    # ------------------------------------------------------------------
    # 3. Suffix array + LCP over WORD IDs
    # ------------------------------------------------------------------
    sa = pydivsufsort.divsufsort(sequence)
    lcp = pydivsufsort.kasai(sequence, sa)

    # Reverse vocabulary for reconstructing phrases
    id_to_word = {v: k for k, v in word_to_id.items()}

    candidate_ngrams = set()

    # ------------------------------------------------------------------
    # 4. Find LCP intervals
    #
    # lcp[i] is now the number of COMPLETE WORDS shared by adjacent
    # suffixes in suffix-array order.
    # ------------------------------------------------------------------
    i = 0

    while i < len(lcp):

        if lcp[i] < min_length:
            i += 1
            continue

        l_val = int(lcp[i])
        start_i = i

        while i < len(lcp) and lcp[i] >= l_val:
            i += 1

        # Number of suffixes sharing this prefix
        freq = (i - start_i) + 1

        if freq < min_freq:
            continue

        # Representative suffix
        suffix_pos = int(sa[start_i])

        # Extract exactly l_val COMPLETE WORD TOKENS
        phrase_ids = sequence[
            suffix_pos:suffix_pos + l_val
        ]

        # Do not allow a phrase to cross a description boundary.
        if 0 in phrase_ids:
            continue

        phrase = " ".join(
            id_to_word[int(word_id)]
            for word_id in phrase_ids
        )

        candidate_ngrams.add(phrase)

    return candidate_ngrams

    # ------------------------------------------------------------------
    # 5. Keep only maximal phrases.
    #
    # If both:
    #
    #     A B
    #     A B C
    #
    # occur frequently enough, retain only "A B C".
    # ------------------------------------------------------------------
#    candidates = sorted(
#        candidate_ngrams,
#        key=lambda x: len(x.split()),
#        reverse=True,
#    )
#
#    maximal = set()

    for phrase in candidates:
        # A candidate is non-maximal if it is a prefix of an already
        # retained longer candidate.
        if any(
            longer.startswith(phrase + " ")
            for longer in maximal
        ):
            continue

        maximal.add(phrase)

    return maximal

def get_representative_descriptions(
    FULL_DF: pd.DataFrame,
    view_name: str,
    project=None,
    top_k: int = 10,
    min_length: int = 3,
    min_freq: int = 2,
    verbose: bool = False,
) -> list[tuple[str, float]]:
    """
    Extracts maximal n-grams and selects representative descriptions deterministically.
    """
    # ------------------------------------------------------------
    # 1. Background Corpus Preprocessing & N-Gram Extraction
    # ------------------------------------------------------------
    bg_texts = (
        FULL_DF["description"]
        .dropna()
        .astype(str)
        .map(daliscope.analysis.desc_stemming.trim_first_word)
        .map(preprocess_label)
        .map(deduplicate_by_displacement_peak)
    )

    # Step 1-4: Extract maximal n-grams and SORT to guarantee determinism
    raw_ngrams = extract_maximal_ngrams(
        bg_texts, min_length=min_length, min_freq=min_freq
    )
    # Sorting ensures identical processing order across all runs
    candidate_ngrams = sorted(raw_ngrams)

    if verbose:
        print(f"Extracted {len(candidate_ngrams)} Maximal N-Grams:")
        for ng in candidate_ngrams:
            print(f"  - '{ng}'")
        print("\n" + "=" * 50 + "\n")

    # ------------------------------------------------------------
    # 2. Target View Preprocessing
    # ------------------------------------------------------------
    if project is None:
        project = globals()["project"]

    module_df = project.views[view_name].copy()
    module_df["description"] = module_df["description"].map(
        deduplicate_by_displacement_peak
    )

    original_texts = module_df["description"].dropna().astype(str).tolist()

    module_texts = (
        module_df["description"]
        .dropna()
        .astype(str)
        .map(daliscope.analysis.desc_stemming.trim_first_word)
        .map(preprocess_label)
        .map(deduplicate_by_displacement_peak)
    )

    processed_texts = module_texts.tolist()

    # ------------------------------------------------------------
    # 3. Deterministic TF-IDF Computation & Rescoring
    # ------------------------------------------------------------
    # Step 5: Compute TF-IDF
    tfidf_scores = compute_tfidf(candidate_ngrams, module_texts, bg_texts)

    if verbose:
        print("N-Gram TF-IDF Scores:")
        for ng, score in sorted(
            tfidf_scores.items(), key=lambda x: (x[1], x[0]), reverse=True
        ):
            print(f"  - '{ng}': {score:.4f}")
        print("\n" + "=" * 50 + "\n")

    # Step 7: Rescore descriptions
    rescored_descs = rescore_descriptions_tfidf(module_texts, tfidf_scores)
    # Tie-break by text string to ensure strict ordering
    rescored_descs.sort(key=lambda x: (x[1], x[0]), reverse=True)

    if verbose:
        print("Description Scores (Mean Max-TF-IDF per Word):")
        for desc, score in rescored_descs:
            print(f"  - [{score:.4f}] {desc}")
        print("\n" + "=" * 50 + "\n")

    # ------------------------------------------------------------
    # 4. Representative Selection
    # ------------------------------------------------------------
    # Step 8: Select representative descriptions
    reps = select_representatives(
        original_texts,
        tfidf_scores,
        top_k=top_k,
        processed_descriptions=processed_texts,
    )

    print(f"Selected Representative Descriptions for '{view_name}' with Relevance Score:")
    for rank, (desc, score) in enumerate(reps, 1):
        print(f"  {rank}. [{score:.4f}] {desc}")

    return reps


def _score_single_description(processed_desc: str, tfidf_dict: dict[str, float]) -> float:
    """Helper to compute mean max-TF-IDF for a single processed description."""
    formatted_desc = " " + processed_desc.strip()
    words = formatted_desc.split()
    if not words:
        return 0.0

    char_scores = np.zeros(len(formatted_desc))

    for ngram, val in tfidf_dict.items():
        if val <= 0:
            continue
        start = 0
        while True:
            idx = formatted_desc.find(ngram, start)
            if idx == -1:
                break
            char_scores[idx : idx + len(ngram)] = np.maximum(
                char_scores[idx : idx + len(ngram)], val
            )
            start = idx + 1

    word_scores = []
    curr_pos = 0
    for word in words:
        w_idx = formatted_desc.find(word, curr_pos)
        w_len = len(word)
        if w_idx != -1:
            word_scores.append(char_scores[w_idx : w_idx + w_len].max())
            curr_pos = w_idx + w_len

    return float(np.mean(word_scores)) if word_scores else 0.0

# Default stopword set to fall back on if none provided

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
    "unplaced",
    "unplaced", "genomic", "scaffold", "supercont", "whole", 
    "genome", "shotgun", "sequence", "uncharacterized", "protein",
    "putative", "probable", "hypothetical", "predicted",
}

def preprocess_label(label: str, stopwords=None) -> str:
    """Clean one module_df['label'] entry and return a normalized string.

    Expected label format is roughly:
        target_id  DESCRIPTION

    Examples of removed boilerplate:
        c9kdA
        RAT:AF-M0RDF1-F1
        swissprot/AFDB:AF-Q9XXXX-F1

    Operations:
    1. Strip target_id and database accession boilerplate.
    2. Split on runs of punctuation ([^\\w\\s]+).
    3. Reject any token containing a digit.
    4. Remove stopwords.
    5. Return a normalized space-delimited string.

    Returns
    -------
    str
        Cleaned text string.
    """
    if pd.isna(label):
        return ""

    text = str(label).strip().lower()

    # ------------------------------------------------------------
    # 1. Remove target_id (first whitespace-delimited field)
    # ------------------------------------------------------------
    text = re.sub(r"^\S+\s+", "", text)

    # ------------------------------------------------------------
    # 2. Remove database/accession prefixes
    # ------------------------------------------------------------
    text = re.sub(
        r"\b[a-z0-9_./-]+:af-[a-z0-9]+-f\d+\b",
        " ",
        text,
        flags=re.IGNORECASE,
    )

    text = re.sub(
        r"\b(?:swissprot|afdb|afdb/|uniprot|alphafold)[_/]?(?:af-)?"
        r"[a-z0-9]+-f\d+\b",
        " ",
        text,
        flags=re.IGNORECASE,
    )

    # ------------------------------------------------------------
    # 3. Split on all punctuation (replaces runs of non-alphanumeric with space)
    # ------------------------------------------------------------
    text = re.sub(r"[^\w\s]+", " ", text)

    # ------------------------------------------------------------
    # 4. Tokenize & filter out digits and stopwords
    # ------------------------------------------------------------
    if stopwords is None:
        stopwords = DEFAULT_STOPWORDS
    else:
        stopwords = set(stopwords)

    words = text.split()

    filtered_words = [
        w for w in words
        if not re.search(r"\d", w) and w not in stopwords
    ]

    # ------------------------------------------------------------
    # 5. Return cleaned, normalized space-delimited string
    # ------------------------------------------------------------
    return " ".join(filtered_words)

def compute_tfidf(
    ngrams: set[str],
    target_module_descs: list[str],
    background_corpus: list[str],
) -> dict[str, float]:
    """
    Compute module-specific TF-IDF weights for word n-grams.

    For each n-gram, module TF is the number of target-module
    descriptions containing that n-gram. This rewards phrases that
    recur across the module while avoiding inflation from repeated
    occurrences within a single description.

    IDF is computed from the background corpus using smoothed IDF.
    Sublinear scaling is applied to module-level prevalence.
    """
    # 1. Compute Term Frequency (TF) within the target module
    tf_counts = Counter()
    for doc in target_module_descs:
        #formatted_doc = " " + doc.strip()
        formatted_doc = doc.strip()
        for ngram in ngrams:
            if ngram in formatted_doc:
                tf_counts[ngram] += 1

    # 2. Compute Document Frequency (DF) across the background corpus
    total_bg_docs = len(background_corpus)
    df_counts = Counter()
    for doc in background_corpus:
        #formatted_doc = " " + doc.strip()
        formatted_doc = doc.strip()
        for ngram in ngrams:
            if ngram in formatted_doc:
                df_counts[ngram] += 1
#    print(target_module_descs.value_counts())
#    print("\ntf_counts:", len(tf_counts), tf_counts)
#    print("\ndf_counts:", len(df_counts), df_counts)

    # 3. Combine into Sublinear TF-IDF
    tfidf_scores = {}
    for ngram in ngrams:
        tf = tf_counts[ngram]
        if tf == 0:
            tfidf_scores[ngram] = 0.0
            continue

        # Sublinear TF scaling: 1 + log(tf)
        tf_weight = 1.0 + math.log(tf)

        # Standard smoothed IDF
        df = df_counts[ngram]
        idf = math.log((total_bg_docs + 1) / (df + 1)) + 1.0

        tfidf_scores[ngram] = tf_weight * idf

    return tfidf_scores

def rescore_descriptions_tfidf(
    descriptions: list[str],
    tfidf_dict: dict[str, float]
) -> list[tuple[str, float]]:
    """
    Step 7 (TF-IDF Version):
    Scores words using the max (or sum) TF-IDF of covering n-grams, 
    then averages over words in the description.
    """
    scored_results = []

    for desc in descriptions:
        formatted_desc = " " + desc.strip()
        words = formatted_desc.split()
        if not words:
            scored_results.append((desc, 0.0))
            continue

        char_scores = np.zeros(len(formatted_desc))

        for ngram, tfidf_val in tfidf_dict.items():
            start = 0
            while True:
                idx = formatted_desc.find(ngram, start)
                if idx == -1:
                    break
                # Assign maximum TF-IDF value to character span covered by n-gram
                char_scores[idx : idx + len(ngram)] = np.maximum(
                    char_scores[idx : idx + len(ngram)], tfidf_val
                )
                start = idx + 1

        # Map character scores to words
        word_scores = []
        curr_pos = 0
        for word in words:
            word_idx = formatted_desc.find(word, curr_pos)
            w_len = len(word)
            word_max_score = char_scores[word_idx : word_idx + w_len].max()
            word_scores.append(word_max_score)
            curr_pos = word_idx + w_len

        #desc_score = float(np.mean(word_scores)) if word_scores else 0.0
        desc_score = float(np.sum(word_scores)) if word_scores else 0.0
        scored_results.append((desc, desc_score))

    return scored_results

def deduplicate_by_displacement_peak(text: str, min_peak_size: int = 3) -> str:
    """
    Detects repeated/truncated text blocks by finding off-diagonal 
    word position displacement peaks (k-mer dot-plot logic).
    """
    words = text.split()
    n = len(words)
    if n < 4:
        return text

    # 1. Map lowercased words to their positional indices
    pos_map = defaultdict(list)
    for idx, word in enumerate(words):
        pos_map[word.lower()].append(idx)

    # 2. Compute pairwise index displacements (Pos2 - Pos1)
    displacement_counts = Counter()

    for word, positions in pos_map.items():
        if len(positions) > 1:
            # Check all forward offsets for this word
            for i in range(len(positions)):
                for j in range(i + 1, len(positions)):
                    shift = positions[j] - positions[i]
                    displacement_counts[shift] += 1

    if not displacement_counts:
        return text

    # 3. Find the dominant off-diagonal peak
    most_common_shift, peak_count = displacement_counts.most_common(1)[0]

    # If the peak is significant (e.g., at least 3 matching displaced words),
    # slice the string at the peak shift offset (where the repeat starts)
    if peak_count >= min_peak_size and most_common_shift > 0:
        return " ".join(words[:most_common_shift])

    return text
