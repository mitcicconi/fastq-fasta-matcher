"""K-mer containment matching, generalized from *pipeline_kmer.py.

The original script used fixed 20bp junction barcodes specific to one TadA
insertion assay. This version generalizes that idea to arbitrary references:
build a k-mer set per reference, then score each read by how much it overlaps
that reference's k-mer set (its "containment"). This is orientation-independent
(checks both strands) and doesn't require any assay-specific anchor sequence.

Containment is normalized by the SMALLER of the read's/reference's k-mer
counts, not always the read's. Always dividing by the read's k-mer count
would make containment collapse toward zero whenever the read is much longer
than the reference (e.g. a 2kb amplicon read with only a 600bp target
region), even for a perfect match in the overlapping part.
"""

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def kmer_set(seq, k):
    seq = seq.upper()
    if len(seq) < k:
        return set()
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def build_reference_kmer_index(fasta_records, k):
    """fasta_records: list of (name, sequence). Returns {name: kmer_set}."""
    return {name: kmer_set(seq, k) for name, seq in fasta_records}


def match_reads_by_kmer(fastq_records, ref_kmers, k, min_containment=0.3):
    """fastq_records: list of (name, sequence, quality).

    For each read, finds the reference with the highest containment
    (fraction of the read's k-mers found in that reference), checking
    both the read's forward and reverse-complement strand.

    Returns a DataFrame with one row per read.
    """
    import pandas as pd

    rows = []
    for name, seq, _qual in fastq_records:
        read_kmers = kmer_set(seq, k) | kmer_set(revcomp(seq), k)

        if not read_kmers:
            rows.append(
                {"read_id": name, "reference": None, "containment": 0.0, "matched": False}
            )
            continue

        best_ref, best_score = None, 0.0
        for ref_name, rk in ref_kmers.items():
            if not rk:
                continue
            shared = len(read_kmers & rk)
            containment = shared / min(len(read_kmers), len(rk))
            if containment > best_score:
                best_score, best_ref = containment, ref_name

        matched = best_score >= min_containment
        rows.append(
            {
                "read_id": name,
                "reference": best_ref if matched else None,
                "containment": round(best_score, 4),
                "matched": matched,
            }
        )

    return pd.DataFrame(rows)


def summarize_by_reference(per_read_df, ref_names):
    """Collapse per-read k-mer results into a per-reference matched-read count table."""
    import pandas as pd

    counts = per_read_df[per_read_df["matched"]]["reference"].value_counts()
    return pd.DataFrame(
        [
            {
                "reference": r,
                "reads_matched_kmer": int(counts.get(r, 0)),
                "status": "Matched" if counts.get(r, 0) > 0 else "No confident match",
            }
            for r in ref_names
        ]
    ).sort_values("reads_matched_kmer", ascending=False).reset_index(drop=True)
