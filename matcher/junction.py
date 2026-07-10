"""Junction barcode matching, generalized from *pipeline_kmer.py.

The original script hardcoded a single payload (TadA, 614bp identical across
all variants) between two IscB flanks, and used fixed offsets derived from
that exact length. This version works against any "junction reference" FASTA
built the same way — [flank][payload][flank] per record, same payload in
every record of a file — regardless of what the payload is or how long it
is, so the same algorithm serves a TadA-inserted database, a different
construct's database (e.g. an intermediate insertion-site placeholder), or
anything else built with the same flank structure.

Read assignment logic is unchanged from *pipeline_kmer.py: reads are scanned
for both the left and right junction barcodes (20bp of flank immediately
adjacent to the payload) via a 1-mismatch-tolerant hash table, and a read is
counted if at least one side gives an unambiguous call (concordant + left_only
+ right_only) — see project memory for why this is the safe metric to use,
and why 1-mismatch tolerance is the right ceiling, not an arbitrary choice.
"""

from Bio import SeqIO

_BASES = "ACGT"
_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def mismatch_variants(seq, max_mm=1):
    """All sequences within max_mm substitutions of seq (including itself)."""
    result = {seq}
    for i in range(len(seq)):
        for b in _BASES:
            if b != seq[i]:
                result.add(seq[:i] + b + seq[i + 1 :])
    return result


def find_anchor(read_seq, anchor, max_mm=1):
    alen = len(anchor)
    return [
        i
        for i in range(len(read_seq) - alen + 1)
        if hamming(read_seq[i : i + alen], anchor) <= max_mm
    ]


def load_junction_reference(fasta_path, flank=150, bc_len=20):
    """Loads a [flank][payload][flank] FASTA and derives per-variant barcodes.

    Every record must have the same total length (same payload length across
    the file — this holds for a single pool's junction database). Returns
    (left_barcodes, right_barcodes, payload_seq, payload_len) where the two
    barcode dicts map variant id -> 20bp barcode sequence.
    """
    recs = list(SeqIO.parse(fasta_path, "fasta"))
    if not recs:
        raise ValueError(f"No sequences found in {fasta_path}")

    payload_len = len(recs[0].seq) - 2 * flank
    if payload_len <= 0:
        raise ValueError(f"Flank ({flank}bp x2) is longer than the record itself")

    left_barcodes, right_barcodes = {}, {}
    payload_seq = None
    for rec in recs:
        seq = str(rec.seq).upper()
        if len(seq) - 2 * flank != payload_len:
            raise ValueError(
                f"Record {rec.id} has a different payload length than the rest of {fasta_path}"
            )
        left_barcodes[rec.id] = seq[flank - bc_len : flank]
        right_barcodes[rec.id] = seq[flank + payload_len : flank + payload_len + bc_len]
        if payload_seq is None:
            payload_seq = seq[flank : flank + payload_len]

    return left_barcodes, right_barcodes, payload_seq, payload_len


def _build_mismatch_table(barcodes, mm_allow=1):
    """Maps every mm_allow-mismatch variant of every non-conflicting barcode
    to its variant id; collisions between two different barcodes' variant
    sequences are marked None (ambiguous) rather than guessed at.
    """
    inv = {}
    conflicts = set()
    for vid, bc in barcodes.items():
        if bc in inv:
            conflicts.add(bc)
        else:
            inv[bc] = vid

    table = {}
    for vid, bc in barcodes.items():
        if bc in conflicts:
            continue
        for variant_seq in mismatch_variants(bc, mm_allow):
            if variant_seq in table:
                table[variant_seq] = None
            else:
                table[variant_seq] = vid
    return table, conflicts


def match_reads_by_junction(
    fastq_records, left_barcodes, right_barcodes, payload_seq, payload_len, bc_len=20, mm_allow=1
):
    """fastq_records: list of (name, seq, qual).

    Returns (per_variant_rows, qc_dict). per_variant_rows has one row per
    variant with concordant/left_only/right_only counts and their sum
    (mapped_reads) — the metric actually used downstream, matching
    *pipeline_kmer.py's implemented behavior (not its docstring).
    """
    left_mm, _left_conflicts = _build_mismatch_table(left_barcodes, mm_allow)
    right_mm, _right_conflicts = _build_mismatch_table(right_barcodes, mm_allow)

    left_anchor = payload_seq[:bc_len]
    right_anchor = payload_seq[-bc_len:]

    concordant, left_only, right_only = {}, {}, {}
    n_reads = n_concordant = n_left_only = n_right_only = n_conflict = n_no_match = 0

    for name, seq, _qual in fastq_records:
        n_reads += 1
        seq_fwd = seq.upper()
        seq_rev = revcomp(seq_fwd)
        best = None

        for s in (seq_fwd, seq_rev):
            slen = len(s)

            for lpos in find_anchor(s, left_anchor, mm_allow):
                if lpos < bc_len:
                    continue
                left_iscb = s[lpos - bc_len : lpos]
                left_vid = left_mm.get(left_iscb)

                r_start = lpos + payload_len
                right_vid = None
                if r_start + bc_len <= slen:
                    right_iscb = s[r_start : r_start + bc_len]
                    right_vid = right_mm.get(right_iscb)

                if left_vid is not None or right_vid is not None:
                    best = (left_vid, right_vid)
                    break

            if best is None:
                for rpos in find_anchor(s, right_anchor, mm_allow):
                    r_iscb_start = rpos + bc_len
                    if r_iscb_start + bc_len > slen:
                        continue
                    right_iscb = s[r_iscb_start : r_iscb_start + bc_len]
                    right_vid = right_mm.get(right_iscb)

                    l_iscb_end = rpos - (payload_len - bc_len)
                    l_iscb_start = l_iscb_end - bc_len
                    left_vid = None
                    if l_iscb_start >= 0:
                        left_iscb = s[l_iscb_start:l_iscb_end]
                        left_vid = left_mm.get(left_iscb)

                    if right_vid is not None or left_vid is not None:
                        if best is None:
                            best = (left_vid, right_vid)
                        break

            if best is not None:
                break

        if best is None:
            n_no_match += 1
            continue

        lv, rv = best
        if lv is not None and rv is not None:
            if lv == rv:
                concordant[lv] = concordant.get(lv, 0) + 1
                n_concordant += 1
            else:
                n_conflict += 1
        elif lv is not None:
            left_only[lv] = left_only.get(lv, 0) + 1
            n_left_only += 1
        else:
            right_only[rv] = right_only.get(rv, 0) + 1
            n_right_only += 1

    all_variants = set(left_barcodes.keys())
    total_assigned = (
        sum(concordant.values()) + sum(left_only.values()) + sum(right_only.values())
    ) or 1

    rows = []
    for vid in all_variants:
        count = concordant.get(vid, 0) + left_only.get(vid, 0) + right_only.get(vid, 0)
        rows.append(
            {
                "variant": vid,
                "mapped_reads": count,
                "concordant": concordant.get(vid, 0),
                "left_only": left_only.get(vid, 0),
                "right_only": right_only.get(vid, 0),
                "frequency_pct": round(count / total_assigned * 100, 4),
                "status": "Matched" if count > 0 else "No confident match",
            }
        )

    n_assigned = n_concordant + n_left_only + n_right_only
    qc = {
        "n_reads": n_reads,
        "concordant": n_concordant,
        "left_only": n_left_only,
        "right_only": n_right_only,
        "conflicting": n_conflict,
        "no_match": n_no_match,
        "assigned": n_assigned,
        "pct_assigned": round(100 * n_assigned / n_reads, 1) if n_reads else 0,
    }
    return rows, qc
