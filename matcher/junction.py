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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

POOLS = ["T1", "T2", "T3", "T4", "T5", "T6"]
CONSTRUCTS = {"TadA": "tada", "Insertion site": "insertion_site"}
REFDB_DIR = "reference_databases"

_BASES = "ACGT"
_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def load_payload_sequence(path):
    """Reads a new payload sequence from either a FASTA file or a plain-text file."""
    if path.lower().endswith((".fasta", ".fa", ".fna")):
        recs = list(SeqIO.parse(path, "fasta"))
        if not recs:
            raise ValueError(f"No sequence found in {path}")
        return str(recs[0].seq).upper()
    with open(path) as fh:
        return "".join(line.strip() for line in fh if not line.startswith(">")).upper()


def _detect_shared_payload_len(seqs, flank):
    """Finds the payload length shared by every sequence in `seqs`, starting
    right after `flank`bp. The payload is identical across every record in a
    single-pool junction file (same construct inserted at every position),
    while the flanking IscB sequence differs by variant — so the payload is
    exactly the longest prefix (after the left flank) shared by every record.
    Robust to the right flank not being the same length as `flank` (only
    true when the payload happens to be perfectly centered — it often isn't,
    e.g. when a few bp of what looks like flank is actually still payload).
    """
    min_len = min(len(s) for s in seqs)
    if flank >= min_len:
        raise ValueError(f"Flank ({flank}bp) is longer than the shortest record ({min_len}bp)")

    after_flank = [s[flank:] for s in seqs]
    payload_len = min(len(s) for s in after_flank)
    reference = after_flank[0]
    for s in after_flank[1:]:
        limit = min(payload_len, len(s))
        i = 0
        while i < limit and reference[i] == s[i]:
            i += 1
        payload_len = min(payload_len, i)

    if payload_len <= 0:
        raise ValueError(
            "Could not detect a shared payload — is this really a single pool's "
            "junction file, all with the same construct inserted?"
        )
    return payload_len


def swap_payload(in_fasta, new_payload, out_fasta, flank=150, old_payload_len=None):
    """Rebuilds a junction reference FASTA with a new payload, keeping the flanks
    (and therefore each variant's position/codon mapping) identical.

    old_payload_len defaults to auto-detected (the prefix after `flank`bp
    shared by every record — see _detect_shared_payload_len). Pass it
    explicitly to override that.

    Returns (n_records, old_payload_len, new_payload_len).
    """
    new_payload = new_payload.upper()
    recs = list(SeqIO.parse(in_fasta, "fasta"))
    if not recs:
        raise ValueError(f"No sequences found in {in_fasta}")

    seqs = [str(rec.seq).upper() for rec in recs]
    ref_len = len(seqs[0])
    if any(len(s) != ref_len for s in seqs):
        raise ValueError(f"Not every record in {in_fasta} has the same length — is this really a single pool's junction file?")

    if old_payload_len is None:
        old_payload_len = _detect_shared_payload_len(seqs, flank)
    if old_payload_len <= 0:
        raise ValueError(f"Flank ({flank}bp x2) is longer than the record itself")
    if flank + old_payload_len > ref_len:
        raise ValueError(
            f"flank + old_payload_len ({flank + old_payload_len}bp) exceeds the record length ({ref_len}bp)"
        )

    out_recs = []
    for rec, seq in zip(recs, seqs):
        left = seq[:flank]
        right = seq[flank + old_payload_len :]
        new_seq = left + new_payload + right
        out_recs.append(SeqRecord(Seq(new_seq), id=rec.id, description=""))

    SeqIO.write(out_recs, out_fasta, "fasta")
    return len(out_recs), old_payload_len, len(new_payload)


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
    """Loads a [left flank][payload][right flank] FASTA and derives per-variant
    barcodes. The left flank is exactly `flank`bp; the right flank is whatever's
    left after the payload — NOT assumed to also be `flank`bp, since that's only
    true when the payload happens to be perfectly centered.

    Payload length is auto-detected rather than assumed: the payload is
    identical across every record in a single-pool junction file (same
    construct inserted at every position), while the flanking IscB sequence
    differs by variant. So the payload is exactly the longest prefix (starting
    right after the left flank) shared by every record — this finds the
    payload/right-flank boundary without needing to know the payload's length
    or content in advance, and works whether or not the right flank happens to
    equal `flank`bp.

    Returns (left_barcodes, right_barcodes, payload_seq, payload_len) where the
    two barcode dicts map variant id -> 20bp barcode sequence.
    """
    recs = list(SeqIO.parse(fasta_path, "fasta"))
    if not recs:
        raise ValueError(f"No sequences found in {fasta_path}")

    seqs = [str(r.seq).upper() for r in recs]
    try:
        payload_len = _detect_shared_payload_len(seqs, flank)
    except ValueError as exc:
        raise ValueError(f"{exc} ({fasta_path})") from exc

    left_barcodes, right_barcodes = {}, {}
    payload_seq = seqs[0][flank : flank + payload_len]
    for rec, seq in zip(recs, seqs):
        left_barcodes[rec.id] = seq[flank - bc_len : flank]
        right_barcodes[rec.id] = seq[flank + payload_len : flank + payload_len + bc_len]

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


def detect_best_pool(fastq_records, construct_suffix, refdb_dir, pools=POOLS, sample_size=500, mm_allow=1):
    """Tests a sample of reads against every pool's database for a given
    construct, to catch a pool/file mismatch (a FASTQ file's name/label not
    actually matching the pool its reads came from) before running the full
    match. Returns (pool, pct_assigned, variants_matched, n_variants) rows,
    sorted best-first.
    """
    import os

    sample = fastq_records[:sample_size]
    results = []
    for pool in pools:
        path = os.path.join(refdb_dir, f"{pool}_{construct_suffix}.fasta")
        if not os.path.exists(path):
            continue
        left_bc, right_bc, payload_seq, payload_len = load_junction_reference(path)
        rows, qc = match_reads_by_junction(sample, left_bc, right_bc, payload_seq, payload_len, mm_allow=mm_allow)
        matched = sum(1 for r in rows if r["mapped_reads"] > 0)
        results.append({"pool": pool, "pct_assigned": qc["pct_assigned"], "variants_matched": matched, "n_variants": len(rows)})
    results.sort(key=lambda r: r["pct_assigned"], reverse=True)
    return results
