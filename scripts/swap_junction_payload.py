#!/usr/bin/env python3
"""Regenerates a junction reference database with the payload swapped.

Takes an existing junction FASTA (built as [flank][old payload][flank] per
record — e.g. a TadA junction database from *pipeline_junction.py) and a new
payload sequence, and writes out a new junction FASTA with
[flank][new payload][flank] for every record. The flanks — and therefore
each variant's position/codon mapping — are preserved exactly; only the
payload in the middle changes.

Usage:
    python3 scripts/swap_junction_payload.py T1_junction.fasta insertion_site.fasta T1_junction_insertion_site.fasta
    python3 scripts/swap_junction_payload.py T1_junction.fasta insertion_site.fasta T1_junction_insertion_site.fasta --flank 150
"""

import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_payload_sequence(path):
    if path.lower().endswith((".fasta", ".fa", ".fna")):
        recs = list(SeqIO.parse(path, "fasta"))
        if not recs:
            raise ValueError(f"No sequence found in {path}")
        return str(recs[0].seq).upper()
    with open(path) as fh:
        return "".join(line.strip() for line in fh if not line.startswith(">")).upper()


def swap_payload(in_fasta, new_payload, out_fasta, flank=150):
    new_payload = new_payload.upper()
    recs = list(SeqIO.parse(in_fasta, "fasta"))
    if not recs:
        raise ValueError(f"No sequences found in {in_fasta}")

    old_payload_len = len(recs[0].seq) - 2 * flank
    if old_payload_len <= 0:
        raise ValueError(f"Flank ({flank}bp x2) is longer than the record itself")

    out_recs = []
    for rec in recs:
        seq = str(rec.seq).upper()
        if len(seq) - 2 * flank != old_payload_len:
            raise ValueError(
                f"Record {rec.id} has a different length than the rest of {in_fasta} "
                "— is this really a single pool's junction file?"
            )
        left = seq[:flank]
        right = seq[len(seq) - flank :]
        new_seq = left + new_payload + right
        out_recs.append(SeqRecord(Seq(new_seq), id=rec.id, description=""))

    SeqIO.write(out_recs, out_fasta, "fasta")
    return len(out_recs), old_payload_len, len(new_payload)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_fasta", help="Existing junction reference FASTA (e.g. the TadA version)")
    parser.add_argument(
        "new_payload", help="New payload sequence: a FASTA file, or a plain-text file with just the sequence"
    )
    parser.add_argument("out_fasta", help="Output path for the new junction reference FASTA")
    parser.add_argument(
        "--flank", type=int, default=150, help="Flank length on each side (default 150, matches *pipeline_junction.py)"
    )
    args = parser.parse_args()

    payload = load_payload_sequence(args.new_payload)
    n, old_len, new_len = swap_payload(args.in_fasta, payload, args.out_fasta, flank=args.flank)
    print(f"Wrote {n} records to {args.out_fasta}")
    print(f"  old payload length: {old_len}bp -> new payload length: {new_len}bp")
