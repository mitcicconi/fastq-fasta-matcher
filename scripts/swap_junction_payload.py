#!/usr/bin/env python3
"""Regenerates a junction reference database with the payload swapped.

Takes an existing junction FASTA (built as [flank][old payload][flank] per
record — e.g. a TadA junction database from *pipeline_junction.py) and a new
payload sequence, and writes out a new junction FASTA with
[flank][new payload][flank] for every record. The flanks — and therefore
each variant's position/codon mapping — are preserved exactly; only the
payload in the middle changes.

CLI wrapper over matcher.junction.swap_payload — for a UI instead, run
`streamlit run app.py` and use the "Build Junction Database" page.

Usage:
    python3 scripts/swap_junction_payload.py T1_junction.fasta insertion_site.fasta T1_junction_insertion_site.fasta
    python3 scripts/swap_junction_payload.py T1_junction.fasta insertion_site.fasta T1_junction_insertion_site.fasta --flank 150
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from matcher.junction import load_payload_sequence, swap_payload

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
    parser.add_argument(
        "--old-payload-len",
        type=int,
        default=None,
        help="Override the old payload's length instead of deriving it from --flank. Use this when the old "
        "payload actually extends further than the nominal flank boundary (e.g. a few bp of linker got bucketed "
        "into what looked like flank) — the right flank then comes out shorter than --flank, whatever's left.",
    )
    args = parser.parse_args()

    payload = load_payload_sequence(args.new_payload)
    n, old_len, new_len = swap_payload(
        args.in_fasta, payload, args.out_fasta, flank=args.flank, old_payload_len=args.old_payload_len
    )
    print(f"Wrote {n} records to {args.out_fasta}")
    print(f"  old payload length: {old_len}bp -> new payload length: {new_len}bp")
