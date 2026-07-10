# Bundled junction-barcode reference databases

This is where the app looks for pre-built junction-barcode reference FASTAs
to power the "Bundled junction database" mode in `app.py` — an alternative
to uploading a FASTA every time, where you instead pick a pool + construct
from a dropdown.

## Files

One FASTA per pool x construct, named exactly:

```
T1_tada.fasta            T1_insertion_site.fasta
T2_tada.fasta            T2_insertion_site.fasta
T3_tada.fasta            T3_insertion_site.fasta
T4_tada.fasta            T4_insertion_site.fasta
T5_tada.fasta            T5_insertion_site.fasta
T6_tada.fasta            T6_insertion_site.fasta
```

The `*_tada.fasta` files are direct copies of `*pipeline_junction.py`'s
output (`*Junction Intermediate/{pool}_junction.fasta` in the main project
directory). The `*_insertion_site.fasta` files were derived from those with
`scripts/swap_junction_payload.py`, which keeps the flanks — and therefore
each variant's position/codon mapping — identical, and only swaps the
payload in the middle (TadA, 614bp → the insertion-site placeholder, 24bp).

Each file is a "junction reference": every record is
`[150bp flank][payload][150bp flank]`, with every record in a file sharing
the same payload length.

## Regenerating

If the insertion-site sequence changes, rebuild with:

```bash
python3 scripts/swap_junction_payload.py \
    reference_databases/T1_tada.fasta \
    /path/to/new_insertion_site.fasta \
    reference_databases/T1_insertion_site.fasta
```

Repeat for T2–T6.
