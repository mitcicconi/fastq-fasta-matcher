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
payload in the middle (TadA+linkers, 618bp → the insertion-site placeholder,
24bp).

Each file is a "junction reference": every record is
`[150bp left flank][payload][right flank]`. The payload length is
auto-detected from the file itself (`matcher.junction.load_junction_reference`
and `swap_payload` both do this — see the "618, not 614" note below for why
this matters), not assumed fixed, so the right flank isn't necessarily also
150bp.

## A real bug this project hit, worth knowing about

The payload actually replaced when swapping TadA for the insertion site is
**618bp, not 614bp**. The original `*pipeline_kmer.py`'s `TADA_LINKERS`
constant (614bp) is missing 4bp (`TGGA`) that are genuinely part of the
TadA+linkers block being replaced — those 4bp were silently absorbed into
what the original pipeline treated as "150bp right flank" (harmless there,
since they're identical across all variants). But naively reusing that
614bp boundary when building the insertion-site database left those 4bp
stuck onto the wrong side of the junction, causing a 4bp frameshift in
every single barcode and a 0% match rate.

Fixed by making `swap_payload`/`load_junction_reference` auto-detect the
payload boundary instead of assuming a fixed length: since the payload is
identical across every variant in a pool (same construct inserted at every
position) while the flanking IscB sequence differs by variant, the payload
is exactly the longest prefix (after the left flank) shared by every
record in the file. This finds the true boundary without needing to know
the payload's exact length in advance, and self-corrects this exact class
of error.

## Regenerating

If the insertion-site sequence changes, rebuild with:

```bash
python3 scripts/swap_junction_payload.py \
    reference_databases/T1_tada.fasta \
    /path/to/new_insertion_site.fasta \
    reference_databases/T1_insertion_site.fasta
```

Repeat for T2–T6. (No `--flank`/`--old-payload-len` flags needed — the
payload boundary is auto-detected. Only pass `--old-payload-len` if you
need to force a specific boundary for some reason.)

## If a FASTQ file gets 0% matched against every pool

Don't assume the file's name tells you which pool it's from. This project
hit a real case where a new sequencing run's `Var1`–`Var6` file labels
turned out to be numbered in the *opposite* order from the original
`T1`–`T6` pools (`Var1`→`T6`, `Var2`→`T5`, ... `Var6`→`T1`) — completely
consistent, just reversed, for reasons not documented at the point of
handoff. Use the **"Detect which pool this file belongs to"** button in the
app (samples the first 500 reads against every pool's database for the
selected construct and reports the best match) before concluding a file's
data is bad.
