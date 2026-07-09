# FASTQ ↔ FASTA Matcher

A small tool that answers one question: **given the sequences you expected
and what you actually sequenced (FASTQ), which references show up in your
reads, and which don't?**

Upload your reference sequences (FASTA, GenBank, or SnapGene — one or many
expected sequences/variants, mixed formats are fine) and one or more FASTQ
files, and get a per-reference and per-read match report.

**Supported reference formats:** FASTA (`.fasta`/`.fa`/`.fna`), GenBank flat
files (`.gb`/`.gbk`/`.genbank`), and SnapGene (`.dna`). GenBank/SnapGene
records that don't carry a usable name fall back to the file's own name.

There are two versions of this tool in the repo:

| | Browser version | Python/Streamlit version |
|---|---|---|
| **Run it** | **[Live demo →](https://mitcicconi.github.io/fastq-fasta-matcher/)** | `streamlit run app.py` |
| Hosting | Static, GitHub Pages, zero server | Needs a Python process (local or Streamlit Cloud) |
| Alignment engine | Hand-written JS local aligner (approximate) | Real [minimap2](https://github.com/lh3/minimap2) via [`mappy`](https://pypi.org/project/mappy/) |
| Where your data goes | Never leaves your browser | Sent to the Streamlit process (local or cloud) |

Use the browser version for a quick check with no install. Use the Python
version when you want results computed by the actual minimap2 aligner
(same engine used by most long-read sequencing pipelines).

## Methods

Both versions run two independent matching methods side by side so you can
cross-check:

1. **MAPQ-thresholded alignment** — aligns every read against your
   reference set and keeps only reads whose best-alignment confidence is at
   or above a threshold you set (default 10).
   - *Python version:* real minimap2 MAPQ score.
   - *Browser version:* a from-scratch vanilla-JS seed-and-extend aligner
     (k-mer seed → banded Needleman-Wunsch) producing a 0–60 confidence
     score analogous to MAPQ — high when the best hit clearly beats the
     runner-up, low when hits are ambiguous. This exists so the whole app
     can run with zero server and zero external dependencies; it is **not**
     minimap2 and won't reproduce its exact numbers.
2. **K-mer containment matching** — breaks reads and references into
   overlapping k-mers and scores what fraction of a read's k-mers are
   contained in a given reference, checking both strands. Identical logic
   in both versions. Needs no external aligner and catches cases (short
   reads, high error rate, highly similar references) where alignment MAPQ
   alone is unreliable.

A reference confirmed by both methods is a high-confidence match; a
reference with zero matching reads by either method is flagged as
**no confident match**.

Try either version immediately with the bundled synthetic example in
`sample_data/` (`expected_references.fasta` + `sequencing_results.fastq`) —
two of the four references have matching reads, two don't, so you can see
both outcomes.

## Browser version

Nothing to install — open **[the live page](https://mitcicconi.github.io/fastq-fasta-matcher/)**,
upload your files, click Run. Everything (parsing, matching, alignment)
happens in your browser via plain JavaScript; no data is uploaded anywhere.
It's `index.html`, one self-contained file, served straight off the `main`
branch by GitHub Pages.

## Python/Streamlit version

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

Then open the URL Streamlit prints (usually http://localhost:8501).

### Deploying the Python version

This app has no external binary dependencies (`mappy` ships minimap2 as a
compiled Python extension, and GenBank/SnapGene parsing is pure-Python via
Biopython), so it deploys as-is on
[Streamlit Community Cloud](https://streamlit.io/cloud): point it at this
repo's `app.py` and it builds from `requirements.txt` directly.

## Project layout

```
index.html               Browser version — self-contained, served via GitHub Pages
app.py                   Streamlit UI (Python version)
matcher/
  alignment.py            MAPQ alignment matching (mappy / real minimap2)
  kmer.py                 K-mer containment matching
  io_utils.py              File upload / FASTA / GenBank / SnapGene / FASTQ helpers
sample_data/              Tiny synthetic example for demoing either version
scripts/generate_sample_data.py   Regenerates sample_data/
```
