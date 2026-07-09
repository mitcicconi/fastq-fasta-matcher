# FASTQ ↔ FASTA Matcher

A small Streamlit app that answers one question: **given the sequences you
expected (FASTA) and what you actually sequenced (FASTQ), which references
show up in your reads, and which don't?**

Upload a reference FASTA (one or many expected sequences/variants) and one
or more FASTQ files, and get a per-reference and per-read match report —
no command line, no local aligner install required.

## Methods

Two independent matching methods are run side by side so you can cross-check:

1. **MAPQ-thresholded alignment** — aligns every read against your
   reference set with [minimap2](https://github.com/lh3/minimap2) (via its
   Python binding, [`mappy`](https://pypi.org/project/mappy/)) and keeps
   only reads whose best-alignment MAPQ is at or above a threshold you set
   (default 10).
2. **K-mer containment matching** — breaks reads and references into
   overlapping k-mers and scores what fraction of a read's k-mers are
   contained in a given reference, checking both strands. This needs no
   external aligner and catches cases (short reads, high error rate,
   highly similar references) where alignment MAPQ alone is unreliable.

A reference confirmed by both methods is a high-confidence match; a
reference with zero matching reads by either method is flagged as
**no confident match**.

## Running locally

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

Then open the URL Streamlit prints (usually http://localhost:8501).

Try it immediately with the bundled synthetic example in `sample_data/`
(`expected_references.fasta` + `sequencing_results.fastq`) — two of the
four references have matching reads, two don't, so you can see both
outcomes.

## Deploying

This app has no external binary dependencies (`mappy` ships minimap2 as a
compiled Python extension), so it deploys as-is on
[Streamlit Community Cloud](https://streamlit.io/cloud): point it at this
repo's `app.py` and it builds from `requirements.txt` directly.

## Project layout

```
app.py                  Streamlit UI
matcher/
  alignment.py           MAPQ alignment matching (mappy)
  kmer.py                K-mer containment matching
  io_utils.py             File upload / FASTA / FASTQ helpers
sample_data/             Tiny synthetic example for demoing the app
scripts/generate_sample_data.py   Regenerates sample_data/
```
