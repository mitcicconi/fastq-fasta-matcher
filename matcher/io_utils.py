"""Helpers for turning Streamlit-uploaded files into on-disk paths / in-memory records."""

import io
import os
import tempfile

import mappy as mp
from Bio import SeqIO

_FORMAT_BY_EXT = {
    "dna": "snapgene",
    "gb": "genbank",
    "gbk": "genbank",
    "genbank": "genbank",
}
_GENERIC_IDS = {"", ".", "none", "<unknown id>", "unknown", "unnamed", "exported"}


def save_uploaded_files(uploaded_files, suffix):
    """Write one or more Streamlit UploadedFile objects to a single concatenated temp file.

    FASTA/FASTQ are line-oriented, so concatenating multiple uploads into one file
    is equivalent to treating them as one combined input.
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix, mode="wb")
    for uf in uploaded_files:
        tmp.write(uf.getbuffer())
        if not uf.getbuffer().tobytes().endswith(b"\n"):
            tmp.write(b"\n")
    tmp.close()
    return tmp.name


def normalize_references_to_fasta(uploaded_files):
    """Accepts reference files in FASTA, GenBank (.gb/.gbk/.genbank), or SnapGene
    (.dna) format — mixed together, if you like — and writes them all out as one
    combined FASTA file for the matchers to use. GenBank and SnapGene records often
    carry a generic/empty id, so those fall back to the file's name.
    """
    records = []
    used_ids = set()

    def unique_id(candidate):
        candidate = candidate or "sequence"
        name, i = candidate, 2
        while name in used_ids:
            name = f"{candidate}_{i}"
            i += 1
        used_ids.add(name)
        return name

    for uf in uploaded_files:
        ext = uf.name.rsplit(".", 1)[-1].lower() if "." in uf.name else ""
        fmt = _FORMAT_BY_EXT.get(ext, "fasta")
        stem = os.path.splitext(uf.name)[0]
        raw = uf.getvalue()
        handle = io.BytesIO(raw) if fmt == "snapgene" else io.StringIO(raw.decode("utf-8", errors="replace"))

        try:
            parsed = list(SeqIO.parse(handle, fmt))
        except Exception as exc:
            raise ValueError(f"Could not parse {uf.name} as {fmt}: {exc}") from exc
        if not parsed:
            raise ValueError(f"No sequences found in {uf.name}")

        for i, rec in enumerate(parsed):
            base_id = stem if rec.id.strip().lower() in _GENERIC_IDS else rec.id
            if len(parsed) > 1 and rec.id.strip().lower() in _GENERIC_IDS:
                base_id = f"{stem}_{i + 1}"
            rec.id = unique_id(base_id)
            rec.description = ""
            records.append(rec)

    fd, path = tempfile.mkstemp(suffix=".fasta")
    os.close(fd)
    SeqIO.write(records, path, "fasta")
    return path


def read_fasta_records(path):
    """Return list of (name, sequence) tuples from a FASTA file."""
    return [(name, seq) for name, seq, _ in mp.fastx_read(path)]


def read_fastq_records(path):
    """Return list of (name, sequence, quality) tuples from a FASTQ file."""
    return list(mp.fastx_read(path))


def cleanup(*paths):
    for p in paths:
        try:
            os.unlink(p)
        except OSError:
            pass
