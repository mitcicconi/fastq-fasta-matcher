"""Helpers for turning Streamlit-uploaded files into on-disk paths / in-memory records."""

import os
import tempfile

import mappy as mp


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
