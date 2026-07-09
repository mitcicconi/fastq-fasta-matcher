"""FASTQ ↔ FASTA Matcher
A user-friendly UI for checking which sequencing reads match which reference
sequences, using two independent methods:

  1. MAPQ-thresholded alignment (minimap2 via the mappy binding)
  2. K-mer containment matching (orientation-independent, no external tools)

Upload your reference FASTA(s) and read FASTQ(s), pick your thresholds, and
get a per-reference and per-read match/no-match report you can download.
"""

import pandas as pd
import streamlit as st

from matcher.alignment import run_alignment_match
from matcher.io_utils import cleanup, normalize_references_to_fasta, read_fastq_records, save_uploaded_files
from matcher.kmer import build_reference_kmer_index, match_reads_by_kmer, summarize_by_reference

st.set_page_config(page_title="FASTQ ↔ FASTA Matcher", layout="wide")

st.title("FASTQ ↔ FASTA Matcher")
st.caption(
    "Upload the reference sequences you expect (FASTA, GenBank, or SnapGene) and "
    "the FASTQ reads you sequenced. This tool tells you which references were "
    "actually found in your reads, and which reads didn't confidently match anything."
)

with st.sidebar:
    st.header("Inputs")
    fasta_files = st.file_uploader(
        "Reference FASTA / GenBank / SnapGene — the sequences you expect",
        type=["fasta", "fa", "fna", "gb", "gbk", "genbank", "dna"],
        accept_multiple_files=True,
    )
    fastq_files = st.file_uploader(
        "Read FASTQ(s) — your sequencing results",
        type=["fastq", "fq"],
        accept_multiple_files=True,
    )

    st.header("Alignment method (MAPQ)")
    run_alignment = st.checkbox("Run MAPQ alignment matching", value=True)
    mapq_threshold = st.slider("MAPQ threshold", min_value=0, max_value=60, value=10)
    preset = st.selectbox(
        "minimap2 preset",
        options=["map-ont", "map-pb", "map-hifi", "sr", "asm5", "asm10", "asm20"],
        index=0,
        help="map-ont = Nanopore reads, sr = short accurate reads (Illumina), "
        "asm* = assembly-to-reference comparison.",
    )

    st.header("K-mer method")
    run_kmer = st.checkbox("Run k-mer containment matching", value=True)
    kmer_k = st.slider("K-mer length", min_value=11, max_value=41, value=17, step=2)
    min_containment = st.slider(
        "Minimum containment to call a match",
        min_value=0.02,
        max_value=1.0,
        value=0.1,
        step=0.02,
        help="Shared k-mers as a fraction of whichever is smaller — the read's "
        "k-mer count or the reference's. Normalizing by the smaller side keeps "
        "this meaningful even when a read is much longer than the reference "
        "(e.g. a long amplicon read with only a short target region).",
    )

    go = st.button("Run matching", type="primary", use_container_width=True)

if not go:
    st.info("Upload files on the left, choose your settings, then click **Run matching**.")
    st.markdown(
        "**How it works**\n\n"
        "- *MAPQ alignment* aligns each read to your reference set with minimap2 "
        "and keeps only reads whose best alignment quality (MAPQ) is at or above "
        "your threshold — the same approach as this project's `*pipeline.py`.\n"
        "- *K-mer containment* breaks reads and references into short k-mers and "
        "scores how much of each read's k-mer content is contained in a given "
        "reference, checking both strands — a generalization of the barcode "
        "matching in `*pipeline_kmer.py` that works for any sequences, not just "
        "one specific assay.\n"
        "- Running both gives you a cross-check: a reference confirmed by both "
        "methods is a high-confidence match."
    )
    st.stop()

if not fasta_files:
    st.error("Upload at least one reference file (FASTA, GenBank, or SnapGene).")
    st.stop()
if not fastq_files:
    st.error("Upload at least one FASTQ read file.")
    st.stop()
if not run_alignment and not run_kmer:
    st.error("Enable at least one matching method in the sidebar.")
    st.stop()

try:
    fasta_path = normalize_references_to_fasta(fasta_files)
except ValueError as exc:
    st.error(str(exc))
    st.stop()
fastq_path = save_uploaded_files(fastq_files, ".fastq")

try:
    with st.spinner("Reading input files..."):
        from matcher.io_utils import read_fasta_records

        fasta_records = read_fasta_records(fasta_path)
        ref_names = [name for name, _ in fasta_records]
        fastq_records = read_fastq_records(fastq_path)

    st.success(
        f"Loaded {len(ref_names)} reference sequence(s) and {len(fastq_records)} read(s)."
    )

    align_per_ref = align_per_read = None
    kmer_per_ref = kmer_per_read = None

    col_a, col_b = st.columns(2)

    if run_alignment:
        with st.spinner("Running MAPQ alignment matching..."):
            align_per_ref, align_per_read, total_reads = run_alignment_match(
                fasta_path, fastq_path, preset=preset, mapq_threshold=mapq_threshold
            )
        pass_col = f"reads_mapq_ge_{mapq_threshold}"
        n_matched_refs = (align_per_ref[pass_col] > 0).sum()
        n_matched_reads = int(align_per_read["matched"].sum())
        with col_a:
            st.subheader("MAPQ alignment results")
            st.metric(
                f"References matched (MAPQ ≥ {mapq_threshold})",
                f"{n_matched_refs} / {len(ref_names)}",
            )
            st.metric(
                "Reads matched",
                f"{n_matched_reads} / {total_reads}"
                f"  ({100 * n_matched_reads / total_reads:.1f}%)"
                if total_reads
                else "0 / 0",
            )

    if run_kmer:
        with st.spinner("Running k-mer containment matching..."):
            ref_kmers = build_reference_kmer_index(fasta_records, kmer_k)
            kmer_per_read = match_reads_by_kmer(
                fastq_records, ref_kmers, kmer_k, min_containment=min_containment
            )
            kmer_per_ref = summarize_by_reference(kmer_per_read, ref_names)
        n_matched_refs_k = (kmer_per_ref["reads_matched_kmer"] > 0).sum()
        n_matched_reads_k = int(kmer_per_read["matched"].sum())
        total_reads_k = len(kmer_per_read)
        with col_b:
            st.subheader("K-mer matching results")
            st.metric(
                f"References matched (containment ≥ {min_containment})",
                f"{n_matched_refs_k} / {len(ref_names)}",
            )
            st.metric(
                "Reads matched",
                f"{n_matched_reads_k} / {total_reads_k}"
                f"  ({100 * n_matched_reads_k / total_reads_k:.1f}%)"
                if total_reads_k
                else "0 / 0",
            )

    st.divider()
    st.header("Per-reference report")

    if run_alignment and run_kmer:
        combined = align_per_ref.merge(
            kmer_per_ref, on="reference", suffixes=("_mapq", "_kmer")
        )
        combined["status"] = combined.apply(
            lambda r: "Matched (both methods)"
            if r["status_mapq"] == "Matched" and r["status_kmer"] == "Matched"
            else (
                "Matched (one method only)"
                if r["status_mapq"] == "Matched" or r["status_kmer"] == "Matched"
                else "No confident match"
            ),
            axis=1,
        )
        combined = combined.drop(columns=["status_mapq", "status_kmer"])
        report_df = combined
    elif run_alignment:
        report_df = align_per_ref
    else:
        report_df = kmer_per_ref

    status_filter = st.multiselect(
        "Filter by status",
        options=sorted(report_df["status"].unique()),
        default=list(report_df["status"].unique()),
    )
    filtered = report_df[report_df["status"].isin(status_filter)]
    st.dataframe(filtered, use_container_width=True, hide_index=True)

    unmatched_refs = report_df[report_df["status"] == "No confident match"]
    if len(unmatched_refs):
        st.warning(
            f"{len(unmatched_refs)} reference sequence(s) had no confidently "
            f"matching reads: {', '.join(unmatched_refs['reference'].astype(str).tolist())}"
        )
    else:
        st.success("Every reference sequence had at least one confidently matching read.")

    numeric_cols = [c for c in report_df.columns if "reads" in c]
    if numeric_cols:
        chart_df = report_df.set_index("reference")[numeric_cols].sort_values(
            numeric_cols[0], ascending=False
        )
        st.bar_chart(chart_df)

    st.download_button(
        "Download per-reference report (CSV)",
        report_df.to_csv(index=False).encode("utf-8"),
        file_name="per_reference_report.csv",
        mime="text/csv",
    )

    st.divider()
    st.header("Per-read report")
    tabs = st.tabs(
        [name for name, enabled in [("MAPQ alignment", run_alignment), ("K-mer", run_kmer)] if enabled]
    )
    tab_idx = 0
    if run_alignment:
        with tabs[tab_idx]:
            st.dataframe(align_per_read, use_container_width=True, hide_index=True)
            unmatched = align_per_read[~align_per_read["matched"]]
            st.caption(f"{len(unmatched)} of {len(align_per_read)} reads did not pass MAPQ threshold.")
            st.download_button(
                "Download per-read MAPQ results (CSV)",
                align_per_read.to_csv(index=False).encode("utf-8"),
                file_name="per_read_mapq_results.csv",
                mime="text/csv",
                key="dl_mapq_reads",
            )
        tab_idx += 1
    if run_kmer:
        with tabs[tab_idx]:
            st.dataframe(kmer_per_read, use_container_width=True, hide_index=True)
            unmatched = kmer_per_read[~kmer_per_read["matched"]]
            st.caption(f"{len(unmatched)} of {len(kmer_per_read)} reads did not pass the containment threshold.")
            st.download_button(
                "Download per-read k-mer results (CSV)",
                kmer_per_read.to_csv(index=False).encode("utf-8"),
                file_name="per_read_kmer_results.csv",
                mime="text/csv",
                key="dl_kmer_reads",
            )

finally:
    cleanup(fasta_path, fastq_path)
