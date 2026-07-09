"""MAPQ-thresholded alignment matching, generalized from *pipeline.py.

Uses mappy (the minimap2 Python binding) so no external minimap2/samtools
binary is required — this is what makes the tool deployable on a hosted
server instead of only on a machine with those tools installed locally.
"""

import mappy as mp
import pandas as pd


def run_alignment_match(fasta_path, fastq_path, preset="map-ont", mapq_threshold=10):
    """Align every read in fastq_path against every reference in fasta_path.

    Returns (per_reference_df, per_read_df, total_reads).
    """
    aligner = mp.Aligner(fasta_path, preset=preset)
    if not aligner:
        raise RuntimeError(
            "Could not build an index from the reference FASTA — check the file is valid."
        )

    ref_names = list(aligner.seq_names)
    total_aligned = {r: 0 for r in ref_names}
    total_passed = {r: 0 for r in ref_names}
    read_rows = []
    total_reads = 0

    for name, seq, _qual in mp.fastx_read(fastq_path):
        total_reads += 1
        best_hit = None
        for hit in aligner.map(seq):
            if not hit.is_primary:
                continue
            if best_hit is None or hit.mapq > best_hit.mapq:
                best_hit = hit

        if best_hit is None:
            read_rows.append(
                {"read_id": name, "reference": None, "mapq": None, "matched": False}
            )
            continue

        total_aligned[best_hit.ctg] = total_aligned.get(best_hit.ctg, 0) + 1
        passed = best_hit.mapq >= mapq_threshold
        if passed:
            total_passed[best_hit.ctg] = total_passed.get(best_hit.ctg, 0) + 1

        read_rows.append(
            {
                "read_id": name,
                "reference": best_hit.ctg,
                "mapq": best_hit.mapq,
                "matched": passed,
            }
        )

    pass_col = f"reads_mapq_ge_{mapq_threshold}"
    per_reference = pd.DataFrame(
        [
            {
                "reference": r,
                "reads_aligned_any_mapq": total_aligned.get(r, 0),
                pass_col: total_passed.get(r, 0),
                "status": "Matched" if total_passed.get(r, 0) > 0 else "No confident match",
            }
            for r in ref_names
        ]
    ).sort_values(pass_col, ascending=False).reset_index(drop=True)

    per_read = pd.DataFrame(read_rows)
    return per_reference, per_read, total_reads
