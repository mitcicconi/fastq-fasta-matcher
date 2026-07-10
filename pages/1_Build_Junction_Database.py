"""Build Junction Database
A small UI over matcher.junction.swap_payload — the same operation
scripts/swap_junction_payload.py does from the command line, but with file
pickers and immediate feedback instead of typing out paths. Still fully
local: this is just another page of the same `streamlit run app.py` process.
"""

import os
import re
import tempfile

import streamlit as st

from matcher.io_utils import cleanup
from matcher.junction import CONSTRUCTS, POOLS, REFDB_DIR, load_payload_sequence, swap_payload

st.set_page_config(page_title="Build Junction Database", layout="wide")

st.title("Build Junction Database")
st.caption(
    "Derive a new junction reference database from an existing one by swapping "
    "the payload in the middle — the flanks (and therefore each variant's "
    "position/codon mapping) stay identical. Same operation as "
    "`scripts/swap_junction_payload.py`, with a UI."
)

st.header("1. Source junction FASTA")
source_mode = st.radio(
    "Where's the existing junction reference?",
    ["Existing pool in reference_databases/", "Upload a junction FASTA"],
    horizontal=True,
)

in_fasta_path = None
default_out_name = None

if source_mode == "Existing pool in reference_databases/":
    col1, col2 = st.columns(2)
    with col1:
        pool = st.selectbox("Pool", POOLS)
    with col2:
        source_construct = st.selectbox("Source construct", list(CONSTRUCTS.keys()), index=0)
    candidate_path = os.path.join(REFDB_DIR, f"{pool}_{CONSTRUCTS[source_construct]}.fasta")
    if os.path.exists(candidate_path):
        in_fasta_path = candidate_path
        st.success(f"Found `{candidate_path}`")
    else:
        st.error(f"No file at `{candidate_path}` — see {REFDB_DIR}/README.md")
    default_out_name = f"{pool}_"
else:
    uploaded_ref = st.file_uploader("Junction FASTA", type=["fasta", "fa", "fna"])
    if uploaded_ref is not None:
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
        tmp.write(uploaded_ref.getvalue())
        tmp.close()
        in_fasta_path = tmp.name
    default_out_name = ""

st.header("2. New payload")
payload_mode = st.radio("How are you providing the new payload?", ["Paste sequence", "Upload a file"], horizontal=True)

payload_seq = None
if payload_mode == "Paste sequence":
    pasted = st.text_area("Payload sequence (letters other than ACGT are ignored)", height=100)
    if pasted.strip():
        payload_seq = re.sub(r"[^ACGTacgt]", "", pasted).upper()
        st.caption(f"{len(payload_seq)}bp after cleanup")
else:
    uploaded_payload = st.file_uploader("Payload file (FASTA or plain text)", type=["fasta", "fa", "fna", "txt"])
    if uploaded_payload is not None:
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(uploaded_payload.name)[1])
        tmp.write(uploaded_payload.getvalue())
        tmp.close()
        try:
            payload_seq = load_payload_sequence(tmp.name)
            st.caption(f"{len(payload_seq)}bp")
        finally:
            cleanup(tmp.name)

st.header("3. Output")
flank = st.number_input(
    "Flank length on each side",
    min_value=1,
    value=150,
    help="Matches *pipeline_junction.py's FLANK setting — leave at 150 unless you know it's different.",
)

if source_mode == "Existing pool in reference_databases/" and in_fasta_path:
    construct_label = st.text_input(
        "New construct name (used for the filename and the dropdown label)",
        value="insertion_site",
        help="Output will be saved as reference_databases/{pool}_{this}.fasta",
    )
    out_filename = f"{pool}_{re.sub(r'[^a-z0-9_]', '_', construct_label.lower())}.fasta"
else:
    out_filename = st.text_input("Output filename", value="junction_database_new_payload.fasta")

build = st.button("Build database", type="primary")

if build:
    if not in_fasta_path:
        st.error("No source junction FASTA available yet.")
    elif not payload_seq:
        st.error("No new payload sequence provided yet.")
    else:
        tmp_out = None
        try:
            fd, tmp_out = tempfile.mkstemp(suffix=".fasta")
            os.close(fd)
            n, old_len, new_len = swap_payload(in_fasta_path, payload_seq, tmp_out, flank=flank)
            with open(tmp_out, "rb") as fh:
                result_bytes = fh.read()
            # st.button's True only lasts for the one rerun it triggered — stash the
            # result in session_state so it survives the rerun that clicking Save or
            # Download themselves cause, instead of vanishing before Save can run.
            st.session_state["built_result"] = {
                "bytes": result_bytes,
                "filename": out_filename,
                "stats": (n, old_len, new_len),
            }
        except ValueError as exc:
            st.error(str(exc))
            st.session_state.pop("built_result", None)
        finally:
            if tmp_out:
                cleanup(tmp_out)

result = st.session_state.get("built_result")
if result:
    n, old_len, new_len = result["stats"]
    st.success(f"Built `{result['filename']}` — {n} record(s), payload {old_len}bp → {new_len}bp")
    st.download_button(
        "Download FASTA",
        result["bytes"],
        file_name=result["filename"],
        mime="text/x-fasta",
    )

    dest_path = os.path.join(REFDB_DIR, result["filename"])
    if os.path.exists(dest_path):
        st.warning(f"`{dest_path}` already exists.")
        overwrite = st.checkbox(f"Overwrite {dest_path}")
    else:
        overwrite = True

    if st.button(f"Save to {REFDB_DIR}/ (so it appears in the matcher's dropdown)", disabled=not overwrite):
        with open(dest_path, "wb") as fh:
            fh.write(result["bytes"])
        st.success(f"Saved to `{dest_path}` — it'll show up next time you pick a pool/construct on the main page.")
        del st.session_state["built_result"]
