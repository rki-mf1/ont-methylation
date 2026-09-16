import subprocess
from datetime import datetime
from pathlib import Path

import pandas as pd
import streamlit as st

REPO_ROOT = Path(__file__).resolve().parent.parent
RUNS_DIR = REPO_ROOT / "gui_runs"
SAMPLE_NAME = "sample"

st.set_page_config(page_title="ONT-bmod", layout="wide")
st.title("ONT-bmod")
st.caption("Bacterial methylation motif extraction from Nanopore data — main flow")

st.header("1. Input files")
fasta_file = st.file_uploader("Reference FASTA", type=["fasta", "fa"])
bam_file = st.file_uploader("Basecalled BAM (Dorado output, with MM/ML tags)", type=["bam"])

st.header("2. Parameters")
col1, col2 = st.columns(2)
with col1:
    auto_threshold = st.checkbox("Automatic modkit threshold", value=False)
    filter_threshold = st.number_input(
        "Modkit filter threshold", min_value=0.0, max_value=1.0, value=0.75, disabled=auto_threshold
    )
with col2:
    percent_cutoff = st.number_input(
        "Modification table % cutoff", min_value=0.0, max_value=1.0, value=0.5
    )

run_button = st.button("Run pipeline", type="primary", disabled=not (fasta_file and bam_file))

log_placeholder = st.empty()
results_placeholder = st.container()

if run_button:
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    # Save uploads under a shared basename: the pipeline requires --fasta
    # and --bam to have matching basenames, regardless of what the user's
    # original files were called.
    fasta_path = run_dir / f"{SAMPLE_NAME}.fasta"
    bam_path = run_dir / f"{SAMPLE_NAME}.bam"
    fasta_path.write_bytes(fasta_file.getvalue())
    bam_path.write_bytes(bam_file.getvalue())

    outdir = run_dir / "results"
    workdir = run_dir / "work"

    cmd = [
        "nextflow", "run", str(REPO_ROOT / "ont-bmod.nf"),
        "--fasta", str(fasta_path),
        "--bam", str(bam_path),
        "--percent_cutoff_modification_table", str(percent_cutoff),
        "-w", str(workdir),
        "--outdir", str(outdir),
    ]
    if auto_threshold:
        cmd += ["--automatic_threshold_modkit", "true"]
    else:
        cmd += ["--filter_threshold_modkit", str(filter_threshold)]

    log_lines = []
    with st.spinner("Running..."):
        process = subprocess.Popen(
            cmd, cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        for line in process.stdout:
            log_lines.append(line.rstrip())
            log_placeholder.code("\n".join(log_lines[-40:]))
        process.wait()

    if process.returncode == 0:
        st.success("Pipeline completed successfully.")
    else:
        st.error(f"Pipeline failed (exit code {process.returncode}). See log above.")

    with results_placeholder:
        st.header("3. Results")
        sample_outdir = outdir / SAMPLE_NAME
        if sample_outdir.exists():
            stats_dir = sample_outdir / "methylation_statistics"
            if stats_dir.exists():
                st.subheader("Methylation statistics")
                for csv_file in sorted(stats_dir.glob("*.csv")):
                    df = pd.read_csv(csv_file, sep="\t")
                    st.write(csv_file.stem.replace("_", " "))
                    st.dataframe(df)

            motifs_file = sample_outdir / "modkit_motifs.tsv"
            if motifs_file.exists():
                st.subheader("Motifs")
                st.dataframe(pd.read_csv(motifs_file, sep="\t"))

            st.subheader("All output files")
            for f in sorted(sample_outdir.rglob("*")):
                if f.is_file():
                    st.text(str(f.relative_to(sample_outdir)))
        else:
            st.warning("No output directory found.")
