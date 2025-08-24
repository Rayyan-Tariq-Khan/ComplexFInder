import streamlit as st
import os
import time
import pandas as pd
from analysis import run_complex_finder
from corum_analysis import run_corum_analysis
from stringdb_analysis import generate_stringdb_url
from plotting import plot_protein_profiles
from Complexes import run_complex_analysis
import shutil

# -------------------------
# Initialize Session State
# -------------------------
if 'analysis_done' not in st.session_state:
    st.session_state.analysis_done = False
if 'demo_mode' not in st.session_state:
    st.session_state.demo_mode = False
if 'result_files' not in st.session_state:
    st.session_state.result_files = {}
if 'corum_output' not in st.session_state:
    st.session_state.corum_output = pd.DataFrame()
if 'stringdb_url' not in st.session_state:
    st.session_state.stringdb_url = None
if 'mode' not in st.session_state:
    st.session_state.mode = None  # "analysis", "demo", "reanalyze"

# -------------------------
# UI Layout
# -------------------------
st.image("logo/logo.png")
st.title("ComplexFinder")
st.markdown("**ComplexFinder searches for complexes by using the abundance values in your ProteinGroups file!**")
st.markdown("Upload or reuse results to analyze protein complexes, or checkout the demo!")

# Output file names (internal paths)
output_pdf = "WholeSheet_FinalPlots_Grouped.pdf"
output_tsv = "WholeSheet_Results.tsv"
peak_groupings_file = "WholeSheet_PeakMatchedGroupings.tsv"

# -------------------------
# Mode Selection Buttons
# -------------------------
col1, col2, col3 = st.columns(3)
with col1:
    if st.button("Run ComplexFinder Analysis"):
        st.session_state.mode = "analysis"
with col2:
    if st.button("Run Demo"):
        st.session_state.mode = "demo"
with col3:
    if st.button("Reanalyze past output"):
        st.session_state.mode = "reanalyze"

# -------------------------
# Parameters (only relevant for fresh analysis)
# -------------------------
st.sidebar.header("⚙️ Parameters")
peak_prominence_factor = st.sidebar.slider("Peak Prominence Factor", 0.1, 3.0, 1.0, 0.1, help="Modifies the prominency of a peak while it is being used to match with other peaks. Decreasing the factor can make it easier for peaks to be matched.")
user_threshold = st.sidebar.slider("OGSS Threshold", 0.0, 1.0, 0.95, 0.01, help="Threshold for OGSS (Observed Group Similarity Score). Higher values mean stricter complex matching.")
flatten_width = st.sidebar.slider("Flatten Width", 0, 5, 1, 1, help="Determines how many slices next to either side of the peak's apex will be flattened as alternate combinations of peaks in a profile are used for iterative matching. Less slices being flattened adjacent to wider peaks can lead to artifact peaks being formed.")
num_sig_peaks = st.sidebar.number_input("Significant Peaks Per Protein", min_value=1, max_value=10, value=3, help="How many of the most prominent peaks to consider for each protein during the top most level complex detection.")

# -------------------------
# MODE: Fresh Analysis
# -------------------------
if st.session_state.mode == "analysis":
    uploaded_file = st.file_uploader("Upload an Abundances.tsv file (a minimized verison of the the ProteinGroups file with the Protein IDs, and Gene names columns, plus the column with abundance values, minimally numbered; check demo for reference", type=["tsv"])
    if st.button("Start Analysis"):
        if not uploaded_file:
            st.error("Please upload an abundance TSV file before running analysis.")
        else:
            input_file = "uploaded_data.tsv"
            with open(input_file, "wb") as f:
                f.write(uploaded_file.read())

            st.session_state.demo_mode = False
            st.session_state.analysis_done = False

            start_time = time.time()
            run_complex_finder(
                input_file=input_file,
                output_pdf=output_pdf,
                output_tsv=output_tsv,
                peak_groupings_file=peak_groupings_file,
                peak_prominence_factor=peak_prominence_factor,
                user_threshold=user_threshold,
                flatten_width=flatten_width,
                num_sig_peaks=num_sig_peaks
            )
            elapsed = time.time() - start_time

            st.session_state.analysis_done = True
            st.session_state.result_files = {
                "Results TSV": output_tsv,
                "Plots PDF": output_pdf,
                "Peak Groupings": peak_groupings_file,
                "Abundances": input_file
            }
            st.success(f"Analysis completed in {elapsed:.2f} seconds.")

# -------------------------
# MODE: Demo
# -------------------------
elif st.session_state.mode == "demo":
    demo_folder = "demo"
    demo_files = {
        "Results TSV": os.path.join(demo_folder, "WholeSheet_Results.tsv"),
        "Plots PDF": os.path.join(demo_folder, "WholeSheet_FinalPlots_Grouped.pdf"),
        "Peak Groupings": os.path.join(demo_folder, "WholeSheet_PeakMatchedGroupings.tsv"),
        "Abundances": os.path.join(demo_folder, "Abundances.tsv")
    }
    missing = [name for name, path in demo_files.items() if not os.path.exists(path)]
    if missing:
        st.error(f"Demo files are missing: {', '.join(missing)}")
    else:
        shutil.copy(demo_files["Abundances"], "uploaded_data.tsv")
        shutil.copy(demo_files["Peak Groupings"], peak_groupings_file)
        st.session_state.analysis_done = True
        st.session_state.demo_mode = True
        st.session_state.result_files = demo_files
        st.success("Demo mode activated. You can now explore results.")

# -------------------------
# MODE: Reanalyze past output
# -------------------------
elif st.session_state.mode == "reanalyze":
    st.markdown("### Reanalyze previous output")
    st.info("Upload your **Abundances.tsv** and **WholeSheet_PeakMatchedGroupings.tsv** from a previous run.")
    abundances_file = st.file_uploader("Upload Abundances file", type=["tsv"], key="reanalyze_abund")
    peak_file = st.file_uploader("Upload Peak Groupings file", type=["tsv"], key="reanalyze_peak")

    if abundances_file and peak_file:
        with open("uploaded_data.tsv", "wb") as f:
            f.write(abundances_file.read())
        with open(peak_groupings_file, "wb") as f:
            f.write(peak_file.read())

        st.session_state.analysis_done = True
        st.session_state.demo_mode = False
        st.session_state.result_files = {
            "Abundances": "uploaded_data.tsv",
            "Peak Groupings": peak_groupings_file
        }
        st.success("Files loaded. You can now run CORUM, STRING-db, plotting, or complex analysis.")

# -------------------------
# Download Results
# -------------------------
if st.session_state.analysis_done:
    st.subheader("📂 Download Results")
    for label, path in st.session_state.result_files.items():
        if os.path.exists(path):
            with open(path, "rb") as f:
                st.download_button(label=f"Download {label}", data=f, file_name=os.path.basename(path))

# -------------------------
# Load available protein names
# -------------------------
protein_list = []
if os.path.exists("uploaded_data.tsv"):
    try:
        df_abundance = pd.read_csv("uploaded_data.tsv", sep="\t")
        if 'Gene names' in df_abundance.columns:
            protein_list = sorted(df_abundance['Gene names'].dropna().unique().tolist())
    except Exception as e:
        st.warning(f"Could not read abundance file for protein list: {e}")

# -------------------------
# Suggestion helper
# -------------------------
def show_suggestions(query, all_items):
    q = query.strip().lower()
    starts = [p for p in all_items if p.lower().startswith(q)]
    contains = [p for p in all_items if q in p.lower() and not p.lower().startswith(q)]
    suggestions = (starts + contains)[:20]
    if suggestions:
        st.caption("Suggestions (type the exact name to select):")
        st.markdown("\n".join([f"- {s}" for s in suggestions]))
    else:
        st.caption("No suggestions found.")

# -------------------------
# CORUM Analysis
# -------------------------
st.subheader("CORUM Complex Lookup")
st.markdown("**This section lets you check a protein and it's interactors against known CORUM complexes. This is only possible for **Human** complexes**")
if 'corum_query' not in st.session_state:
    st.session_state.corum_query = ""
if 'corum_query_submitted' not in st.session_state:
    st.session_state.corum_query_submitted = False

def _mark_corum_submitted():
    st.session_state.corum_query_submitted = True

st.text_input("Enter a protein name for CORUM analysis (Press enter or click away to get search suggestions)",
              key="corum_query", on_change=_mark_corum_submitted)
if st.session_state.corum_query.strip() and protein_list and st.session_state.corum_query_submitted:
    show_suggestions(st.session_state.corum_query, protein_list)

corum_file = "corum_humanComplexes.txt"
if st.button("Run CORUM Analysis"):
    if not os.path.exists(peak_groupings_file):
        st.error("Peak groupings file not found.")
    elif not st.session_state.corum_query.strip():
        st.error("Please enter a protein name.")
    else:
        with st.spinner("Searching CORUM database..."):
            st.session_state.corum_output = run_corum_analysis(
                query_protein=st.session_state.corum_query.strip(),
                corum_file=corum_file,
                peak_groupings_file=peak_groupings_file
            )
if not st.session_state.corum_output.empty:
    st.dataframe(st.session_state.corum_output)

# -------------------------
# STRING-db Analysis
# -------------------------
st.subheader("STRING-db Interaction Network")
st.markdown("**This section lets you check a protein and it's interactors as a STRING-db network**")
if 'stringdb_query' not in st.session_state:
    st.session_state.stringdb_query = ""
if 'stringdb_query_submitted' not in st.session_state:
    st.session_state.stringdb_query_submitted = False

def _mark_stringdb_submitted():
    st.session_state.stringdb_query_submitted = True

st.text_input("Enter a protein name for STRING-db analysis (Press enter or click away to get search suggestions)",
              key="stringdb_query", on_change=_mark_stringdb_submitted)
if st.session_state.stringdb_query.strip() and protein_list and st.session_state.stringdb_query_submitted:
    show_suggestions(st.session_state.stringdb_query, protein_list)

organism = st.text_input("Enter organism name", value="Homo sapiens")
if st.button("Run STRING-db Analysis"):
    if not os.path.exists(peak_groupings_file):
        st.error("Peak groupings file not found.")
    elif not st.session_state.stringdb_query.strip():
        st.error("Please enter a protein name.")
    else:
        with st.spinner("Generating STRING-db URL..."):
            st.session_state.stringdb_url = generate_stringdb_url(
                st.session_state.stringdb_query.strip(), organism
            )
if st.session_state.stringdb_url:
    st.markdown(f"[🔗 Open STRING-db Network]({st.session_state.stringdb_url})", unsafe_allow_html=True)

# -------------------------
# Profile Plotting
# -------------------------
st.subheader("Plot Protein & Interactor Profiles")
st.markdown("**This section lets you check plot the profiles of a protein and it's interactors**")
if 'plot_query' not in st.session_state:
    st.session_state.plot_query = ""
if 'plot_query_submitted' not in st.session_state:
    st.session_state.plot_query_submitted = False

def _mark_plot_query_submitted():
    st.session_state.plot_query_submitted = True

st.text_input("Enter a protein name to plot profiles (Press enter or click away to get search suggestions)",
              key="plot_query", on_change=_mark_plot_query_submitted)
if st.session_state.plot_query.strip() and protein_list and st.session_state.plot_query_submitted:
    show_suggestions(st.session_state.plot_query, protein_list)

peak_count_filter = st.number_input("Peak match filter (optional)", min_value=0, step=1, value=0, help="Filter interactors based on how many peaks matched between them")
if st.button("Plot Profiles"):
    if not os.path.exists(peak_groupings_file) or not os.path.exists("uploaded_data.tsv"):
        st.error("Required files not found.")
    elif not st.session_state.plot_query.strip():
        st.error("Please enter a protein name.")
    else:
        fig, error_msg = plot_protein_profiles(
            peak_groupings_file,
            st.session_state.plot_query.strip(),
            "uploaded_data.tsv",
            peak_count_filter if peak_count_filter > 0 else None
        )
        if fig:
            st.pyplot(fig)
        else:
            st.error(error_msg)

# -------------------------
# Complex Analysis
# -------------------------
st.subheader("Complex Analysis")
st.markdown("**This section lets you check a known human complex and whether the interactors in the complex were found together by ComplexFinder**")
corum_txt_file = "corum.txt"
complex_list = []
if os.path.exists(corum_txt_file):
    try:
        df_complex = pd.read_csv(corum_txt_file, sep="\t")
        if "Complex" in df_complex.columns:
            complex_list = sorted(df_complex["Complex"].dropna().unique().tolist())
    except Exception as e:
        st.warning(f"Could not read corum.txt: {e}")

if 'complex_query' not in st.session_state:
    st.session_state.complex_query = ""
if 'complex_query_submitted' not in st.session_state:
    st.session_state.complex_query_submitted = False

def _mark_complex_submitted():
    st.session_state.complex_query_submitted = True

st.text_input("Enter a complex name (Press enter or click away to get search suggestions)", key="complex_query", on_change=_mark_complex_submitted)
if st.session_state.complex_query.strip() and complex_list and st.session_state.complex_query_submitted:
    show_suggestions(st.session_state.complex_query, complex_list)

peak_filter_val = st.number_input("Peak match filter (optional)", min_value=0, step=1, value=0,
                                  help="Filter interactors in the complex based on how many peaks matched between them")
 
if st.button("Run Complex Analysis"):
    if not os.path.exists(peak_groupings_file):
        st.error("Peak groupings file not found.")
    elif not os.path.exists(corum_txt_file):
        st.error("corum.txt not found.")
    elif not st.session_state.complex_query.strip():
        st.error("Please enter a complex name.")
    else:
        with st.spinner("Running complex analysis..."):
            results = run_complex_analysis(
                st.session_state.complex_query.strip(),
                corum_txt_file,
                peak_groupings_file,
                peak_filter_val if peak_filter_val > 0 else None
            )
            if results["error"]:
                st.error(results["error"])
            else:
                st.write(f"**Confidence:** {results['confidence']:.2f}%")
                if results["graph"]:
                    st.plotly_chart(results["graph"], use_container_width=True)




# -------------------------
# Reset App
# -------------------------
if st.button("🔁 Reset App"):
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    st.rerun()
