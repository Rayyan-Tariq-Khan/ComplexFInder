import streamlit as st
import os
import time
import pandas as pd
from analysis import run_complex_finder
from corum_analysis import run_corum_analysis
from stringdb_analysis import generate_stringdb_url
from plotting import plot_protein_profiles
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

# -------------------------
# UI Layout
# -------------------------
st.image("logo/logo.png")
st.title("💪🥩🧬 ComplexFinder")
st.markdown("Upload a TSV file to begin analyzing protein complexes.")

uploaded_file = st.file_uploader("Upload your TSV file", type=["tsv"])

# Parameters
st.sidebar.header("⚙️ Parameters")

peak_prominence_factor = st.sidebar.slider(
    "Peak Prominence Factor", 0.1, 3.0, 1.0, 0.1,
    help="Modifies the prominency of a peak while it is being used to match with other peaks. Decreasing the factor can make it easier for peaks to be matched."
)

user_threshold = st.sidebar.slider(
    "OGSS Threshold", 0.0, 1.0, 0.95, 0.01,
    help="Threshold for OGSS (Observed Group Similarity Score). Higher values mean stricter complex matching."
)

flatten_width = st.sidebar.slider(
    "Flatten Width", 0, 5, 1, 1,
    help="Determines how many slices next to either side of the peak's apex will be flattened as alternate combinations of peaks in a profile are used for iterative matching. Less slices being flattened adjacent to wider peaks can lead to artifact peaks being formed"
)

num_sig_peaks = st.sidebar.number_input(
    "Significant Peaks Per Protein", min_value=1, max_value=10, value=3,
    help="How many of the most prominent peaks to consider for each protein during the top most level complex detection."
)


# Output file names
output_pdf = "WholeSheet_FinalPlots_Grouped.pdf"
output_tsv = "WholeSheet_Results.tsv"
peak_groupings_file = "WholeSheet_PeakMatchedGroupings.tsv"

# -------------------------
# Run Main Analysis
# -------------------------
col1, col2 = st.columns(2)

with col1:
    if st.button("Run ComplexFinder Analysis"):
        if not uploaded_file:
            st.error("Please upload an abundance TSV file before running analysis.")
        else:
            input_file = "uploaded_data.tsv"
            with open(input_file, "wb") as f:
                f.write(uploaded_file.read())

            st.session_state.demo_mode = False
            st.session_state.analysis_done = False

            start_time = time.time()

            # Run actual analysis
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

with col2:
    if st.button("Run Demo"):
        demo_folder = "demo"
        demo_files = {
            "Results TSV": os.path.join(demo_folder, "WholeSheet_Results.tsv"),
            "Plots PDF": os.path.join(demo_folder, "WholeSheet_FinalPlots_Grouped.pdf"),
            "Peak Groupings": os.path.join(demo_folder, "WholeSheet_PeakMatchedGroupings.tsv"),
            "Abundances": os.path.join(demo_folder, "Abundances.tsv")
        }

        # Check if all demo files exist
        missing = [name for name, path in demo_files.items() if not os.path.exists(path)]
        if missing:
            st.error(f"Demo files are missing: {', '.join(missing)}")
        else:
            shutil.copy(demo_files["Abundances"], "uploaded_data.tsv")

            st.session_state.analysis_done = True
            st.session_state.demo_mode = True
            st.session_state.result_files = demo_files
            st.success("Demo mode activated. You can now explore results.")

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
# Further Analysis Section
# -------------------------
st.markdown("---")
st.header("🧪 Further Analysis")

# CORUM analysis
st.subheader("🔍 CORUM Complex Lookup")
st.markdown("Check which known CORUM complexes a protein and its interactors belong to.")
corum_protein_query = st.text_input("Enter a protein name for CORUM analysis")
corum_file = "corum_humanComplexes.txt"

if st.button("Run CORUM Analysis"):
    if not os.path.exists(peak_groupings_file):
        st.error("Peak groupings file not found. Please run the main analysis first.")
    elif not corum_protein_query.strip():
        st.error("Please enter a protein name.")
    else:
        with st.spinner("Searching CORUM database..."):
            st.session_state.corum_output = run_corum_analysis(
                query_protein=corum_protein_query.strip(),
                corum_file=corum_file,
                peak_groupings_file=peak_groupings_file
            )

if not st.session_state.corum_output.empty:
    st.success("CORUM analysis complete.")
    st.dataframe(st.session_state.corum_output)

# STRING-db analysis
st.subheader("🌐 STRING-db Interaction Network")
st.markdown("Generate a STRING-db network view for a specific protein.")
stringdb_protein_query = st.text_input("Enter a protein name for STRING-db analysis")
organism = st.text_input("Enter organism name (e.g., Homo sapiens)", value="Homo sapiens")

if st.button("Run STRING-db Analysis"):
    if not os.path.exists(peak_groupings_file):
        st.error("Peak groupings file not found. Please run the main analysis first.")
    elif not stringdb_protein_query.strip():
        st.error("Please enter a protein name.")
    else:
        with st.spinner("Generating STRING-db URL..."):
            st.session_state.stringdb_url = generate_stringdb_url(stringdb_protein_query.strip(), organism)

if st.session_state.stringdb_url:
    st.success("STRING-db URL generated.")
    st.markdown(f"[🔗 Open STRING-db Network]({st.session_state.stringdb_url})", unsafe_allow_html=True)

# -------------------------
# Profile Plotting
# -------------------------
st.subheader("📈 Plot Protein & Interactor Profiles")
plot_query = st.text_input("Enter a protein name to plot profiles")
if st.button("Plot Profiles"):
    if not os.path.exists(peak_groupings_file) or not os.path.exists("uploaded_data.tsv"):
        st.error("Required files not found. Please run the analysis or load demo.")
    elif not plot_query.strip():
        st.error("Please enter a protein name.")
    else:
        abundance_file = "uploaded_data.tsv"
        fig, error_msg = plot_protein_profiles(peak_groupings_file, plot_query.strip(), abundance_file)
        if fig:
            st.pyplot(fig)
        else:
            st.error(error_msg)

# -------------------------
# Reset App
# -------------------------
if st.button("🔁 Reset App"):
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    st.rerun()

