# analysis.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import pearsonr, ConstantInputWarning
from matplotlib.backends.backend_pdf import PdfPages
import warnings
from itertools import combinations
import time
import csv
from tqdm import tqdm

warnings.filterwarnings("ignore", category=ConstantInputWarning)

def auto_peak_threshold(profile, peak_prominence_factor):
    values = profile[profile > 0]
    if len(values) == 0:
        return 0
    median = np.median(values)
    std = np.std(values)
    return median + peak_prominence_factor * std

def compute_similarity(prof1, prof2, peak_thresh1, peak_thresh2):
    peaks1, _ = find_peaks(prof1, prominence=peak_thresh1)
    peaks2, _ = find_peaks(prof2, prominence=peak_thresh2)
    set1 = set(peaks1)
    set2 = set(peaks2)
    jaccard = len(set1 & set2) / len(set1 | set2) if len(set1 | set2) > 0 else 0
    dot = np.dot(prof1, prof2)
    norm1 = np.linalg.norm(prof1)
    norm2 = np.linalg.norm(prof2)
    cosine = dot / (norm1 * norm2) if norm1 > 0 and norm2 > 0 else 0
    try:
        pearson_corr, _ = pearsonr(prof1, prof2)
    except:
        pearson_corr = 0
    ogss = (jaccard + cosine + pearson_corr) / 3.0
    return ogss, jaccard, cosine, pearson_corr

def flatten_peak(profile, peak_idx, flatten_width):
    start = max(0, peak_idx - flatten_width)
    end = min(len(profile), peak_idx + flatten_width + 1)
    profile[start:end] = 0
    return profile

def plot_profiles(ax, title, profiles, slice_cols):
    slices = np.arange(1, len(slice_cols) + 1)
    for label, prof in profiles.items():
        ax.plot(slices, prof, label=label, linewidth=2)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("BN-PAGE Slice")
    ax.set_ylabel("Abundance")
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1, 1))

def run_complex_finder(
    input_file,
    output_pdf,
    output_tsv,
    peak_groupings_file,
    peak_prominence_factor,
    user_threshold,
    flatten_width,
    num_sig_peaks
):
    df = pd.read_csv(input_file, sep="\t")
    id_col, gene_col = df.columns[:2]
    slice_cols = df.columns[2:]

    start_global = time.time()
    results_rows = []
    peak_groupings_data = []
    pdf_pages = PdfPages(output_pdf)
    total_proteins = len(df)

    for idx, query_row in tqdm(df.iterrows(), total=total_proteins, desc="Processing proteins"):
        query_label = query_row[gene_col] if pd.notna(query_row[gene_col]) else query_row[id_col]
        query_profile_orig = query_row[slice_cols].astype(float).to_numpy()

        threshold = auto_peak_threshold(query_profile_orig, peak_prominence_factor)
        all_peaks, _ = find_peaks(query_profile_orig, prominence=threshold)

        if len(all_peaks) == 0:
            continue

        peak_heights = query_profile_orig[all_peaks]
        sorted_indices = np.argsort(peak_heights)[::-1]
        sorted_peaks = all_peaks[sorted_indices]

        selected_peaks = sorted_peaks[:num_sig_peaks]
        if len(sorted_peaks) > num_sig_peaks:
            cutoff_height = peak_heights[sorted_indices][num_sig_peaks - 1]
            tied_indices = [i for i in sorted_indices if peak_heights[i] == cutoff_height]
            tied_peaks = all_peaks[tied_indices]
            selected_peaks = np.unique(np.concatenate((sorted_peaks[:num_sig_peaks], tied_peaks)))

        peak_to_prots = {i: set() for i in range(1, len(selected_peaks) + 1)}
        protein_scores = {}

        for _, row in df.iterrows():
            if row.equals(query_row):
                continue
            prot_label = row[gene_col] if pd.notna(row[gene_col]) else row[id_col]
            test_profile = row[slice_cols].astype(float).to_numpy()
            t1 = auto_peak_threshold(query_profile_orig, peak_prominence_factor)
            t2 = auto_peak_threshold(test_profile, peak_prominence_factor)
            ogss, jacc, cos, pearson = compute_similarity(query_profile_orig, test_profile, t1, t2)
            if ogss >= user_threshold:
                key = (prot_label, len(selected_peaks))
                if key not in protein_scores:
                    protein_scores[key] = {"ogss": [], "jaccard": [], "cosine": [], "pearson": []}
                protein_scores[key]["ogss"].append(ogss)
                protein_scores[key]["jaccard"].append(jacc)
                protein_scores[key]["cosine"].append(cos)
                protein_scores[key]["pearson"].append(pearson)
                peak_to_prots[len(selected_peaks)].add(prot_label)

        for r in range(len(selected_peaks) - 1, 0, -1):
            for comb in combinations(selected_peaks, r):
                current_profile = query_profile_orig.copy()
                for p in comb:
                    current_profile = flatten_peak(current_profile, p, flatten_width)

                for _, row in df.iterrows():
                    if row.equals(query_row):
                        continue
                    prot_label = row[gene_col] if pd.notna(row[gene_col]) else row[id_col]
                    test_profile = row[slice_cols].astype(float).to_numpy()
                    t1 = auto_peak_threshold(current_profile, peak_prominence_factor)
                    t2 = auto_peak_threshold(test_profile, peak_prominence_factor)
                    ogss, jacc, cos, pearson = compute_similarity(current_profile, test_profile, t1, t2)
                    if ogss >= user_threshold:
                        key = (prot_label, r)
                        if key not in protein_scores:
                            protein_scores[key] = {"ogss": [], "jaccard": [], "cosine": [], "pearson": []}
                        protein_scores[key]["ogss"].append(ogss)
                        protein_scores[key]["jaccard"].append(jacc)
                        protein_scores[key]["cosine"].append(cos)
                        protein_scores[key]["pearson"].append(pearson)
                        peak_to_prots[r].add(prot_label)

        for peak_count in sorted(peak_to_prots.keys(), reverse=True):
            prots = peak_to_prots[peak_count]
            if len(prots) == 0:
                continue
            combined = {query_label: query_profile_orig}
            for p in prots:
                row = df[(df[gene_col] == p) | (df[id_col] == p)].iloc[0]
                profile = row[slice_cols].astype(float).to_numpy()
                combined[p] = profile

            fig, ax = plt.subplots(figsize=(10, 4))
            title = f"{idx} - {peak_count} {query_label} (Found {len(prots)} interactors)"
            plot_profiles(ax, title, combined, slice_cols)
            pdf_pages.savefig(fig)
            plt.close(fig)

        if protein_scores:
            for (prot_label, peak_count), scores in protein_scores.items():
                peak_groupings_data.append({
                    "Query Protein": query_label,
                    "Peak Count": peak_count,
                    "Interactor": prot_label,
                    "OGSS": np.mean(scores["ogss"]),
                    "Pearson": np.mean(scores["pearson"]),
                    "Cosine": np.mean(scores["cosine"]),
                    "Jaccard": np.mean(scores["jaccard"])
                })

            grouped_by_prot = {}
            for (prot_label, _), _ in protein_scores.items():
                if prot_label not in grouped_by_prot:
                    grouped_by_prot[prot_label] = [k for k in protein_scores if k[0] == prot_label]

            row = {
                "ID": idx,
                "Protein": query_label,
                "Interactors": ",".join(grouped_by_prot.keys()),
                "Interactors Found In Iterations": ",".join(str(len(grouped_by_prot[k])) for k in grouped_by_prot),
                "Avg OGSS Per Interactor": ",".join(f"{np.mean(protein_scores[(k, max([p for p in range(1, num_sig_peaks + 1) if (k, p) in protein_scores]) )]['ogss']):.4f}" for k in grouped_by_prot),
                "Avg Pearson Per Interactor": ",".join(f"{np.mean(protein_scores[(k, max([p for p in range(1, num_sig_peaks + 1) if (k, p) in protein_scores]) )]['pearson']):.4f}" for k in grouped_by_prot),
                "Avg Cosine Per Interactor": ",".join(f"{np.mean(protein_scores[(k, max([p for p in range(1, num_sig_peaks + 1) if (k, p) in protein_scores]) )]['cosine']):.4f}" for k in grouped_by_prot),
                "Avg Jaccard Per Interactor": ",".join(f"{np.mean(protein_scores[(k, max([p for p in range(1, num_sig_peaks + 1) if (k, p) in protein_scores]) )]['jaccard']):.4f}" for k in grouped_by_prot),
                "Final Score": f"{np.mean([np.mean(protein_scores[k]['ogss']) for k in protein_scores])*100:.2f}"
            }
            results_rows.append(row)

    pdf_pages.close()

    if results_rows:
        pd.DataFrame(results_rows).to_csv(output_tsv, sep="\t", index=False)

    with open(peak_groupings_file, "w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=[
            "Query Protein", "Peak Count", "Interactor", "OGSS", "Pearson", "Cosine", "Jaccard"
        ])
        writer.writeheader()
        for row in peak_groupings_data:
            writer.writerow(row)
            

