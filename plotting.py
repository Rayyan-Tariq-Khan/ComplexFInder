import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_protein_profiles(peak_groupings_file, query_protein, abundance_file="uploaded_data.tsv", peak_count_filter=None):
    try:
        df_abundance = pd.read_csv(abundance_file, sep="\t")
        df_peaks = pd.read_csv(peak_groupings_file, sep="\t")

        # Filter by query protein
        df_query = df_peaks[df_peaks['Query Protein'] == query_protein]

        # Apply peak count filter if provided
        if peak_count_filter is not None:
            df_query = df_query[df_query['Peak Count'] == peak_count_filter]

        # Get interactors
        interactors = df_query['Interactor'].unique().tolist()
        proteins_to_plot = [query_protein] + interactors

        # Filter abundance data
        df_filtered = df_abundance[df_abundance['Gene names'].isin(proteins_to_plot)].copy()
        if df_filtered.empty:
            return None, f"No abundance data found for {query_protein} or its interactors."

        # Melt for plotting
        slice_cols = df_abundance.columns[1:]  # Assuming first col is 'Gene names'
        df_melted = df_filtered.melt(
            id_vars='Gene names',
            value_vars=slice_cols,
            var_name='Slice',
            value_name='Abundance'
        )

        # Ensure numeric sort of slices
        df_melted['Slice'] = pd.to_numeric(df_melted['Slice'], errors='coerce')
        df_melted = df_melted.dropna(subset=['Slice'])
        df_melted.sort_values(by='Slice', inplace=True)

        # Plotting
        plt.figure(figsize=(12, 6))
        sns.lineplot(data=df_melted, x='Slice', y='Abundance', hue='Gene names', marker='o', ci=None)
        plt.title(f"Abundance Profile: {query_protein} and Interactors")
        plt.xlabel("BN-PAGE Slice")
        plt.ylabel("Abundance")
        plt.legend(title="Protein", bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        plt.tight_layout()

        return plt, None

    except Exception as e:
        return None, f"Error during plotting: {str(e)}"
