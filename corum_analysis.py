import pandas as pd

def run_corum_analysis(query_protein, corum_file="corum_humanComplexes.txt", peak_groupings_file="WholeSheet_PeakMatchedGroupings.tsv"):
    try:
        # Load input TSVs
        df_peaks = pd.read_csv(peak_groupings_file, sep='\t')
        df_corum = pd.read_csv(corum_file, sep='\t')

        # Check required columns
        required_peak_cols = {'Query Protein', 'Interactor', 'Peak Count'}
        required_corum_cols = {'complex_name', 'subunits_gene_name'}

        if not required_peak_cols.issubset(df_peaks.columns):
            raise ValueError(f"Missing columns in peak file: {required_peak_cols - set(df_peaks.columns)}")

        if not required_corum_cols.issubset(df_corum.columns):
            raise ValueError(f"Missing columns in CORUM file: {required_corum_cols - set(df_corum.columns)}")

        # Filter for selected query protein
        df_filtered = df_peaks[df_peaks['Query Protein'] == query_protein]

        if df_filtered.empty:
            return pd.DataFrame([{
                'Query Protein': query_protein,
                'Interactor': 'N/A',
                'Peak Count': 'N/A',
                'CORUM Complex Match': 'Query protein not found in data'
            }])

        results = []

        # Iterate over each interactor
        for _, row in df_filtered.iterrows():
            interactor = row['Interactor']
            peak_count = row['Peak Count']

            matched_complexes = []

            # Iterate over CORUM entries
            for _, corum_row in df_corum.iterrows():
                try:
                    genes = str(corum_row['subunits_gene_name']).split(';')
                    if query_protein in genes and interactor in genes:
                        matched_complexes.append(str(corum_row['complex_name']))
                except Exception:
                    continue  # skip malformed rows

            # Prepare output row
            results.append({
                'Query Protein': query_protein,
                'Interactor': interactor,
                'Peak Count': peak_count,
                'CORUM Complex Match': '; '.join(matched_complexes) if matched_complexes else 'Possibly novel'
            })

        return pd.DataFrame(results)

    except Exception as e:
        return pd.DataFrame([{
            'Query Protein': query_protein,
            'Interactor': 'Error',
            'Peak Count': 'N/A',
            'CORUM Complex Match': f'Error: {str(e)}'
        }])
