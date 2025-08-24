# stringdb_analysis.py

import pandas as pd
import urllib.parse

def generate_stringdb_url(query_protein, organism, peak_groupings_file="WholeSheet_PeakMatchedGroupings.tsv"):
    try:
        df = pd.read_csv(peak_groupings_file, sep='\t')

        # Filter rows for the query protein
        df_filtered = df[df['Query Protein'] == query_protein]

        # Get all unique interactors for the query protein
        interactors = set(df_filtered['Interactor'].dropna().unique())
        interactors.add(query_protein)  # include the query protein itself

        if not interactors:
            return "No interactors found for this protein."

        # Build STRING URL
        base_url = "https://string-db.org/cgi/network.pl"
        params = {
            "identifiers": "%0d".join(interactors),  # newline separator for proteins
            "species": organism
        }
        query_string = urllib.parse.urlencode(params)
        return f"{base_url}?{query_string}"

    except FileNotFoundError:
        return f"Error: Peak groupings file '{peak_groupings_file}' not found."
    except KeyError as e:
        return f"Error: Missing expected column in data: {e}"
    except Exception as e:
        return f"Error: {str(e)}"
