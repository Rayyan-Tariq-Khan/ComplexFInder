# Complexes.py
import os
import re
import itertools
import functools
import pandas as pd
import networkx as nx
import plotly.graph_objects as go


def _split_aliases(cell: str):
    """Split a peak file cell that may contain aliases like 'TMEM120A;TMPIT'."""
    if pd.isna(cell):
        return []
    return [tok.strip() for tok in str(cell).split(";") if tok.strip()]


def _split_interactors_field(s: str):
    """
    Split the Interactors field from corum.txt.
    Accepts ';' or ':' (or commas) as separators just in case.
    """
    return [tok.strip() for tok in re.split(r"[;,:]", str(s)) if tok.strip()]


@functools.lru_cache(maxsize=8)
def _present_pairs_cache(peaks_path: str, peak_filter: int | None, mtime: float):
    """
    Build a set of unordered alias pairs present in the peak file, optionally
    filtered by Peak Count == peak_filter. Each edge is stored as
    frozenset({ALIAS_A_UPPER, ALIAS_B_UPPER}).
    Cached by (path, filter, mtime) for speed.
    """
    usecols = ["Query Protein", "Interactor", "Peak Count"]
    df = pd.read_csv(peaks_path, sep="\t", usecols=usecols, dtype=str, low_memory=False)

    if peak_filter is not None:
        # Keep exact matches only
        # coerce to Int64 to safely compare while handling missing/str
        df["_pc"] = pd.to_numeric(df["Peak Count"], errors="coerce").astype("Int64")
        df = df[df["_pc"] == peak_filter].drop(columns=["_pc"])

    present = set()

    # itertuples is fast; expand aliases cross-product per row
    for row in df.itertuples(index=False, name=None):
        q_val, i_val, _pc = row
        q_aliases = _split_aliases(q_val)
        i_aliases = _split_aliases(i_val)
        if not q_aliases or not i_aliases:
            continue
        for qa in q_aliases:
            qa_u = qa.upper()
            for ia in i_aliases:
                ia_u = ia.upper()
                if qa_u and ia_u and qa_u != ia_u:
                    present.add(frozenset({qa_u, ia_u}))
    return present


def run_complex_analysis(complex_name, corum_file, peak_groupings_file, peak_filter_val=None):
    # Load corum
    corum = pd.read_csv(corum_file, sep="\t")
    if "Complex" not in corum or "Interactors" not in corum:
        return {"error": "Invalid corum.txt format.", "graph": None, "confidence": 0.0}

    row = corum[corum["Complex"].str.strip().str.upper() == complex_name.strip().upper()]
    if row.empty:
        return {"error": f"Complex '{complex_name}' not found in corum file.", "graph": None, "confidence": 0.0}

    interactors_raw = row.iloc[0]["Interactors"]
    interactors = [p.strip().upper() for p in re.split(r"[;:]", interactors_raw) if p.strip()]
    if len(interactors) < 2:
        return {"error": f"Complex '{complex_name}' has fewer than 2 interactors.", "graph": None, "confidence": 0.0}

    # All unique unordered pairs
    expected_pairs = {frozenset(x) for x in itertools.combinations(interactors, 2)}

    # Load peaks
    df = pd.read_csv(peak_groupings_file, sep="\t")
    if not {"Query Protein", "Interactor", "Peak Count"}.issubset(df.columns):
        return {"error": "Invalid peak groupings file format.", "graph": None, "confidence": 0.0}

    # Apply filter if provided
    if peak_filter_val is not None:
        df = df[df["Peak Count"] == peak_filter_val]

    present_pairs = set()
    for _, row in df.iterrows():
        q = row["Query Protein"].strip().upper()
        i = row["Interactor"].strip().upper()
        if q != i:
            present_pairs.add(frozenset([q, i]))

    # Find matches
    found_pairs = expected_pairs & present_pairs
    confidence = (len(found_pairs) / len(expected_pairs)) * 100 if expected_pairs else 0

    # Build network graph
    G = nx.Graph()
    G.add_nodes_from(interactors)
    for u, v in expected_pairs:
        G.add_edge(u, v, found=(frozenset([u, v]) in found_pairs))

    pos = nx.spring_layout(G, seed=42)
    edge_traces = []
    for u, v, d in G.edges(data=True):
        color = "green" if d["found"] else "lightgray"
        edge_traces.append(go.Scatter(
            x=[pos[u][0], pos[v][0]], y=[pos[u][1], pos[v][1]],
            mode="lines", line=dict(width=2, color=color), hoverinfo="none"
        ))

    node_trace = go.Scatter(
        x=[pos[n][0] for n in G.nodes()],
        y=[pos[n][1] for n in G.nodes()],
        mode="markers+text",
        text=list(G.nodes()),
        textposition="top center",
        marker=dict(size=15, color="skyblue", line=dict(width=2, color="black"))
    )

    fig = go.Figure(data=edge_traces + [node_trace])
    fig.update_layout(
        showlegend=False, hovermode="closest",
        margin=dict(b=20, l=5, r=5, t=40),
        title=f"Complex: {complex_name} (Green line indicates valid connection)",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
    )

    return {"graph": fig, "confidence": confidence, "error": None}

