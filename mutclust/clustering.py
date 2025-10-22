"""
MutClust: Clustering Utilities

Functions to filter mutual rank pairs, apply exponential decay weighting, and perform clustering (via ClusterONE).
"""
import pandas as pd
import numpy as np
import os
import tempfile
import shutil
import sys


def check_clusterone():
    if shutil.which("clusterone") is None:
        print(
            "WARNING: 'clusterone' executable not found in PATH. Please install it (e.g. conda install -c bioconda clusterone) "
            "and ensure it is available from the command line. Exiting without running clustering."
        )
        sys.exit(0)

def filter_and_apply_decay(mr_df, e_val):
    mr_df = mr_df[mr_df["MR"] > 0]  # Why do I have zero values?
    mr_df = mr_df.copy()
    mr_df['ED'] = np.exp(-(mr_df['MR'] - 1.0) / e_val)
    mr_df = mr_df[mr_df['ED'] >= 0.01]
    mr_df = mr_df[["Gene1", "Gene2", "ED"]]
    return mr_df

def run_clusterone_on_edges(edge_df, min_weight=0.1, pval_threshold=0.1):
    check_clusterone()
    # Prepare edge file as source, target, weight for clusterone
    with tempfile.TemporaryDirectory() as tmpdir:
        edge_path = os.path.join(tmpdir, "clusterone.tsv")
        cluster_path = os.path.join(tmpdir, "clusters.csv")
        edge_df.to_csv(edge_path, index=False, header=False, sep="\t")
        os.system(f"cp {edge_path} test.tsv")
        os.system(f'clusterone {edge_path} -f edge_list -F csv > {cluster_path}')
        clusters = pd.read_csv(cluster_path)
        clusters = clusters[clusters["P-value"] < pval_threshold]
        rows = []
        for idx, row in clusters.iterrows():
            cluster_id = f"c{idx+1}"
            pval = row["P-value"]
            for gene in str(row["Members"]).split():
                rows.append([cluster_id, gene, pval])
        result = pd.DataFrame(rows, columns=["clusterID", "geneID", "pval"])
        return result