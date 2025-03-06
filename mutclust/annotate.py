import pandas as pd

def add_gene_annotations(cluster_df, annotations):
    if isinstance(annotations, str):  # If a file path is provided
        annotations = pd.read_csv(annotations, sep="\t")
    cluster_df = cluster_df.merge(annotations, on="geneID", how="left")
    return cluster_df