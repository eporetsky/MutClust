import pandas as pd

def add_gene_annotations(cluster_df, path_annot):
    annot = pd.read_csv(path_annot, sep="\t")

    cluster_df = cluster_df.merge(annot, on="geneID", how="left")

    return(cluster_df)