import pandas as pd
import numpy as np
from igraph import Graph


def filter_to_long_array(mr_df, threshold=100):
    mr_df = mr_df.astype(float)
    mr_df.values[np.tril_indices_from(mr_df, k=0)] = np.inf
    rows, cols = np.where(mr_df < threshold)
    values = mr_df.values[rows, cols]
    long_array = pd.DataFrame({'Gene1': rows, 'Gene2': cols, 'MR': values})
    return long_array


def filter_and_apply_decay(gene_id_mapping, long_array, e_val):
    long_array['ED'] = np.exp(-(long_array['MR'] - 1.0) / e_val)
    long_array = long_array[long_array['ED'] >= 0.01].copy()
    long_array['Gene1'] = long_array['Gene1'].map(gene_id_mapping)
    long_array['Gene2'] = long_array['Gene2'].map(gene_id_mapping)
    return long_array


def create_graph_from_dataframe(long_array):
    graph = Graph.DataFrame(long_array.drop("MR", axis=1), directed=False, use_vids=False)
    return graph


def leiden_clustering(graph, resolution_parameter):
    partition = graph.community_leiden(weights='ED', resolution_parameter=resolution_parameter)
    return partition


def convert_clusters_to_gene_ids(graph, clusters):
    gene_id_clusters = []
    for cluster in clusters:
        if len(cluster) == 1:
            continue
        gene_ids = [graph.vs[vertex_id]["name"] for vertex_id in cluster]
        gene_id_clusters.append(gene_ids)
    return gene_id_clusters