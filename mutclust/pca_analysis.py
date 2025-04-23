import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import multiprocessing
from functools import partial

def _process_cluster(args):
    """
    Helper function to process a single cluster in parallel.
    """
    cluster_id, cluster_genes, expression_data = args
    
    # Get expression data for genes in the cluster
    cluster_expression = expression_data.loc[cluster_genes]
    
    # For single-gene clusters, the eigen-gene is just the gene's expression
    if len(cluster_genes) == 1:
        pc1 = cluster_expression.values
    else:
        # Perform PCA
        pca = PCA(n_components=1)
        pc1 = pca.fit_transform(cluster_expression.T)  # Transpose to get samples as rows
    
    # Return the results with the cluster ID from the file
    return pd.Series(pc1.flatten(), 
                    index=expression_data.columns,
                    name=cluster_id)

def calculate_eigen_genes(expression_data, gene_clusters):
    """
    Calculate the first principal component (PC1) for each gene cluster.
    
    Args:
        expression_data (pd.DataFrame): Gene expression matrix with genes as rows and samples as columns
        gene_clusters (list): List of lists containing gene IDs for each cluster
        
    Returns:
        pd.DataFrame: DataFrame with clusters as rows and samples as columns, containing PC1 values
        
    Raises:
        ValueError: If any cluster is empty
    """
    # Check for empty clusters
    if any(len(cluster) == 0 for cluster in gene_clusters):
        raise ValueError("Cannot perform PCA on empty clusters")
    
    # Prepare arguments for parallel processing
    args = [(f"c{i+1}", genes, expression_data) 
            for i, genes in enumerate(gene_clusters)]
    
    # Use multiprocessing to process clusters in parallel
    num_processes = multiprocessing.cpu_count()
    with multiprocessing.Pool(processes=num_processes) as pool:
        eigen_genes = pool.map(_process_cluster, args)
    
    # Convert to DataFrame
    eigen_genes_df = pd.DataFrame(eigen_genes).T
    return eigen_genes_df 