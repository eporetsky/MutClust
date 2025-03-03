import pandas as pd
import numpy as np
from pynetcor.cor import corrcoef


def calculate_correlation_matrix(file_path, threads):
    """
    Calculates the correlation matrix for a given RNA-seq dataset.

    Args:
        file_path (str): The path to the RNA-seq dataset file in TSV format.
        threads (int): The number of threads to use for the correlation calculation.

    Returns:
        tuple: A tuple containing:
            - corr_np (numpy.ndarray): The correlation matrix as a NumPy array.
            - list: A list of gene IDs from the dataset.
    """

    df = pd.read_csv(file_path, sep="\t", index_col="geneID")
    corr_np = corrcoef(df, threads=threads)
    return corr_np, df.index.tolist()


def calculate_mutual_rank(corr_np):
    """
    Calculates the mutual rank matrix for a given correlation matrix.

    Args:
        corr_np (numpy.ndarray): The correlation matrix as a NumPy array.

    Returns:
        numpy.ndarray: The mutual rank matrix as a NumPy array.
    """
    row_ranks = np.argsort(-corr_np, axis=1).argsort(axis=1) + 1
    col_ranks = np.argsort(-corr_np, axis=0).argsort(axis=0) + 1
    mr_np = np.sqrt(row_ranks * col_ranks)
    return mr_np