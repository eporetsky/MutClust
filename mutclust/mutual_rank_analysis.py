import pandas as pd
import numpy as np
from pynetcor.cor import corrcoef
from scipy import sparse
import gc
import pigz

def calculate_mutual_rank(file_path, threads, mr_threshold, output_path):
    df = pd.read_csv(file_path, sep="\t", index_col="geneID")
    # Zero variance row filter
    df = df.loc[df.var(axis=1) > 0]
    # log2(cpm+1) transformation
    df = np.log2(df + 1)

    print("calculating correlation matrix")
    corr_np = corrcoef(df, threads=threads)

    # Save gene IDs before deleting df
    gene_ids = df.index.to_numpy()
    del df
    gc.collect()

    print("Converting to float32")
    corr_np = corr_np.astype(np.float32)

    # Mask diagonal before ranking (self-correlation should not affect ranks)
    #corr_np_masked = corr_np.copy()
    #np.fill_diagonal(corr_np_masked, -np.inf)

    #row_ranks = np.argsort(-corr_np_masked, axis=1).argsort(axis=1) + 1
    #col_ranks = np.argsort(-corr_np_masked, axis=0).argsort(axis=0) + 1
    
    import concurrent.futures
    def rank_row(row):
        # Clip and convert to uint16 immediately to save memory
        ranks = np.argsort(-row).argsort() + 1
        return np.clip(ranks, 1, 250).astype(np.uint16)
    
    print("Ranking rows")
    # Cap ranks at 250 to keep product <= 62500 (fits in uint16)
    with concurrent.futures.ThreadPoolExecutor() as executor:
        row_ranks = np.array(list(executor.map(rank_row, corr_np)), dtype=np.uint16)

    print("Ranking columns")
    with concurrent.futures.ThreadPoolExecutor() as executor:
        col_ranks = np.array(list(executor.map(rank_row, corr_np.T)), dtype=np.uint16).T
 
    del corr_np
    gc.collect()

    print("calculating MR")
    mr_np = np.sqrt(row_ranks * col_ranks).astype(np.float32)
    del row_ranks, col_ranks
    gc.collect()

    print("filtering MR and converting to sparse")
    # Set values above threshold to 0 for sparsity
    mr_np[mr_np > mr_threshold] = 0
    # Convert to sparse matrix (only stores nonzero values)
    mr_np = sparse.csr_matrix(mr_np)
    
    # Extract upper triangle only (k=1 excludes diagonal)
    mr_np = sparse.triu(mr_np, k=1)
    
    # Get nonzero indices and values (only the filtered pairs!)
    i, j = mr_np.nonzero()
    mr_np = mr_np.data
    
    print(f"Writing {len(mr_np)} gene pairs to file")
    with pigz.open(output_path, "wt") as f:
        f.write("Gene1\tGene2\tMR\n")
        for a, b, v in zip(gene_ids[i], gene_ids[j], mr_np):
            f.write(f"{a}\t{b}\t{v:.6g}\n")
    del i, j, mr_np
    gc.collect()

    return gene_ids.tolist()
