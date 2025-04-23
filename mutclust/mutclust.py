import pandas as pd

# Save eigen-genes to file
eigen_genes_df = pd.DataFrame(eigen_genes, index=gene_names)
eigen_genes_df.index.name = 'geneID'
eigen_genes_df.columns = [f'cluster_{i}' for i in range(eigen_genes_df.shape[1])]
eigen_genes_df.to_csv(output_eigen_genes, sep='\t') 