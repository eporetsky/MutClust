import pytest
import pandas as pd
import numpy as np
from mutclust.pca_analysis import calculate_eigen_genes

def test_empty_cluster_error():
    """Test that empty clusters raise an error."""
    expression_data = pd.DataFrame({
        'Sample1': [1.0],
        'Sample2': [2.0]
    }, index=['Gene1'])
    
    gene_clusters = [
        ['Gene1'],
        []  # Empty cluster
    ]
    
    with pytest.raises(ValueError, match="Cannot perform PCA on empty clusters"):
        calculate_eigen_genes(expression_data, gene_clusters)

def test_single_gene_clusters():
    """Test that single-gene clusters work correctly."""
    expression_data = pd.DataFrame({
        'Sample1': [1.0, 2.0],
        'Sample2': [3.0, 4.0]
    }, index=['Gene1', 'Gene2'])
    
    gene_clusters = [
        ['Gene1'],
        ['Gene2']
    ]
    
    eigen_genes = calculate_eigen_genes(expression_data, gene_clusters)
    
    # Check that each cluster's values match the original expression
    assert np.allclose(eigen_genes['Cluster_0'].values, 
                      expression_data.loc['Gene1'].values, rtol=1e-5)
    assert np.allclose(eigen_genes['Cluster_1'].values, 
                      expression_data.loc['Gene2'].values, rtol=1e-5)

def test_large_dataset_parallel():
    """Test that the function can handle larger datasets in parallel."""
    # Create a larger dataset with 100 genes and 50 samples
    np.random.seed(42)
    n_genes = 100
    n_samples = 50
    
    # Generate random expression data
    expression_data = pd.DataFrame(
        np.random.randn(n_genes, n_samples),
        index=[f'Gene{i}' for i in range(n_genes)],
        columns=[f'Sample{i}' for i in range(n_samples)]
    )
    
    # Create 10 clusters with 10 genes each
    gene_clusters = [
        [f'Gene{i}' for i in range(j*10, (j+1)*10)]
        for j in range(10)
    ]
    
    # Calculate eigen-genes
    eigen_genes = calculate_eigen_genes(expression_data, gene_clusters)
    
    # Check the output format
    assert isinstance(eigen_genes, pd.DataFrame)
    assert eigen_genes.shape == (n_samples, 10)  # 50 samples, 10 clusters
    assert all(eigen_genes.columns == [f'Cluster_{i}' for i in range(10)])
    
    # Check that each cluster's eigen-gene has the right dimensions
    for i in range(10):
        assert len(eigen_genes[f'Cluster_{i}']) == n_samples
        # Check that the values are not all zeros
        assert not np.allclose(eigen_genes[f'Cluster_{i}'], 0) 