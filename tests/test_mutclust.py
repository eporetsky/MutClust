import pytest
import pandas as pd
from mutclust.gene_clustering import filter_to_long_array, filter_and_apply_decay
from mutclust.annotate import add_gene_annotations

def test_filter_to_long_array():
    # Create a mock mutual rank DataFrame
    data = {
        "GeneA": [0, 50, 200],
        "GeneB": [50, 0, 75],
        "GeneC": [200, 75, 0],
    }
    mr_df = pd.DataFrame(data, index=["GeneA", "GeneB", "GeneC"])
    long_array = filter_to_long_array(mr_df, threshold=100)
    assert len(long_array) == 2  # Only two pairs should pass the threshold

def test_filter_and_apply_decay():
    # Create a mock long array
    long_array = pd.DataFrame({
        "Gene1": [0, 1],
        "Gene2": [1, 2],
        "MR": [10, 20]  # Adjusted MR values to ensure ED > 0.01
    })
    gene_id_mapping = {0: "GeneA", 1: "GeneB", 2: "GeneC"}
    filtered_array = filter_and_apply_decay(gene_id_mapping, long_array, e_val=10)
    assert "ED" in filtered_array.columns
    assert len(filtered_array) == 2  # Both pairs should pass the decay filter

def test_add_gene_annotations():
    # Create mock cluster and annotation data
    cluster_df = pd.DataFrame({"geneID": ["GeneA", "GeneB"]})
    annotations = pd.DataFrame({
        "geneID": ["GeneA", "GeneB"],
        "description": ["Protein A", "Protein B"]
    })
    annotated_df = add_gene_annotations(cluster_df, annotations)
    assert "description" in annotated_df.columns
    assert annotated_df.loc[0, "description"] == "Protein A"