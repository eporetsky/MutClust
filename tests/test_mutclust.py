import pytest
import pandas as pd
from mutclust.clustering import filter_and_apply_decay, run_clusterone_on_edges

def test_filter_and_apply_decay():
    # Create a mock long array
    mr_df = pd.DataFrame({
        "Gene1": ["GeneA", "GeneB"],
        "Gene2": ["GeneB", "GeneC"],
        "MR": [10, 20]  # Adjusted MR values to ensure ED > 0.01
    })
    mr_df = filter_and_apply_decay(mr_df, e_val=10)
    assert "ED" in mr_df.columns
    assert len(mr_df) == 2  # Both pairs should pass the decay filter