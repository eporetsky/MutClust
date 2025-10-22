import os
import sys
import click
import pandas as pd

# Load the MutClust modules
from mutclust.mutual_rank import calculate_mutual_rank
from mutclust.clustering import filter_and_apply_decay, run_clusterone_on_edges

@click.group()
def cli():
    """MutClust: Mutual rank-based coexpression and clustering analysis."""
    pass

@cli.command()
@click.option('--input', '-i', required=True, help='Path to RNA-seq dataset (TSV format).')
@click.option('--output', '-o', required=True, help='Path to output MR file (TSV).')
@click.option('--mr-threshold', '-m', type=float, default=100, help='Mutual rank threshold.')
@click.option('--threads', '-t', type=int, default=4, help='Number of threads for correlation calculation.')
@click.option('--log2', is_flag=True, default=False, help='Apply log2(x+1) transform before correlation.')
def mr(input, output, mr_threshold, threads, log2):
    """Calculate mutual rank from an expression dataset."""
    calculate_mutual_rank(input, threads, mr_threshold, output, log2=log2)
    click.echo(f"Saved filtered and decayed MR pairs to {output}")

@cli.command()
@click.option('--input', '-i', required=True, help='Path to Mutual Rank (MR) table (TSV format).')
@click.option('--output', '-o', required=True, help='Path to output clusters file (TSV).')
@click.option('--e_value', '-e', type=float, default=10, help='Exponential decay constant (default: 10).')
def cls(input, output, e_value):
    """Run clustering analysis on a given MR table using ClusterONE."""
    mr_df = pd.read_csv(input, sep="\t")
    # Apply exponential decay
    # Assume all geneIDs are string-mapped, or map as needed
    long_array = filter_and_apply_decay(mr_df, e_val=e_value)
    # Run ClusterONE
    cluster_df = run_clusterone_on_edges(long_array)
    # Ensure output directory exists
    output_dir = os.path.dirname(output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    cluster_df.to_csv(output, sep="\t", index=False)
    click.echo(f"Saved clusters to {output}")

if __name__ == "__main__":
    cli()