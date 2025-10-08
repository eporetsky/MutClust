import os
import sys
import click
import pandas as pd
from contextlib import redirect_stdout

# Load the MutClust modules
from mutclust.mutual_rank_analysis import calculate_correlation_matrix, calculate_mutual_rank
from mutclust.gene_clustering import filter_to_long_array, filter_and_apply_decay, create_graph_from_dataframe, leiden_clustering, convert_clusters_to_gene_ids
from mutclust.go_enrichment import run_go_enrichment
from mutclust.annotate import add_gene_annotations
from mutclust.pca_analysis import calculate_eigen_genes

@click.group()
def cli():
    """MutClust: Mutual rank-based coexpression, clustering, and GO enrichment analysis."""
    pass

@cli.command()
@click.option('--input', '-i', required=True, help='Path to RNA-seq dataset (TSV format).')
@click.option('--output', '-o', required=True, help='Path to output MR file (TSV).')
@click.option('--mr-threshold', '-m', type=float, default=100, help='Mutual rank threshold.')
@click.option('--e-value', '-e', type=float, default=10, help='Exponential decay e-value.')
@click.option('--threads', '-t', type=int, default=4, help='Number of threads for correlation calculation.')
@click.option('--save-intermediate', is_flag=True, help='Save intermediate files (PCC, MR, filtered pairs).')
def mr(input, output, mr_threshold, e_value, threads, save_intermediate):
    """Calculate mutual rank from an expression dataset."""
    expression_data = pd.read_csv(input, sep="\t", index_col=0)
    corr_np, gene_ids = calculate_correlation_matrix(input, threads)
    mr_np = calculate_mutual_rank(corr_np)
    mr_df = pd.DataFrame(mr_np, index=gene_ids, columns=gene_ids)
    click.echo("Completed calculating correlation matrix and mutual rank.")

    if save_intermediate:
        pd.DataFrame(corr_np, index=gene_ids, columns=gene_ids).to_csv(os.path.splitext(output)[0] + ".corr.tsv", sep="\t")
        mr_df.to_csv(os.path.splitext(output)[0] + ".mr.tsv", sep="\t")

    long_array = filter_to_long_array(mr_df, threshold=mr_threshold)
    gene_id_mapping = mr_df.index.to_series().reset_index(drop=True)
    long_array = filter_and_apply_decay(gene_id_mapping, long_array, e_val=e_value)
    long_array.to_csv(output, sep="\t", index=None)
    click.echo(f"Saved filtered and decayed MR pairs to {output}")

@cli.command()
@click.option('--input', '-i', required=True, help='Path to Mutual Rank (MR) table (TSV format).')
@click.option('--output', '-o', required=True, help='Path to output clusters file (TSV).')
@click.option('--annotations', '-a', help='Path to gene annotation file.')
@click.option('--resolution', '-r', type=float, default=0.1, help='Leiden clustering resolution parameter.')
@click.option('--eigengene', help='Path to RNA-seq dataset (TSV format) for eigen-gene calculation. If provided, eigen-genes will be calculated and saved alongside clusters.')
def cls(input, output, annotations, resolution, eigengene):
    """Run clustering analysis on a given MR table."""
    long_array = pd.read_csv(input, sep="\t")
    graph = create_graph_from_dataframe(long_array)
    click.echo("Completed creating graph from MR table.")

    clusters = leiden_clustering(graph, resolution_parameter=resolution)
    gene_id_clusters = convert_clusters_to_gene_ids(graph, clusters)
    click.echo("Completed Leiden clustering.")

    cluster_df = []
    for cluster_id, cluster_genes in enumerate(gene_id_clusters, 1):
        for gene in cluster_genes:
            cluster_df.append([f"c{cluster_id}", gene])
    cluster_df = pd.DataFrame(cluster_df, columns=["clusterID", "geneID"])
    if annotations:
        cluster_df = add_gene_annotations(cluster_df, annotations)
    cluster_df.to_csv(output, sep="\t", index=False)
    click.echo(f"Saved clusters to {output}")

    if eigengene:
        expression_data = pd.read_csv(eigengene, sep="\t", index_col=0)
        clusters_df = pd.read_csv(output, sep="\t")
        gene_clusters = clusters_df.groupby('clusterID')['geneID'].apply(list).tolist()
        eigen_genes = calculate_eigen_genes(expression_data, gene_clusters)
        eigen_genes.index.name = "geneID"
        eigen_out = os.path.splitext(output)[0] + ".eigen.tsv"
        eigen_genes.to_csv(eigen_out, sep="\t")
        click.echo(f"Saved eigen-genes to {eigen_out}")

@cli.command()
@click.option('--clusters', '-c', required=True, help='Path to clusters file (TSV format).')
@click.option('--go-obo', '-go', required=True, help='Path to Gene Ontology (GO) OBO file.')
@click.option('--go-gaf', '-gf', required=True, help='Path to GO annotation file (GAF format).')
@click.option('--output', '-o', required=True, help='Output prefix for results.')
@click.option('--expression', help='Path to RNA-seq dataset (TSV format) for background gene set.')
def enr(clusters, go_obo, go_gaf, output, expression):
    """Run GO enrichment analysis on clusters."""
    clusters_df = pd.read_csv(clusters, sep="\t")
    gene_clusters = clusters_df.groupby('clusterID')['geneID'].apply(list).tolist()
    if expression:
        expression_data = pd.read_csv(expression, sep="\t", index_col=0)
        background_genes = set(expression_data.index)
    else:
        background_genes = set(clusters_df['geneID'])
    with open(os.devnull, 'w') as fnull:
        with redirect_stdout(fnull):
            run_go_enrichment(gene_clusters, go_obo, go_gaf, background_genes, output)
    click.echo(f"Saved GO enrichment results to {output}_go_enrichment_results.tsv")

if __name__ == "__main__":
    cli()