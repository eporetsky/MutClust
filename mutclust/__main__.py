import os
import sys
import argparse
import pandas as pd
from contextlib import redirect_stdout

# Load the MutClust modules
from mutclust.mutual_rank_analysis import calculate_correlation_matrix, calculate_mutual_rank
from mutclust.gene_clustering import filter_to_long_array, filter_and_apply_decay, create_graph_from_dataframe, leiden_clustering, convert_clusters_to_gene_ids
from mutclust.go_enrichment import run_go_enrichment
from mutclust.annotate import add_gene_annotations
from mutclust.pca_analysis import calculate_eigen_genes

def main():
    """
    Main entry point for the MutClust package.

    This function parses command-line arguments for MutClust and performs the following steps:

    1. Calculate correlation matrix and mutual rank from RNA-seq dataset.
    2. Filter, calculate exponential decay, and create an graph from mutual rank matrix.
    3. Perform Leiden clustering to group coexpressed genes.
    4. Calculate eigen-genes for each cluster.
    5. Perform GO enrichment analysis on the resulting clusters.
    6. Save results to specified output files.
    """
    
    parser = argparse.ArgumentParser(description="MutClust: Mutual rank-based coexpression analysis, Leiden clustering and GO enrichment analysis.")
    
    # Create a mutually exclusive group for input files
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--expression", "-ex", help="Path to RNA-seq dataset (TSV format).")
    group.add_argument("--mutual_rank", "-mr", help="Path to Mutual Rank file (TSV format).")
    
    # Add the rest of the input arguments
    parser.add_argument("--annotations", "-a", help="Path to gene annotation file.")
    parser.add_argument("--go_obo", "-go", help="Path to Gene Ontology (GO) OBO file.")
    parser.add_argument("--go_gaf", "-gf", help="Path to GO annotation file (GAF format).")
    parser.add_argument("--output", "-o", required=True, help="Output prefix for results.")
    parser.add_argument("--mr_threshold", "-m", type=float, default=100, help="Mutual rank threshold.")
    parser.add_argument("--e_value", "-e", type=float, default=10, help="Exponential decay e-value.")
    parser.add_argument("--resolution", "-r", type=float, default=0.1, help="Leiden clustering resolution parameter.")
    parser.add_argument("--threads", "-t", type=int, default=4, help="Number of threads for correlation calculation.")
    parser.add_argument("--save_intermediate", action="store_true", help="Save intermediate files (PCC, MR, filtered pairs).")
    parser.add_argument("--eigengene", action="store_true", default=True, help="Calculate eigen-genes for clusters.")
    args = parser.parse_args()

    # Step 1: Calculate correlation matrix and mutual rank
    if args.expression:
        expression_data = pd.read_csv(args.expression, sep="\t", index_col=0)
        corr_np, gene_ids = calculate_correlation_matrix(args.expression, args.threads)
        mr_np = calculate_mutual_rank(corr_np)
        mr_df = pd.DataFrame(mr_np, index=gene_ids, columns=gene_ids)
        print("Completed calculating correlation matrix and mutual rank.")

        # Step 2: Filter, decay, and create graph
        long_array = filter_to_long_array(mr_df, threshold=args.mr_threshold)
        gene_id_mapping = mr_df.index.to_series().reset_index(drop=True)
        long_array = filter_and_apply_decay(gene_id_mapping, long_array, e_val=args.e_value)
        long_array.to_csv(f"{args.output}.mrs.tsv", sep="\t", index=None)

    # If a Mutual Rank file is provided, load it
    if args.mutual_rank:
        long_array = pd.read_csv(args.mutual_rank, sep="\t")
        print("Completed loading MR dataframe.")

    # Step 2: Create graph from dataframe
    graph = create_graph_from_dataframe(long_array)
    print("Completed filtering, exonential decay, and creating graph.")

    # Step 3: Leiden clustering
    clusters = leiden_clustering(graph, resolution_parameter=args.resolution)
    gene_id_clusters = convert_clusters_to_gene_ids(graph, clusters)
    print("Completed Leiden clustering.")

    # Step 4: Write clusters to file first
    cluster_df = []
    for cluster_id, cluster_genes in enumerate(gene_id_clusters, 1):
        for gene in cluster_genes:
            cluster_df.append([f"c{cluster_id}", gene])
    cluster_df = pd.DataFrame(cluster_df, columns=["clusterID", "geneID"])
    if args.annotations:
        cluster_df = add_gene_annotations(cluster_df,
                                       args.annotations)
    cluster_df.to_csv(f"{args.output}.clusters.tsv", sep="\t", index=False)
    print("Completed writing clusters to file.")

    # Step 5: Calculate eigen-genes using the clusters file
    if args.eigengene:
        # Read the clusters file
        clusters_df = pd.read_csv(f"{args.output}.clusters.tsv", sep="\t")
        # Group genes by cluster
        gene_clusters = clusters_df.groupby('clusterID')['geneID'].apply(list).tolist()
        eigen_genes = calculate_eigen_genes(expression_data, gene_clusters)
        eigen_genes.index.name = "geneID"
        eigen_genes.to_csv(f"{args.output}.eigen.tsv", sep="\t")
        print("Completed eigen-gene calculation.")

    # Step 6: GO enrichment
    if args.go_obo and args.go_gaf:
        background_genes = set(gene_ids)
        # Suppress the output
        with open(os.devnull, 'w') as fnull:
            with redirect_stdout(fnull):
                run_go_enrichment(gene_id_clusters, 
                                args.go_obo,
                                args.go_gaf, 
                                background_genes,
                                args.output)
        print("Completed GO enrichment analysis.")

    print("Completed MutClust analysis and saving results.")

if __name__ == "__main__":
    main()