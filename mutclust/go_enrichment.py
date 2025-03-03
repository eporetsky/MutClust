from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_gaf
import multiprocessing
import pandas as pd

# Function to calculate enrichment for a single cluster
def calculate_enrichment(item):
    global go_study

    name, genes = item
    
    try:
        results = go_study.run_study(genes)
        return_list = []
        for result in results:
            # Filter significant results
            if result.p_bonferroni < 0.05 and result.enrichment == "e":
                # Calculate fold change
                fc = (result.ratio_in_study[0] / result.ratio_in_study[1]) / (result.ratio_in_pop[0] / result.ratio_in_pop[1])
                if fc < 2:
                    continue
                return_list.append([name, result.goterm.namespace, len(genes), result.GO, result.p_bonferroni, fc, result.goterm.name])
        return return_list
    except Exception as e:
        print(f"Error processing cluster {name}: {e}")
        return []

def run_go_enrichment(gene_id_clusters, 
                      obo_file, 
                      gaf_file, 
                      background_genes,
                      output):
    # go_study cannot be pickled for hyperthreading, using global instead
    global go_study

    #print(gene_id_clusters)

    # Load the OBO file
    go_dag = GODag(obo_file)

    # Load the GAF file
    go_annotations = read_gaf(gaf_file)

    # Prepare GOATOOLS enrichment study
    go_study = GOEnrichmentStudy(
        background_genes,  # Background gene set
        go_annotations,    # Gene-to-GO mapping
        go_dag,            # GO DAG (ontology structure)
        methods=["bonferroni"]  # Multiple testing correction method
    )

    # Parallel enrichment analysis
    num_processes = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=num_processes)
    results = []
    #for result in pool.imap(lambda item: calculate_enrichment(item, go_study), enumerate(gene_id_clusters)):
    for result in pool.imap(calculate_enrichment, enumerate(gene_id_clusters)):
        results += result
    pool.close()
    pool.join()

    # Step 5: Save results
    pd.DataFrame(results, columns=["cluster", "type", "size", "term", "p-val", "FC", "desc"]).to_csv(
        f"{output}_go_enrichment_results.tsv", sep="\t", index=False
    )
        