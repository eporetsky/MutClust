# MutClust: Mutual Rank-Based Clustering and GO Enrichment Analysis

**MutClust** is a Python package designed for RNA-seq gene coexpression analyses. It performs mutual rank (MR)-based clustering of coexpressed genes and identifies enriched Gene Ontology (GO) terms for the resulting clusters. The package is optimized for speed, able to run a whole-genome coexpression analysis in minutes.

---

## Features

- **Mutual Rank Analysis**: Calculates MR from Pearson correlation coefficients to identify coexpressed genes.
- **Leiden Clustering**: Groups genes into clusters based on mutual rank and exponential decay weights.
- **Gene Annotations**: Merge cluster members with gene annotations, if provided.
- **GO Enrichment Analysis**: Identifies enriched GO terms for each cluster using GOATOOLS.
- **Highly Configurable**: Supports adjustable thresholds, resolution parameters, and multi-threading for performance optimization.

---

## Installation

TODO: You can install MutClust directly from PyPI:

```bash
pip install mutclust
```

Alternatively, you can clone the repository and install it locally:

```bash
git clone https://github.com/eporetsky/mutclust.git
cd mutclust
pip install .
```

---

## Usage

MutClust provides a command-line interface (CLI) for running the full pipeline. After installation, you can use the `mutclust` command.

### Command-Line Arguments

| Argument              | Short | Description                                              | Default       |
|-----------------------|-------|----------------------------------------------------------|---------------|
| `--expression`        | `-ex` | Path to the RNA-seq dataset (TSV format).                | **-ex or -mr required**  |
| `--mutual_rank`       | `-mr` | Path to Mutual Rank file (TSV format).                   | **-ex or -mr required**  |
| `--annotations`       | `-a`  | Path to the gene annotation file.                        | **Optional**  |
| `--go_obo`            | `-go` | Path to the Gene Ontology (GO) OBO file.                 | **Optional**  |
| `--go_gaf`            | `-gf` | Path to the GO annotation file (GAF format).             | **Optional**  |
| `--output`            | `-o`  | Output prefix for the results.                           | **Required**  |
| `--mr_threshold`      | `-m`  | Mutual rank threshold for filtering.                     | `100`         |
| `--e_value`           | `-e`  | Exponential decay constant.                              | `10`          |
| `--resolution`        | `-r`  | Resolution parameter for Leiden clustering.              | `0.1`         |
| `--threads`           | `-t`  | Number of threads for correlation calculation.           | `4`           |
| `--save_intermediate` | `-t`  | Number of threads for correlation calculation.           | **Optional**  |

### Example Command

```bash
mutclust --expression data/AtCol-0.cpm.tsv \
         --annotations annotations/AtCol-0.annot.tsv \
         --go_obo go-basic.obo \
         --go_gaf tair.gaf \
         --output results/mutclust_output \
         --mr_threshold 100 \
         --e_value 10 \
         --resolution 0.1 \
         --threads 8
```

---

## Input File Formats

### RNA-seq Dataset
- **Format**: Tab-separated values (TSV).
- **Columns**: Gene IDs as row indices and samples as columns.
- **Example**:
```tsv
geneID    Sample1    Sample2    Sample3
GeneA     1.23       2.34       3.45
GeneB     4.56       5.67       6.78
```

### Gene Annotation File
- **Format**: Tab-separated values (TSV).
- **Columns**: `geneID` and additional annotation fields. 
- **Example**:
```tsv
geneID    description
GeneA     Photosynthesis-related protein
GeneB     Transcription factor
```

### GO OBO File
- **Description**: The Gene Ontology (GO) OBO file contains the ontology structure.
- **Source**: Download from [Gene Ontology](http://geneontology.org/).

### GO GAF File
- **Description**: The Gene Annotation File (GAF) maps genes to GO terms.
- **Source**: Download from [Gene Ontology](http://geneontology.org/).

---

## Output Files

1. **Filtered MR and e-values** (`<output_prefix>.mrs.tsv`):
   - Lists of coexpressed genes with MR and e-values.
   - Columns: `cluster_id`, `geneID`.

   **Example**:
   ```tsv
   Gene1    Gene2    MR    ED
   GeneA    GeneB    10.2  0.39
   GeneB    GeneC    6     0.6
   ```

2. **Clustered Genes** (`<output_prefix>.clusters.tsv`):
   - Lists genes in each cluster.
   - Annotation columns if provided.
   - Columns: `cluster_id`, `geneID`.

   **Example**:
   ```tsv
   cluster_id    geneID    Annotations
   1             GeneA     ...
   1             GeneB     ...
   ```

3. **GO Enrichment Results** (`<output_prefix>_go_enrichment_results.tsv`):
   - Contains enriched GO terms for each cluster.
   - Columns: `cluster`, `type`, `size`, `term`, `p-val`, `FC`, `desc`.

   **Example**:
   ```tsv
   cluster    type    size    term       p-val       FC    desc
   1          BP      25      GO:0008150 0.00123     3.5   Biological Process
   ```
---

## Dependencies

The following Python libraries are required and will be installed automatically:
- `numpy`
- `pandas`
- `pynetcor`
- `python-igraph`
- `goatools`

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Contributing

Contributions, suggestions and issues are welcome!