# MutClust: Mutual Rank-Based Clustering and GO Enrichment Analysis

**MutClust** is a Python package designed for RNA-seq gene coexpression analyses. It performs mutual rank (MR)-based clustering of coexpressed genes and identifies enriched Gene Ontology (GO) terms for the resulting clusters. The package is optimized for speed, able to run a whole-genome coexpression analysis in minutes.

---

## Features

- **Mutual Rank Analysis**: Calculates MR from Pearson correlation coefficients to identify coexpressed genes.
- **Leiden Clustering**: Groups genes into clusters based on mutual rank and exponential decay weights.
- **Gene Annotations**: Merge cluster members with gene annotations, if provided.
- **GO Enrichment Analysis**: Identifies enriched GO terms for each cluster using GOATOOLS.
- **Highly Configurable**: Supports adjustable thresholds, resolution parameters, and multi-threading for performance optimization.
- **Calculate correlation matrix and mutual rank from RNA-seq data**
- **Filter and apply exponential decay to mutual rank values**
- **Perform Leiden clustering to identify co-expressed gene clusters**
- **Calculate eigen-genes for each cluster (first principal component)**
- **Perform GO enrichment analysis on gene clusters**
- **Annotate clusters with gene information**

---

## Installation

You can install MutClust directly from PyPI:

```bash
pip install mutclust
```

Note: Because of a known [dependency issue](https://github.com/01life/pyNetCor/issues/1) with PyNetCor, MutClust is not currently available on MacOS through PyPI but installs properly on Linux. 

Alternatively, you can clone the repository and install it locally:

```bash
git clone https://github.com/eporetsky/mutclust.git
cd mutclust
pip install .
```

### Docker Installation

For users who prefer containerized deployment, MutClust is available as a Docker container:

```bash
# Build the container
docker build -t mutclust .

# Run MutClust with your data
docker run -v /path/to/your/data:/data mutclust --expression /data/your_expression.tsv --output /data/results
```

The container uses Ubuntu 20.04 and includes all necessary dependencies. Mount your data directory to `/data` inside the container to access your files.

---

## Usage

MutClust provides a command-line interface (CLI) for running the full pipeline. After installation, you can use the `mutclust` command.

### Basic Usage

```bash
# Using expression data
mutclust --expression input.tsv --output output_prefix

# Using pre-calculated mutual rank
mutclust --mutual_rank input.mr.tsv --output output_prefix
```

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
| `--save_intermediate` | -     | Save intermediate files (PCC, MR, filtered pairs).       | **Optional**  |
| `--eigengene`         | -     | Calculate eigen-genes for clusters.                      | `True`        |

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

4. **Eigen-gene values** (`<output_prefix>.eigen.tsv`):
   - Eigen-gene values for each cluster.
   - Columns: `geneID` and sample columns.

   **Example**:
   ```tsv
   geneID    Sample1    Sample2    Sample3
   c1        0.707107   0.707107   0.707107
   c2        0.577350   0.577350   0.577350
   c3        0.500000   0.500000   0.500000
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