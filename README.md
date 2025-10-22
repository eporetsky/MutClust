# MutClust: Efficient and Scalable Mutual Rank-Based Coexpression Clustering

**MutClust** is a Python tool for efficient and scalable mutual rank-based gene coexpression analyses. The clustering analysis is conducted using [ClusterONE](https://paccanarolab.org/clusterone/), as described in [Wisecaver *et al.* 2017](https://academic.oup.com/plcell/article-abstract/29/5/944/6099316). MutClust is still under development.

---

## Features
- **Mutual Rank Analysis:** Compute mutual rank (MR) from Pearson correlations on your gene expression matrix.
- **ClusterONE Clustering:** Identify gene coexpression clusters from filtered/weighted MR networks.
- **Fast:** Multi-threaded, sparse matrix operations for speed on large datasets.

---

## Installation

### Recommended

**Install MutClust:** Create the recommended conda environment:
```bash
conda env create -f environment.yml
conda activate mutclust
```

### Alternative

**Step 1:** Make sure that ClusterONE is available from the command line:
```bash
conda install bioconda::clusterone
```

**Step 2a:** Install MutClust from PyPI:
```bash
pip install mutclust
```

**Step 2b:** Or clone the repository from GitHub:
```bash
git clone https://github.com/eporetsky/mutclust.git
cd mutclust
pip install .
```

## Usage

### 1. Calculate Mutual Rank (MR)

```bash
mutclust mr -i expr.tsv -o results.mrs.tsv.gz --mr-threshold 100 --threads 4 [--log2]
```

| Argument         | Short | Description                                         | Default       |
|------------------|-------|-----------------------------------------------------|---------------|
| --input          | -i    | Path to the RNA-seq dataset (.tsv/.tsv.gz)          | **Required**  |
| --output         | -o    | Output file for mutual rank pairs                   | **Required**  |
| --mr-threshold   | -m    | MR threshold for reporting gene pairs               | 100           |
| --threads        | -t    | Number of CPU threads (correlation)                 | 4             |
| --log2           |       | If set, applies log2(x+1) before calculation        | OFF by default|

- Input: Genes as rows, samples as columns (TSV, row index 'geneID').
- Output: Gzipped tab-separated file containing `Gene1`, `Gene2`, `MR`.

### 2. Cluster Genes (with ClusterONE)

```bash
mutclust cls -i results.mrs.tsv.gz -o results.cls.tsv --e_value 10
```

| Argument         | Short | Description                                         | Default       |
|------------------|-------|-----------------------------------------------------|---------------|
| --input          | -i    | Path to Mutual Rank (MR) pairs (.tsv/.tsv.gz)       | **Required**  |
| --output         | -o    | Output file for clusters (.tsv)                     | **Required**  |
| --e_value        | -e    | Exponential decay constant for edge weighting       | 10            |

- The tool filters/weights MR pairs and calls ClusterONE for clustering.
- Output: `clusters.tsv`, listing clusters with p-value < 0.1. Tab-separated file containing `clusterID`, `geneID`, `pval`.

---

## Example Workflow

```bash
mutclust mr -i data/myexpr.tsv -o out.mrs.tsv.gz --mr-threshold 100 --threads 72 --log2
mutclust cls -i out.mrs.tsv.gz -o out.clusters.tsv --e_value 10
```

---

## Input Format

Expression file:
```
geneID\tSample1\tSample2\n...
GeneA \t1.1    \t2.2
GeneB \t4.2    \t3.7
```

Note: MutClust might be limited to linux because of dependency on pynetcor.

---

## Coming Soon
- Generate cluster gene annotation
- Calculate cluster GO term enrichment
- Calculate clusteer eigen-gene data
- Add a MutClust Dockerfile
- Add unit testing

---

## License
MIT License. See LICENSE file for details.

---

## Contributing
Suggestions, pull requests, and issues welcome!