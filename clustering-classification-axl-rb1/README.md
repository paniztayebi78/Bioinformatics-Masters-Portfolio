# Clustering and Classification of AXL and RB1 Gene Sequences: A Comparative Analysis

## Overview
This project investigates the clustering behavior of mRNA sequences from two critical genes with distinct roles in cancer biology: AXL (Receptor Tyrosine Kinase) and RB1 (Cell Cycle Regulator). Using hierarchical clustering and k-mer-based distance matrices, the study reveals fundamental differences in sequence conservation and variability between these genes.

## Biological Context
- **AXL Gene**: Member of the TAM receptor tyrosine kinase family, involved in cell proliferation, migration, and survival. Dysregulated in lung, breast, and pancreatic cancers.
- **RB1 Gene**: Encodes retinoblastoma protein (pRb), a tumor suppressor that regulates cell cycle progression and prevents excessive cell growth.

## Key Findings
- **AXL**: 304 loosely cohesive clusters (mean silhouette = 0.216) - high sequence heterogeneity
- **RB1**: 69 well-separated clusters (mean silhouette = 0.659) - conserved sequence structure
- Results align with biological roles: AXL's functional versatility vs. RB1's conserved tumor-suppressor function

## Technologies Used
- **R** - Primary analysis language
- **LaTeX** - Report formatting
- **Bioconductor packages**: Biostrings, DECIPHER
- **Statistical packages**: cluster, clValid
- **Visualization**: ggplot2, viridis, patchwork, ape

## Dataset
- **Source**: NCBI Nucleotide database (January 13, 2024)
- **AXL**: 358 sequences trimmed to 274 bp
- **RB1**: 161 sequences trimmed to 146 bp
- **Sequence types**: mRNA, coding sequences (CDS), and genomic sequences

## Methods
1. **Data Collection**: Sequences retrieved using NCBI search queries for human AXL and RB1 genes
2. **Quality Control**: Sequences trimmed to uniform lengths for comparison
3. **Distance Calculation**: K-mer-based distance matrices using DECIPHER
4. **Clustering**: Hierarchical clustering with average linkage (height cutoff = 0.2)
5. **Validation**: Silhouette analysis and Dunn indices for cluster quality assessment
6. **Visualization**: PCA plots, dendrograms, and distance distributions

## Results Summary
| Gene | Clusters | Mean Silhouette | Dunn Index | Davies-Bouldin Index | Mean Distance |
|------|----------|----------------|------------|---------------------|---------------|
| AXL  | 304      | 0.216          | 3.952      | 0.02                | 0.35 ± 0.12   |
| RB1  | 69       | 0.659          | 1.176      | 0.082               | 0.18 ± 0.07   |

## Files in this Repository
- `clustering_analysis.R` - Main R script for sequence analysis
- `AXL_RB1_Clustering_Report.pdf` - Complete research report
- `data/` - Processed sequence data files
- `figures/` - Generated plots and visualizations
- `results/` - Clustering outputs and validation metrics

## Usage
```r
# Load required libraries
library(Biostrings)
library(DECIPHER)
library(cluster)
library(ggplot2)

# Run the main analysis script
source("clustering_analysis.R")
```

## Biological Implications
The contrasting clustering patterns suggest that:
- **AXL's sequence heterogeneity** reflects its diverse cellular functions (immune response, metastasis)
- **RB1's sequence conservation** indicates strict evolutionary constraints due to its critical tumor-suppressor role

## Future Directions
- Incorporate protein-level alignment analysis
- Expand to include orthologs from other species
- Apply alternative distance metrics for validation
- Investigate functional implications of identified sequence variants

## Author
Paniz Tayebi, University of Guelph

## References
See full paper for complete bibliography of bioinformatics tools and methods used.
