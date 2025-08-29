###############################################################################
# mRNA Sequence Clustering Pipeline
# Purpose: Compare AXL (Receptor Tyrosine Kinase) vs RB1 (Cell Cycle Regulator) genes
# Input: FASTA files containing mRNA sequences for both genes
# Output: RData file with all results + visualization plots
###############################################################################

# ---------------------------
# 1. PACKAGE INSTALLATION
# ---------------------------

# Check if BiocManager is installed (required for bioinformatics packages)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Define all required packages:
# - Biostrings/DECIPHER: For biological sequence analysis
# - ggplot2/patchwork: For visualization
# - ape: For phylogenetic tree analysis
# - cluster/clValid: For clustering validation
required_packages <- c("Biostrings", "DECIPHER", "ggplot2", "ape", 
                       "patchwork", "cluster", "clValid")

# Install any missing packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("Biostrings", "DECIPHER")) {
      BiocManager::install(pkg)  # Install bioconductor packages
    } else {
      install.packages(pkg)  # Install CRAN packages
    }
  }
}

# ---------------------------
# 2. LIBRARY LOADING
# ---------------------------

# Load all required libraries
library(Biostrings)  # For handling biological sequences
library(DECIPHER)    # For sequence alignment and analysis
library(ggplot2)     # For creating publication-quality plots
library(ape)         # For phylogenetic tree manipulation
library(patchwork)   # For combining ggplot2 plots
library(cluster)     # For clustering algorithms
library(clValid)     # For cluster validation metrics

# ---------------------------
# 3. CORE ANALYSIS FUNCTIONS
# ---------------------------

analyze_gene <- function(seqs, gene_name) {
  # Purpose: Process and cluster sequences for a single gene
  # Input: DNAStringSet of sequences, gene name string
  # Output: List containing all analysis results
  
  # Filter for mRNA sequences only (ignore other types)
  is_mRNA <- grepl("mRNA|cds|transcript", names(seqs), ignore.case = TRUE)
  seqs <- seqs[is_mRNA]
  message("Retained ", length(seqs), " mRNA sequences for ", gene_name)
  
  # Remove length outliers (top/bottom 1%) to avoid extreme sequences
  len <- width(seqs)
  seqs <- seqs[len > quantile(len, 0.01) & len < quantile(len, 0.99)]
  
  # Trim all sequences to common minimum length for fair comparison
  min_len <- min(width(seqs))
  seqs <- subseq(seqs, 1, min_len)
  message("Trimmed all sequences to ", min_len, " bp for consistent comparison")
  
  # Calculate pairwise distance matrix using DECIPHER's optimized algorithm
  message("Calculating distance matrix...")
  dist_mat <- DistanceMatrix(seqs, processors = NULL)  # NULL uses all available cores
  
  # Hierarchical clustering using average linkage method
  hclust <- hclust(as.dist(dist_mat), method = "average")
  clusters <- cutree(hclust, h = 0.2)  # Cut tree at height 0.2 to define clusters
  
  # Calculate cluster validation metrics
  message("Calculating cluster validity indices...")
  sil <- silhouette(clusters, as.dist(dist_mat))  # Measures cluster cohesion/separation
  dunn <- dunn(as.dist(dist_mat), clusters)      # Measures inter-cluster distance
  
  # Create PCA plot for visualization
  pca <- prcomp(dist_mat)$x  # Perform PCA on distance matrix
  pca_plot <- ggplot(data.frame(pca), aes(PC1, PC2, color = factor(clusters))) +
    geom_point(alpha = 0.7, size = 2) +
    ggtitle(paste(gene_name, "Clusters (n =", length(seqs), ")")) +
    theme_minimal(base_size = 12) +
    scale_color_viridis_d(guide = "none") +  # Color by cluster, no legend
    theme(plot.margin = margin(5, 5, 5, 5, "pt"))
  
  # Function to plot dendrogram
  dendro_plot <- function() {
    par(mar = c(1, 4, 2, 1))  # Set plot margins
    plot(hclust, labels = FALSE, main = paste(gene_name, "Dendrogram"), 
         xlab = "", sub = "", cex.main = 0.9)
    rect.hclust(hclust, h = 0.2, border = "red")  # Highlight clusters
  }
  
  # Return all results in a structured list
  return(list(
    sequences = seqs,      # Filtered sequences
    clusters = clusters,   # Cluster assignments
    dist_mat = dist_mat,   # Distance matrix
    pca_plot = pca_plot,   # PCA visualization
    dendro_plot = dendro_plot,  # Dendrogram function
    silhouette = sil,      # Silhouette scores
    dunn_index = dunn,     # Dunn index
    hclust = hclust        # Hierarchical clustering object
  ))
}

compare_genes <- function(axl, rb1) {
  # Purpose: Compare results between AXL and RB1 genes
  # Input: Analysis results for both genes
  # Output: Comparison metrics and combined plots
  
  # Print cluster count comparison to console
  cat("\n===== Gene Comparison Results =====\n")
  cat("AXL clusters:", length(unique(axl$clusters)), "\n")
  cat("RB1 clusters:", length(unique(rb1$clusters)), "\n\n")
  
  # Print cluster quality metrics
  cat("Cluster Quality Indices:\n")
  cat("AXL - Mean Silhouette Width:", mean(axl$silhouette[, 3]), "\n")
  cat("AXL - Dunn Index:", axl$dunn_index, "\n")
  cat("RB1 - Mean Silhouette Width:", mean(rb1$silhouette[, 3]), "\n")
  cat("RB1 - Dunn Index:", rb1$dunn_index, "\n")
  
  # Create combined distance distribution plot
  dist_df <- rbind(
    data.frame(Gene = "AXL", Dist = as.vector(axl$dist_mat)),
    data.frame(Gene = "RB1", Dist = as.vector(rb1$dist_mat))
  )
  
  dist_plot <- ggplot(dist_df, aes(Dist, fill = Gene)) +
    geom_density(alpha = 0.5) +
    ggtitle("Sequence Distance Distributions") +
    theme_minimal(base_size = 12) +
    scale_fill_manual(values = c("#440154", "#21908C")) +  # Purple-teal color scheme
    theme(plot.margin = margin(5, 5, 5, 5, "pt"))
  
  # Combine PCA plots using patchwork
  combined_plots <- (axl$pca_plot + rb1$pca_plot) / 
    dist_plot +
    plot_layout(heights = c(1, 0.7))  # Relative heights of plot rows
  
  # Save combined visualization
  ggsave("combined_plots.png", combined_plots, 
         width = 8, height = 6, dpi = 300, bg = "white")
  
  # Save dendrograms side-by-side
  png("dendrograms.png", width = 8, height = 4, units = "in", res = 300)
  par(mfrow = c(1, 2), mar = c(1, 3, 2, 1))  # Set up 1x2 plot grid
  axl$dendro_plot()
  rb1$dendro_plot()
  dev.off()
  
  # Return comparison results
  return(list(
    cluster_counts = c(AXL = length(unique(axl$clusters)),
                       RB1 = length(unique(rb1$clusters))),
    validity_metrics = data.frame(
      Gene = c("AXL", "RB1"),
      Silhouette = c(mean(axl$silhouette[, 3]), mean(rb1$silhouette[, 3])),
      Dunn = c(axl$dunn_index, rb1$dunn_index)
    ),
    distance_plot = dist_plot
  ))
}

# ---------------------------
# 4. MAIN WORKFLOW
# ---------------------------

# Define input files (update paths as needed)
fasta_files <- c(
  AXL = "/Users/paniztayebi/Library/Mobile Documents/com~apple~CloudDocs/Documents/University of Guelph/CIS*6060/Project/AXL.fasta",  
  RB1 = "/Users/paniztayebi/Library/Mobile Documents/com~apple~CloudDocs/Documents/University of Guelph/CIS*6060/Project/Rb1.fasta"  
)

# Verify files exist before proceeding
for (f in fasta_files) {
  if (!file.exists(f)) stop("File not found: ", f)
}

# Load sequences using Biostrings
message("\nLoading sequence data...")
sequences <- list(
  AXL = readDNAStringSet(fasta_files["AXL"]),  # Load AXL sequences
  RB1 = readDNAStringSet(fasta_files["RB1"])   # Load RB1 sequences
)

# Perform analysis on both genes
message("\nRunning analysis...")
results <- list(
  AXL = analyze_gene(sequences$AXL, "AXL"),  # Analyze AXL
  RB1 = analyze_gene(sequences$RB1, "RB1")   # Analyze RB1
)

# Compare results between genes
comparison <- compare_genes(results$AXL, results$RB1)

# Save all results to RData file (preserves R objects exactly)
save(results, comparison, file = "gene_clustering_results.RData")

# Final status message
message("\nAnalysis complete! Results saved to:")
message("- gene_clustering_results.RData (all R objects)")
message("- combined_plots.png (PCA and distance plots)")
message("- dendrograms.png (cluster dendrograms)")