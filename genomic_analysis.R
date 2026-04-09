## ============================================================================
## Genomic Analysis: PCA & LD Plots with Genotype Group Coloring
## ============================================================================
## Input : Genomic marker data in 0-1-2 format (rows = samples, cols = markers)
##         First column must be sample ID (e.g. Y24_001 … Y24_370)
## Output: PCA scatter plot and LD heatmap, colored / annotated by group
## ============================================================================

# ---- 0. Install / load required packages ------------------------------------

required_packages <- c("ggplot2", "genetics", "LDheatmap", "RColorBrewer",
                        "ggfortify", "reshape2", "grid")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

library(ggplot2)
library(reshape2)

# ---- 1. Read the 0-1-2 genotype data ----------------------------------------
# Adjust the file name / path to match your actual data file.
# Expected layout:
#   SampleID  SNP1  SNP2  SNP3  …
#   Y24_001     0     1     2   …
#   Y24_002     1     0     0   …
#   …

data_file <- "genotype_data.csv"   # <<< CHANGE THIS to your file name

if (!file.exists(data_file)) {
  stop(paste0(
    "Data file '", data_file, "' not found.\n",
    "Please place your 0-1-2 genotype file in the project directory\n",
    "and update the 'data_file' variable in this script."
  ))
}

# Read data – adjust sep if your file uses tab or space delimiters
geno_raw <- read.csv(data_file, header = TRUE, check.names = FALSE)

# First column = sample IDs; remaining columns = markers
sample_ids <- as.character(geno_raw[, 1])
geno_matrix <- as.matrix(geno_raw[, -1])
rownames(geno_matrix) <- sample_ids

cat("Loaded", nrow(geno_matrix), "samples and",
    ncol(geno_matrix), "markers.\n")

# ---- 2. Assign genotype groups & colours ------------------------------------

# Extract numeric part of each sample ID (e.g. "Y24_042" -> 42)
sample_numbers <- as.integer(sub(".*_(\\d+)$", "\\1", sample_ids))

# Define groups
assign_group <- function(num) {
  if (num >= 1 & num <= 326)   return("Winter Wheat")
  if (num >= 327 & num <= 357) return("Spring Wheat")
  if (num >= 358 & num <= 370) return("Spelt")
  return("Unknown")
}

groups <- sapply(sample_numbers, assign_group)

# High-contrast, colour-blind-friendly palette
group_colors <- c(
  "Winter Wheat"  = "#E63946",   # vivid red-coral

  "Spring Wheat"  = "#2A9D8F",   # teal-green
  "Spelt"         = "#457B9D"    # steel-blue
)

cat("Group counts:\n")
print(table(groups))

# ---- 3. Handle missing data -------------------------------------------------
# Replace any NA with the column (marker) mean so PCA can run
na_count <- sum(is.na(geno_matrix))
if (na_count > 0) {
  cat("Imputing", na_count, "missing values with column means.\n")
  col_means <- colMeans(geno_matrix, na.rm = TRUE)
  for (j in seq_len(ncol(geno_matrix))) {
    geno_matrix[is.na(geno_matrix[, j]), j] <- col_means[j]
  }
}

# Remove markers with zero variance (uninformative for PCA / LD)
marker_vars <- apply(geno_matrix, 2, var)
keep <- marker_vars > 0
cat("Keeping", sum(keep), "of", length(keep),
    "markers (removed", sum(!keep), "zero-variance markers).\n")
geno_matrix <- geno_matrix[, keep]

# ---- 4. PCA -----------------------------------------------------------------

pca_result <- prcomp(geno_matrix, center = TRUE, scale. = TRUE)

# Proportion of variance explained
pve <- summary(pca_result)$importance[2, ]   # row 2 = Proportion of Variance

pca_df <- data.frame(
  PC1   = pca_result$x[, 1],
  PC2   = pca_result$x[, 2],
  Group = groups
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_colour_manual(values = group_colors) +
  labs(
    title = "PCA of Genotype Data (0-1-2 Format)",
    x = paste0("PC1 (", round(pve[1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(pve[2] * 100, 1), "% variance)"),
    colour = "Genotype Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(pca_plot)
ggsave("PCA_plot.png", pca_plot, width = 9, height = 7, dpi = 300)
cat("PCA plot saved to PCA_plot.png\n")

# Optional: scree plot
scree_df <- data.frame(
  PC       = seq_along(pve),
  Variance = pve * 100
)
scree_df <- scree_df[1:min(20, nrow(scree_df)), ]   # show first 20 PCs

scree_plot <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_col(fill = "#264653") +
  geom_line(colour = "#E76F51", linewidth = 0.8) +
  geom_point(colour = "#E76F51", size = 2) +
  labs(
    title = "Scree Plot – Variance Explained by Each PC",
    x = "Principal Component",
    y = "% Variance Explained"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(scree_plot)
ggsave("Scree_plot.png", scree_plot, width = 8, height = 5, dpi = 300)
cat("Scree plot saved to Scree_plot.png\n")

# ---- 5. LD (Linkage Disequilibrium) -----------------------------------------
# Computing LD on all markers can be very slow.
# We'll compute pairwise r² on a subset of markers.

max_markers_ld <- 100   # adjust for speed vs. resolution

if (ncol(geno_matrix) > max_markers_ld) {
  set.seed(42)
  ld_idx <- sort(sample(seq_len(ncol(geno_matrix)), max_markers_ld))
  geno_ld <- geno_matrix[, ld_idx]
  cat("Using a random subset of", max_markers_ld,
      "markers for the LD heatmap.\n")
} else {
  geno_ld <- geno_matrix
}

# Compute pairwise r² matrix
n_markers <- ncol(geno_ld)
r2_matrix <- matrix(NA, n_markers, n_markers)

for (i in seq_len(n_markers)) {
  for (j in i:n_markers) {
    if (i == j) {
      r2_matrix[i, j] <- 1
    } else {
      ct <- cor(geno_ld[, i], geno_ld[, j], use = "pairwise.complete.obs")
      r2_matrix[i, j] <- ct^2
      r2_matrix[j, i] <- ct^2
    }
  }
}

colnames(r2_matrix) <- colnames(geno_ld)
rownames(r2_matrix) <- colnames(geno_ld)

# Melt for ggplot
ld_melted <- melt(r2_matrix)
colnames(ld_melted) <- c("Marker1", "Marker2", "r2")

ld_plot <- ggplot(ld_melted, aes(x = Marker1, y = Marker2, fill = r2)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("#F1FAEE", "#A8DADC", "#457B9D", "#1D3557"),
    limits  = c(0, 1),
    name    = expression(r^2)
  ) +
  labs(
    title = "Linkage Disequilibrium (LD) Heatmap"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks   = element_blank(),
    axis.title   = element_blank(),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

print(ld_plot)
ggsave("LD_heatmap.png", ld_plot, width = 8, height = 7, dpi = 300)
cat("LD heatmap saved to LD_heatmap.png\n")

# ---- 6. LD Decay Plot (optional) --------------------------------------------
# If marker positions (bp) are available you can plot LD decay.
# Below is a template assuming columns are named "chr_pos" or similar.

# Uncomment and adapt if you have position information:
# positions <- as.integer(sub(".*_(\\d+)$", "\\1", colnames(geno_ld)))
# dist_mat  <- abs(outer(positions, positions, "-"))
# dist_vec  <- dist_mat[upper.tri(dist_mat)]
# r2_vec    <- r2_matrix[upper.tri(r2_matrix)]
#
# ld_decay_df <- data.frame(Distance = dist_vec, r2 = r2_vec)
#
# ld_decay_plot <- ggplot(ld_decay_df, aes(x = Distance, y = r2)) +
#   geom_point(alpha = 0.15, size = 0.6, colour = "#457B9D") +
#   geom_smooth(method = "loess", colour = "#E63946", se = FALSE) +
#   labs(title = "LD Decay", x = "Distance (bp)", y = expression(r^2)) +
#   theme_minimal(base_size = 14) +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#
# print(ld_decay_plot)
# ggsave("LD_decay.png", ld_decay_plot, width = 8, height = 5, dpi = 300)

cat("\n=== Analysis complete ===\n")
cat("Output files:\n")
cat("  - PCA_plot.png\n")
cat("  - Scree_plot.png\n")
cat("  - LD_heatmap.png\n")
