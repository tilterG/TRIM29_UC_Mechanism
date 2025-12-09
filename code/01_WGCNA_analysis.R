# ============================================================================
# Script: 01_WGCNA_analysis.R
# Project: TRIM29_UC_Mechanism
# Purpose: Weighted Gene Co-expression Network Analysis (WGCNA)
# Steps: 1. Data preparation & filtering
#        2. Sample clustering & outlier detection
#        3. Soft-thresholding power selection
#        4. Network construction & module detection
#        5. Module-trait correlation analysis
#        6. Visualization
# Input:  - WGCNA.fpkm: Expression matrix (genes x samples)
#         - allTraits: Clinical/phenotypic trait data
# Output: - Filtered expression matrix
#         - Module assignment & module-trait correlation results
#         - Diagnostic plots (PDF)
# ============================================================================

# ----------------------------
# 1. Environment Setup
# ----------------------------
# Load required package
library(WGCNA)

# Configure environment
options(stringsAsFactors = FALSE)
enableWGCNAThreads()  # Enable multi-threading

# ----------------------------
# 2. Data Loading & Preparation
# ----------------------------
# Extract expression data and transpose
datExpr0 <- as.data.frame(t(WGCNA.fpkm[, -1]))
colnames(datExpr0) <- WGCNA.fpkm$sample
rownames(datExpr0) <- colnames(WGCNA.fpkm[, -1])

# Check data for missing values and filter
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  # Remove problematic genes and samples
  if (sum(!gsg$goodGenes) > 0) {
    cat("Removing genes:", paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", "), "\n")
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", "), "\n")
  }
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# ----------------------------
# 3. Expression Level Filtering
# ----------------------------
meanFPKM_threshold <- 0.5  # Adjustable threshold
n_samples <- nrow(datExpr0)

# Calculate mean expression per gene
datExpr0[n_samples + 1, ] <- apply(datExpr0[1:n_samples, ], 2, mean)

# Retain genes with mean expression above threshold
datExpr0_filtered <- datExpr0[1:n_samples, datExpr0[n_samples + 1, ] > meanFPKM_threshold]

# Save filtered expression matrix
filtered_fpkm <- data.frame(
  sample = rownames(t(datExpr0_filtered)),
  t(datExpr0_filtered),
  check.names = FALSE
)
write.table(filtered_fpkm,
            file = "results/WGCNA_filtered_expression_matrix.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Update data variable for downstream analysis
datExpr1 <- datExpr0_filtered

# ----------------------------
# 4. Sample Clustering & Outlier Detection
# ----------------------------
sampleTree <- hclust(dist(datExpr1), method = "average")

pdf("results/figures/WGCNA_01_sample_clustering.pdf", width = 22, height = 8)
par(cex = 0.6, mar = c(0, 6, 6, 0))
plot(sampleTree,
     main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "",
     cex.lab = 2, cex.axis = 1.5, cex.main = 2)
dev.off()

# Optional: Automatic outlier cut (adjust height as needed)
# clust <- cutreeStatic(sampleTree, cutHeight = 87, minSize = 10)
# table(clust)

# ----------------------------
# 5. Trait Data Integration
# ----------------------------
# Ensure sample matching between expression and trait data
fpkmSamples <- rownames(datExpr1)
traitSamples <- rownames(allTraits)
traitRows <- match(fpkmSamples, traitSamples)
datTraits <- as.data.frame(allTraits[traitRows, ])
rownames(datTraits) <- fpkmSamples

# Generate sample dendrogram with trait heatmap
sampleTree2 <- hclust(dist(datExpr1), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)

pdf("results/figures/WGCNA_02_sample_dendrogram_trait_heatmap.pdf",
    width = 45, height = 12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(datTraits),
                    main = "Sample Dendrogram and Trait Heatmap",
                    cex.colorLabels = 1.5,
                    cex.dendroLabels = 1,
                    cex.rowText = 2)
dev.off()

# ----------------------------
# 6. Soft-Thresholding Power Selection
# ----------------------------
powers <- c(1:10, seq(from = 5, to = 30, by = 2))
soft_power_results <- pickSoftThreshold(datExpr1,
                                        powerVector = powers,
                                        verbose = 5)

# Plot scale independence and mean connectivity
pdf("results/figures/WGCNA_03_soft_threshold_selection.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale independence plot
plot(soft_power_results$fitIndices[, 1],
     -sign(soft_power_results$fitIndices[, 3]) * soft_power_results$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R²)",
     type = "n",
     main = "Scale Independence")
text(soft_power_results$fitIndices[, 1],
     -sign(soft_power_results$fitIndices[, 3]) * soft_power_results$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red", lty = 2)

# Mean connectivity plot
plot(soft_power_results$fitIndices[, 1],
     soft_power_results$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean Connectivity")
text(soft_power_results$fitIndices[, 1],
     soft_power_results$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")
dev.off()

# Automated power selection
min_mean_conn <- 10
best_power <- NA

for (i in seq_along(powers)) {
  if (soft_power_results$fitIndices[i, 2] > 0.85 &&
      soft_power_results$fitIndices[i, 5] > min_mean_conn) {
    best_power <- powers[i]
    break
  }
}

if (is.na(best_power)) {
  stop("No suitable soft-thresholding power found. Please adjust selection criteria.")
} else {
  cat(sprintf("Selected soft-thresholding power (beta): %d\n", best_power))
}

# ----------------------------
# 7. Network Construction & Module Detection
# ----------------------------
# Calculate adjacency and TOM
adjacency_matrix <- adjacency(datExpr1, power = best_power)
TOM_similarity <- TOMsimilarity(adjacency_matrix)
dissTOM <- 1 - TOM_similarity

# Gene clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")

pdf("results/figures/WGCNA_04_gene_clustering_dendrogram.pdf",
    width = 24, height = 18)
plot(geneTree,
     xlab = "", sub = "",
     main = "Gene Clustering on TOM-based Dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# Dynamic module detection
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree,
                             distM = dissTOM,
                             deepSplit = 2,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)

pdf("results/figures/WGCNA_05_dynamic_tree_cut.pdf", width = 15, height = 6)
plotDendroAndColors(geneTree, dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
dev.off()

# ----------------------------
# 8. Module Merging
# ----------------------------
# Calculate module eigengenes and cluster
MEList <- moduleEigengenes(datExpr1, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot eigengene clustering
pdf("results/figures/WGCNA_06_eigengene_clustering.pdf",
    width = 10, height = 6)
plot(METree,
     main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
MEDissThres <- 0.4  # Adjustable merge threshold
abline(h = MEDissThres, col = "red")
dev.off()

# Merge similar modules
merged <- mergeCloseModules(datExpr1,
                            dynamicColors,
                            cutHeight = MEDissThres,
                            verbose = 3)
mergedColors <- merged$colors
mergedMEs <- merged$newMEs

pdf("results/figures/WGCNA_07_merged_modules.pdf",
    width = 9, height = 6)
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged Dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()

# ----------------------------
# 9. Module-Trait Correlation Analysis
# ----------------------------
moduleColors <- mergedColors
MEs <- mergedMEs

nGenes <- ncol(datExpr1)
nSamples <- nrow(datExpr1)

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Save correlation results
write.csv(moduleTraitCor,
          "results/WGCNA_module_trait_correlation.csv",
          row.names = TRUE)
write.csv(moduleTraitPvalue,
          "results/WGCNA_module_trait_pvalue.csv",
          row.names = TRUE)

# Visualize module-trait relationships
pdf("results/figures/WGCNA_08_module_trait_relationships.pdf",
    width = 10, height = 10)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()

# ----------------------------
# 10. Network Visualization (Optional)
# ----------------------------
# Subset genes for TOM plot visualization
set.seed(10)  # For reproducibility
nSelect <- 400
selectGenes <- sample(nGenes, size = nSelect)
selectTOM <- dissTOM[selectGenes, selectGenes]
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- moduleColors[selectGenes]

plotDiss <- selectTOM^7
diag(plotDiss) <- NA

# Requires 'gplots' for colorpanel function
if (!require("gplots", quietly = TRUE)) {
  install.packages("gplots")
}
library(gplots)

pdf("results/figures/WGCNA_09_network_heatmap_selected_genes.pdf",
    width = 9, height = 9)
mycol <- colorpanel(250, 'red', 'orange', 'lemonchiffon')
TOMplot(plotDiss, selectTree, selectColors,
        col = mycol,
        main = "Network Heatmap Plot (Selected Genes)")
dev.off()

# ----------------------------
# 11. Eigengene Network Visualization
# ----------------------------
pdf("results/figures/WGCNA_10_eigengene_network.pdf",
    width = 5, height = 7.5)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "",
                      marDendro = c(0, 4, 1, 2),
                      marHeatmap = c(3, 4, 1, 2),
                      cex.lab = 0.8,
                      xLabelsAngle = 90)
dev.off()

# Separate plots for dendrogram and heatmap
pdf("results/figures/WGCNA_11_eigengene_dendrogram.pdf",
    width = 6, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MEs,
                      "Eigengene Dendrogram",
                      marDendro = c(0, 4, 2, 0),
                      plotHeatmaps = FALSE)
dev.off()

pdf("results/figures/WGCNA_12_eigengene_adjacency_heatmap.pdf",
    width = 6, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MEs,
                      "Eigengene Adjacency Heatmap",
                      marHeatmap = c(3, 4, 2, 2),
                      plotDendrograms = FALSE,
                      xLabelsAngle = 90)
dev.off()

# ----------------------------
# 12. Save Key R Objects
# ----------------------------
save(datExpr1, datTraits, best_power, geneTree,
     moduleColors, MEs, moduleTraitCor, moduleTraitPvalue,
     file = "results/WGCNA_analysis_results.RData")

cat("WGCNA analysis completed successfully.\n")
cat("Outputs saved to 'results/' directory.\n")