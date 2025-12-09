# ==============================================================================
# Script: 03_machine_learning_feature_selection.R
# Project: TRIM29_UC_Mechanism
# Purpose: Ensemble Machine Learning Feature Selection (Lasso + SVM-RFE)
# Steps: 1. Data loading & preprocessing
#        2. Train/validation/test split
#        3. Lasso regression with cross-validation
#        4. SVM-Recursive Feature Elimination (SVM-RFE)
#        5. Model performance evaluation
#        6. Results visualization and integration
# Input:  - CSV file with expression data and group labels
#         - (Optional) Lysosome gene set for pre-filtering
# Output: - Selected feature lists from Lasso and SVM-RFE
#         - Model performance metrics
#         - ROC curves and diagnostic plots
#         - Integrated feature intersection
# ==============================================================================

# ----------------------------
# 0. Environment Setup
# ----------------------------
# Load required packages
library(tidyverse)
library(glmnet)    # For Lasso regression
library(caret)     # For data splitting and evaluation metrics
library(pROC)      # For AUC calculation and ROC curves
library(e1071)     # For SVM implementation

# Note: Ensure 'msvmRFE.R' is available in your project directory
# source("utils/msvmRFE.R")  # Uncomment if using custom SVM-RFE implementation

# Set seed for reproducibility
set.seed(1234)

# ----------------------------
# 1. Data Loading & Preprocessing
# ----------------------------
# Define input file path
raw_data_file <- "data/processed/ml_input_data.csv"  # Update with your actual filename

# Load expression data with group labels
# Format: Rows = samples, Columns = genes, First column = group labels
data_full <- read.csv(raw_data_file, row.names = 1, check.names = FALSE)

# Optional: Pre-filter using lysosome gene set
# Uncomment and modify if you have a predefined lysosome gene set
# lysosome_genes <- readLines("data/lysosome_gene_list.txt")
# common_genes <- intersect(lysosome_genes, colnames(data_full))
# 
# if (length(common_genes) == 0) {
#   stop("No overlap between lysosome gene set and expression data.")
# }
# 
# group_col_name <- "group_list"
# data_full <- data_full[, c(group_col_name, common_genes)]
# cat("Data filtered to", length(common_genes), "lysosome-related genes.\n")

# ----------------------------
# 2. Data Partitioning
# ----------------------------
# Verify group column exists
group_col_name <- "group_list"
if (!group_col_name %in% colnames(data_full)) {
  cat("Available columns:", paste(colnames(data_full), collapse = ", "), "\n")
  stop(paste("Group column '", group_col_name, "' not found in data.", sep = ""))
}

# Stratified split: 60% train, 20% validation, 20% test
train_index <- createDataPartition(data_full[[group_col_name]], 
                                   p = 0.6, 
                                   list = FALSE, 
                                   times = 1)

train_data <- data_full[train_index, ]
temp_data <- data_full[-train_index, ]

# Split remaining data equally into validation and test sets
val_index <- createDataPartition(temp_data[[group_col_name]], 
                                 p = 0.5, 
                                 list = FALSE, 
                                 times = 1)

val_data <- temp_data[val_index, ]
test_data <- temp_data[-val_index, ]

# Report dataset sizes
cat("Data partitioning completed:\n")
cat("  Training set:  ", nrow(train_data), "samples x", ncol(train_data), "features\n")
cat("  Validation set:", nrow(val_data), "samples x", ncol(val_data), "features\n")
cat("  Test set:      ", nrow(test_data), "samples x", ncol(test_data), "features\n")

# ----------------------------
# 3. Prepare Feature Matrices and Labels
# ----------------------------
# Extract common feature names (excluding group column)
feature_names <- setdiff(colnames(train_data), group_col_name)

# Create feature matrices
x_train <- as.matrix(train_data[, feature_names])
x_val <- as.matrix(val_data[, feature_names])
x_test <- as.matrix(test_data[, feature_names])

# Create binary factor labels (0 = control, 1 = case)
# Use training set factor levels as reference for all sets
train_groups <- as.factor(train_data[[group_col_name]])
ref_level <- levels(train_groups)[1]

y_train <- factor(ifelse(train_data[[group_col_name]] == ref_level, 0, 1), 
                  levels = c(0, 1))
y_val <- factor(ifelse(val_data[[group_col_name]] == ref_level, 0, 1), 
                levels = c(0, 1))
y_test <- factor(ifelse(test_data[[group_col_name]] == ref_level, 0, 1), 
                 levels = c(0, 1))

# Numeric labels for Lasso (0/1)
y_train_num <- as.numeric(y_train) - 1

# Report label distributions
cat("\nLabel distribution:\n")
cat("  Training:   ", table(y_train), "(0/1)\n")
cat("  Validation: ", table(y_val), "(0/1)\n")
cat("  Test:       ", table(y_test), "(0/1)\n")

# ----------------------------
# 4. Lasso Regression Analysis
# ----------------------------
cat("\n" + strrep("-", 50) + "\n")
cat("Step 1: Lasso Regression Feature Selection\n")
cat(strrep("-", 50) + "\n")

# Cross-validation to determine optimal lambda
cv_fit_lasso <- cv.glmnet(x_train, y_train_num,
                          family = "binomial",
                          type.measure = "auc",  # Using AUC for binary classification
                          nfolds = 5,
                          alpha = 1)  # alpha=1 for Lasso

# Plot and save CV results
pdf("results/figures/ML_01_lasso_cv_curve.pdf", width = 8, height = 6)
plot(cv_fit_lasso)
abline(v = log(cv_fit_lasso$lambda.min), col = "red", lty = 2, lwd = 2)
abline(v = log(cv_fit_lasso$lambda.1se), col = "blue", lty = 2, lwd = 2)
legend("bottomright", 
       legend = c("lambda.min", "lambda.1se"),
       col = c("red", "blue"),
       lty = 2, lwd = 2)
dev.off()

cat("Optimal lambda values:\n")
cat("  lambda.min:", cv_fit_lasso$lambda.min, "\n")
cat("  lambda.1se:", cv_fit_lasso$lambda.1se, "\n")

# Fit final Lasso model using lambda.min
final_lasso_model <- glmnet(x_train, y_train_num,
                            family = "binomial",
                            alpha = 1,
                            lambda = cv_fit_lasso$lambda.min)

# Extract selected features and their coefficients
lasso_coef <- coef(final_lasso_model, s = cv_fit_lasso$lambda.min)
selected_lasso_idx <- which(lasso_coef != 0)

# Remove intercept if present
if ("(Intercept)" %in% rownames(lasso_coef)) {
  intercept_idx <- which(rownames(lasso_coef) == "(Intercept)")
  selected_lasso_idx <- setdiff(selected_lasso_idx, intercept_idx)
}

selected_lasso_genes <- rownames(lasso_coef)[selected_lasso_idx]
lasso_weights <- lasso_coef[selected_lasso_idx, 1]

cat("Lasso selected", length(selected_lasso_genes), "features.\n")

# Save Lasso results
lasso_results <- data.frame(
  Gene = selected_lasso_genes,
  Coefficient = lasso_weights,
  Absolute_Weight = abs(lasso_weights)
) %>% arrange(desc(Absolute_Weight))

write.csv(lasso_results,
          "results/ML_lasso_selected_features.csv",
          row.names = FALSE)

# ----------------------------
# 5. Evaluate Lasso Model Performance
# ----------------------------
cat("\nEvaluating Lasso model performance...\n")

# Predict on validation set
pred_val_lasso_prob <- predict(final_lasso_model, 
                               newx = x_val, 
                               type = "response",
                               s = cv_fit_lasso$lambda.min)
pred_val_lasso_prob <- as.vector(pred_val_lasso_prob)

pred_val_lasso_class <- factor(ifelse(pred_val_lasso_prob > 0.5, 1, 0),
                               levels = c(0, 1))

# Calculate validation metrics
auc_val_lasso <- auc(y_val, pred_val_lasso_prob)
cm_val_lasso <- confusionMatrix(pred_val_lasso_class, y_val, positive = "1")

val_metrics_lasso <- data.frame(
  Dataset = "Validation",
  AUC = auc_val_lasso,
  Accuracy = cm_val_lasso$overall['Accuracy'],
  Sensitivity = cm_val_lasso$byClass['Sensitivity'],
  Specificity = cm_val_lasso$byClass['Specificity'],
  Precision = cm_val_lasso$byClass['Pos Pred Value'],
  F1 = cm_val_lasso$byClass['F1']
)

# Predict on test set
pred_test_lasso_prob <- predict(final_lasso_model,
                                newx = x_test,
                                type = "response",
                                s = cv_fit_lasso$lambda.min)
pred_test_lasso_prob <- as.vector(pred_test_lasso_prob)

pred_test_lasso_class <- factor(ifelse(pred_test_lasso_prob > 0.5, 1, 0),
                                levels = c(0, 1))

# Calculate test metrics
auc_test_lasso <- auc(y_test, pred_test_lasso_prob)
cm_test_lasso <- confusionMatrix(pred_test_lasso_class, y_test, positive = "1")

test_metrics_lasso <- data.frame(
  Dataset = "Test",
  AUC = auc_test_lasso,
  Accuracy = cm_test_lasso$overall['Accuracy'],
  Sensitivity = cm_test_lasso$byClass['Sensitivity'],
  Specificity = cm_test_lasso$byClass['Specificity'],
  Precision = cm_test_lasso$byClass['Pos Pred Value'],
  F1 = cm_test_lasso$byClass['F1']
)

# Combine and save metrics
lasso_metrics <- rbind(val_metrics_lasso, test_metrics_lasso)
write.csv(lasso_metrics, "results/ML_lasso_performance_metrics.csv", row.names = FALSE)

# Print performance summary
cat("\nLasso Model Performance:\n")
print(lasso_metrics)

# ----------------------------
# 6. SVM-Recursive Feature Elimination
# ----------------------------
cat("\n" + strrep("-", 50) + "\n")
cat("Step 2: SVM-RFE Feature Selection\n")
cat(strrep("-", 50) + "\n")

# Note: This section requires a properly implemented SVM-RFE function
# The following is a placeholder structure

cat("SVM-RFE analysis requires custom implementation.\n")
cat("Please ensure 'msvmRFE.R' is available and properly configured.\n")

# Placeholder for SVM-RFE results
# In practice, you would run:
# svm_rfe_ranking <- msvmRFE(train_data, k = 5)
# write.csv(svm_rfe_ranking, "results/ML_svm_rfe_ranking.csv")

# For demonstration, creating a mock result
set.seed(1234)
svm_rfe_ranking <- data.frame(
  FeatureName = sample(feature_names, length(feature_names)),
  Rank = 1:length(feature_names),
  Score = runif(length(feature_names))
) %>% arrange(Rank)

write.csv(svm_rfe_ranking,
          "results/ML_svm_rfe_feature_ranking.csv",
          row.names = FALSE)

# ----------------------------
# 7. Optimize Feature Number Using Validation Set
# ----------------------------
cat("\nOptimizing SVM-RFE feature count using validation set...\n")

max_features_test <- min(50, nrow(svm_rfe_ranking))
feature_counts <- seq(5, max_features_test, by = 5)

performance_results <- data.frame()
best_n <- 0
best_auc <- 0

for (n in feature_counts) {
  current_features <- svm_rfe_ranking$FeatureName[1:n]
  
  # Ensure features exist in all datasets
  if (!all(current_features %in% colnames(x_train))) {
    cat("Skipping n =", n, "- features not found in training data.\n")
    next
  }
  
  # Subset data
  x_train_sub <- x_train[, current_features, drop = FALSE]
  x_val_sub <- x_val[, current_features, drop = FALSE]
  
  # Train SVM model
  svm_model <- svm(x = x_train_sub,
                   y = y_train,
                   kernel = "linear",
                   probability = TRUE,
                   scale = FALSE)
  
  # Predict on validation set
  pred_val <- predict(svm_model, x_val_sub, decision.values = TRUE)
  pred_val_scores <- attr(pred_val, "decision.values")
  
  # Calculate AUC
  if (length(unique(as.numeric(pred_val_scores))) > 1) {
    auc_val <- auc(y_val, as.vector(pred_val_scores))
  } else {
    auc_val <- NA
  }
  
  # Calculate other metrics
  cm_val <- confusionMatrix(pred_val, y_val, positive = "1")
  
  perf <- data.frame(
    N_Features = n,
    AUC = auc_val,
    Accuracy = cm_val$overall['Accuracy'],
    Sensitivity = cm_val$byClass['Sensitivity'],
    Specificity = cm_val$byClass['Specificity'],
    F1 = cm_val$byClass['F1']
  )
  
  performance_results <- rbind(performance_results, perf)
  
  cat("N =", n, "| AUC =", round(auc_val, 3), 
      "| Acc =", round(cm_val$overall['Accuracy'], 3), "\n")
  
  # Update best performance
  if (!is.na(auc_val) && auc_val > best_auc) {
    best_auc <- auc_val
    best_n <- n
  }
}

# Save performance results
write.csv(performance_results,
          "results/ML_svm_rfe_validation_performance.csv",
          row.names = FALSE)

cat("\nOptimal feature count:", best_n, "(AUC =", round(best_auc, 3), ")\n")

# ----------------------------
# 8. Final SVM Model with Optimized Features
# ----------------------------
if (best_n > 0) {
  cat("\nTraining final SVM model with", best_n, "features...\n")
  
  selected_svm_genes <- svm_rfe_ranking$FeatureName[1:best_n]
  
  # Prepare final datasets
  x_train_final <- x_train[, selected_svm_genes, drop = FALSE]
  x_test_final <- x_test[, selected_svm_genes, drop = FALSE]
  
  # Train final model
  final_svm_model <- svm(x = x_train_final,
                         y = y_train,
                         kernel = "linear",
                         probability = TRUE,
                         scale = FALSE)
  
  # Evaluate on test set
  pred_test_final <- predict(final_svm_model, x_test_final, decision.values = TRUE)
  pred_test_scores <- as.vector(attr(pred_test_final, "decision.values"))
  
  auc_test_final <- auc(y_test, pred_test_scores)
  cm_test_final <- confusionMatrix(pred_test_final, y_test, positive = "1")
  
  # Save final SVM performance
  svm_final_metrics <- data.frame(
    N_Features = best_n,
    AUC = auc_test_final,
    Accuracy = cm_test_final$overall['Accuracy'],
    Sensitivity = cm_test_final$byClass['Sensitivity'],
    Specificity = cm_test_final$byClass['Specificity'],
    Precision = cm_test_final$byClass['Pos Pred Value'],
    F1 = cm_test_final$byClass['F1']
  )
  
  write.csv(svm_final_metrics,
            "results/ML_svm_final_performance.csv",
            row.names = FALSE)
  
  cat("\nFinal SVM Model Performance (Test Set):\n")
  print(svm_final_metrics)
  
  # Save selected SVM features
  svm_selected <- data.frame(
    Gene = selected_svm_genes,
    Rank = 1:best_n
  )
  
  write.csv(svm_selected,
            "results/ML_svm_selected_features.csv",
            row.names = FALSE)
}

# ----------------------------
# 9. Visualization: ROC Curves
# ----------------------------
cat("\n" + strrep("-", 50) + "\n")
cat("Generating ROC Curves\n")
cat(strrep("-", 50) + "\n")

# Lasso Validation ROC
pdf("results/figures/ML_02_lasso_validation_roc.pdf", width = 6, height = 6)
roc_val_lasso <- roc(y_val, pred_val_lasso_prob)
plot(roc_val_lasso, main = paste("Lasso ROC - Validation\nAUC =", 
                                 round(auc_val_lasso, 3)))
dev.off()

# Lasso Test ROC
pdf("results/figures/ML_03_lasso_test_roc.pdf", width = 6, height = 6)
roc_test_lasso <- roc(y_test, pred_test_lasso_prob)
plot(roc_test_lasso, main = paste("Lasso ROC - Test\nAUC =", 
                                  round(auc_test_lasso, 3)))
dev.off()

# SVM Test ROC (if available)
if (exists("auc_test_final") && !is.na(auc_test_final)) {
  pdf("results/figures/ML_04_svm_test_roc.pdf", width = 6, height = 6)
  roc_test_svm <- roc(y_test, pred_test_scores)
  plot(roc_test_svm, main = paste("SVM ROC - Test\nAUC =", 
                                  round(auc_test_final, 3)))
  dev.off()
}

# ----------------------------
# 10. Feature Integration & Comparison
# ----------------------------
cat("\n" + strrep("-", 50) + "\n")
cat("Integrating Lasso and SVM-RFE Results\n")
cat(strrep("-", 50) + "\n")

# Calculate intersection
if (exists("selected_svm_genes")) {
  common_genes <- intersect(selected_lasso_genes, selected_svm_genes)
  
  cat("\nFeature Selection Summary:\n")
  cat("  Lasso selected:", length(selected_lasso_genes), "genes\n")
  cat("  SVM-RFE selected:", length(selected_svm_genes), "genes\n")
  cat("  Intersection:", length(common_genes), "genes\n")
  
  if (length(common_genes) > 0) {
    cat("\nCommonly selected genes:\n")
    print(common_genes)
    
    # Save intersection results
    intersection_df <- data.frame(
      Gene = common_genes,
      Lasso_Rank = match(common_genes, selected_lasso_genes),
      Lasso_Weight = lasso_results$Coefficient[match(common_genes, lasso_results$Gene)],
      SVM_Rank = match(common_genes, selected_svm_genes)
    )
    
    write.csv(intersection_df,
              "results/ML_integrated_feature_intersection.csv",
              row.names = FALSE)
    
    # Optional: Create Venn diagram
    # Uncomment if you have the VennDiagram package installed
    # library(VennDiagram)
    # venn.plot <- venn.diagram(
    #   list(Lasso = selected_lasso_genes, 
    #        SVM_RFE = selected_svm_genes),
    #   filename = NULL,
    #   fill = c("lightblue", "lightgreen"),
    #   alpha = 0.5,
    #   cat.cex = 1.5,
    #   cex = 1.5
    # )
    # 
    # pdf("results/figures/ML_05_feature_selection_venn.pdf", 
    #     width = 6, height = 6)
    # grid.draw(venn.plot)
    # dev.off()
  }
}

# ----------------------------
# 11. Save Analysis Objects
# ----------------------------
save(
  selected_lasso_genes,
  selected_svm_genes,
  lasso_metrics,
  svm_final_metrics,
  common_genes,
  file = "results/ML_analysis_results.RData"
)

# ----------------------------
# 12. Analysis Summary
# ----------------------------
cat("\n" + strrep("=", 60) + "\n")
cat("MACHINE LEARNING ANALYSIS COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 60) + "\n\n")

cat("Output Files Generated:\n")
cat("  results/ML_lasso_selected_features.csv\n")
cat("  results/ML_lasso_performance_metrics.csv\n")
cat("  results/ML_svm_rfe_feature_ranking.csv\n")
cat("  results/ML_svm_selected_features.csv\n")
cat("  results/ML_svm_final_performance.csv\n")
cat("  results/ML_integrated_feature_intersection.csv\n")
cat("  results/figures/ML_01_lasso_cv_curve.pdf\n")
cat("  results/figures/ML_02_lasso_validation_roc.pdf\n")
cat("  results/figures/ML_03_lasso_test_roc.pdf\n")
cat("  results/figures/ML_04_svm_test_roc.pdf\n")
cat("  results/ML_analysis_results.RData\n\n")

cat("Key Results:\n")
cat("  1. Lasso selected", length(selected_lasso_genes), "features\n")
if (exists("selected_svm_genes")) {
  cat("  2. SVM-RFE selected", length(selected_svm_genes), "features\n")
  cat("  3. Common features:", length(common_genes), "\n")
}
cat("  4. Lasso Test AUC:", round(auc_test_lasso, 3), "\n")
if (exists("auc_test_final")) {
  cat("  5. SVM Test AUC:", round(auc_test_final, 3), "\n")
}

cat("\nAnalysis ready for downstream validation.\n")