# ==============================================================================
# Script: 03b_xgboost_feature_selection.R
# Project: TRIM29_UC_Mechanism
# Purpose: XGBoost-based Feature Selection with Validation Tuning
# Steps: 1. Prepare data for XGBoost
#        2. Hyperparameter tuning using validation set
#        3. Train final XGBoost model with optimal parameters
#        4. Extract feature importance and select top features
#        5. Evaluate model performance on test set
# Input:  - x_train, y_train_factor (training data)
#         - x_val, y_val_factor (validation data)
#         - x_test, y_test_factor (test data)
# Output: - Selected features from XGBoost
#         - Feature importance rankings
#         - Model performance metrics
#         - ROC curves
# ==============================================================================

# ----------------------------
# 0. Environment Setup
# ----------------------------
# Load required packages
library(xgboost)    # Extreme Gradient Boosting
library(caret)      # Confusion matrix and performance metrics
library(pROC)       # ROC curve analysis

# Set seed for reproducibility
set.seed(1234)

# ----------------------------
# 1. Data Preparation Check
# ----------------------------
cat(strrep("=", 60), "\n")
cat("XGBOOST FEATURE SELECTION WITH VALIDATION TUNING\n")
cat(strrep("=", 60), "\n\n")

# Check if required data exists
required_objects <- c("x_train", "y_train_factor", "x_val", "y_val_factor", "x_test", "y_test_factor")
missing_objects <- setdiff(required_objects, ls())

if (length(missing_objects) > 0) {
  stop("Missing required data objects: ", paste(missing_objects, collapse = ", "),
       "\nPlease run the main machine learning analysis script first.")
}

cat("Data availability confirmed:\n")
cat("  Training set:   ", nrow(x_train), "samples x", ncol(x_train), "features\n")
cat("  Validation set: ", nrow(x_val), "samples x", ncol(x_val), "features\n")
cat("  Test set:       ", nrow(x_test), "samples x", ncol(x_test), "features\n\n")

# Convert labels to numeric format for XGBoost
y_train_numeric <- as.numeric(y_train_factor) - 1  # 0/1 format
y_val_numeric <- as.numeric(y_val_factor) - 1
y_test_numeric <- as.numeric(y_test_factor) - 1

# Create XGBoost data matrices
dtrain <- xgb.DMatrix(data.matrix(x_train), label = y_train_numeric)
dval <- xgb.DMatrix(data.matrix(x_val), label = y_val_numeric)
dtest <- xgb.DMatrix(data.matrix(x_test), label = y_test_numeric)

# Define watchlist for monitoring training
watchlist <- list(train = dtrain, validation = dval)

# ----------------------------
# 2. Hyperparameter Tuning
# ----------------------------
cat("Step 1: Hyperparameter Tuning with Validation Set\n")
cat(strrep("-", 40), "\n")

# Define parameter grid for tuning
param_grid <- expand.grid(
  max_depth = c(3, 5, 6),          # Tree depth
  eta = c(0.01, 0.05, 0.1),        # Learning rate
  subsample = c(0.7, 0.8, 0.9),    # Row sampling
  colsample_bytree = c(0.7, 0.8, 0.9),  # Column sampling
  gamma = c(0, 0.1, 0.2)           # Minimum loss reduction
)

cat("Testing", nrow(param_grid), "parameter combinations...\n")

# Initialize tracking variables
best_auc <- 0
best_params <- NULL
best_nrounds <- 0
tuning_results <- data.frame()

# Grid search for optimal parameters
for (i in 1:nrow(param_grid)) {
  
  current_params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = param_grid$max_depth[i],
    eta = param_grid$eta[i],
    subsample = param_grid$subsample[i],
    colsample_bytree = param_grid$colsample_bytree[i],
    gamma = param_grid$gamma[i],
    scale_pos_weight = sum(y_train_numeric == 0) / sum(y_train_numeric == 1)  # Handle class imbalance
  )
  
  # Train model with early stopping
  xgb_model <- xgb.train(
    params = current_params,
    data = dtrain,
    nrounds = 200,
    watchlist = watchlist,
    verbose = 0,
    early_stopping_rounds = 20,
    maximize = TRUE
  )
  
  # Extract validation AUC
  val_auc <- max(xgb_model$evaluation_log$validation_auc)
  best_iter <- xgb_model$best_iteration
  
  # Store results
  tuning_results <- rbind(tuning_results, data.frame(
    Combination = i,
    max_depth = param_grid$max_depth[i],
    eta = param_grid$eta[i],
    subsample = param_grid$subsample[i],
    colsample_bytree = param_grid$colsample_bytree[i],
    gamma = param_grid$gamma[i],
    AUC = val_auc,
    Best_Iteration = best_iter
  ))
  
  # Update best parameters if improvement found
  if (val_auc > best_auc) {
    best_auc <- val_auc
    best_params <- current_params
    best_nrounds <- best_iter
  }
  
  # Progress update
  if (i %% 10 == 0) {
    cat("  Completed", i, "/", nrow(param_grid), "combinations\n")
  }
}

# Save tuning results
write.csv(tuning_results, 
          "results/ML_xgboost_hyperparameter_tuning.csv", 
          row.names = FALSE)

cat("\nHyperparameter tuning completed.\n")
cat("Best validation AUC:", round(best_auc, 4), "\n")
cat("Best iteration count:", best_nrounds, "\n")
cat("Optimal parameters:\n")
print(as.data.frame(best_params[1:6]))

# ----------------------------
# 3. Train Final XGBoost Model
# ----------------------------
cat("\n", strrep("-", 40), "\n", sep = "")
cat("Step 2: Training Final XGBoost Model\n")
cat(strrep("-", 40), "\n")

# Train final model with optimal parameters
final_xgb_model <- xgb.train(
  params = best_params,
  data = dtrain,
  nrounds = best_nrounds,
  watchlist = watchlist,
  verbose = 0
)

cat("Final XGBoost model trained successfully.\n")

# ----------------------------
# 4. Feature Importance Analysis
# ----------------------------
cat("\n", strrep("-", 40), "\n", sep = "")
cat("Step 3: Feature Importance Extraction\n")
cat(strrep("-", 40), "\n")

# Extract feature importance
feature_importance <- xgb.importance(
  model = final_xgb_model,
  feature_names = colnames(x_train)
)

# Display top features
cat("\nTop 20 Features by Importance:\n")
print(head(feature_importance, 20))

# Define selection threshold (top N features or importance score)
n_top_features <- 15  # Adjust based on your needs
importance_threshold <- 0.01  # Minimum relative importance

# Select features using multiple criteria
selected_by_rank <- head(feature_importance$Feature, n_top_features)
selected_by_threshold <- feature_importance$Feature[feature_importance$Gain >= importance_threshold]

# Use union of both selection methods
selected_features <- unique(c(selected_by_rank, selected_by_threshold))

cat("\nFeature Selection Summary:\n")
cat("  Selected by top", n_top_features, "ranking:", length(selected_by_rank), "features\n")
cat("  Selected by threshold (Gain >=", importance_threshold, "):", 
    length(selected_by_threshold), "features\n")
cat("  Total unique features selected:", length(selected_features), "\n")

# Save feature importance results
feature_importance$Rank <- 1:nrow(feature_importance)
write.csv(feature_importance, 
          "results/ML_xgboost_feature_importance.csv", 
          row.names = FALSE)

# Save selected features
selected_features_df <- data.frame(
  Feature = selected_features,
  Rank = match(selected_features, feature_importance$Feature),
  Gain = feature_importance$Gain[match(selected_features, feature_importance$Feature)]
)
write.csv(selected_features_df, 
          "results/ML_xgboost_selected_features.csv", 
          row.names = FALSE)

# ----------------------------
# 5. Model Evaluation on Test Set
# ----------------------------
cat("\n", strrep("-", 40), "\n", sep = "")
cat("Step 4: Test Set Performance Evaluation\n")
cat(strrep("-", 40), "\n")

# Make predictions on test set
test_predictions <- predict(final_xgb_model, dtest)
test_pred_classes <- ifelse(test_predictions > 0.5, 1, 0)

# Calculate performance metrics
conf_matrix <- confusionMatrix(
  factor(test_pred_classes, levels = c(0, 1)),
  factor(y_test_numeric, levels = c(0, 1)),
  positive = "1"
)

# ROC curve analysis
roc_curve <- roc(y_test_numeric, test_predictions)
auc_value <- auc(roc_curve)

# Compile performance metrics
performance_metrics <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity", 
             "Precision", "F1-Score", "Balanced Accuracy"),
  Value = c(
    round(auc_value, 4),
    round(conf_matrix$overall['Accuracy'], 4),
    round(conf_matrix$byClass['Sensitivity'], 4),
    round(conf_matrix$byClass['Specificity'], 4),
    round(conf_matrix$byClass['Pos Pred Value'], 4),
    round(conf_matrix$byClass['F1'], 4),
    round(conf_matrix$byClass['Balanced Accuracy'], 4)
  )
)

cat("\nTest Set Performance:\n")
print(performance_metrics)

# Save performance metrics
write.csv(performance_metrics, 
          "results/ML_xgboost_test_performance.csv", 
          row.names = FALSE)

# ----------------------------
# 6. Visualization
# ----------------------------
cat("\n", strrep("-", 40), "\n", sep = "")
cat("Step 5: Generating Visualizations\n")
cat(strrep("-", 40), "\n")

# 6.1 ROC Curve
pdf("results/figures/ML_05_xgboost_roc_curve.pdf", width = 6, height = 6)
plot(roc_curve, main = paste("XGBoost ROC Curve\nAUC =", round(auc_value, 3)),
     col = "blue", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()
cat("✓ ROC curve saved: results/figures/ML_05_xgboost_roc_curve.pdf\n")

# 6.2 Feature Importance Plot (Top 20)
if (nrow(feature_importance) >= 20) {
  top_features <- head(feature_importance, 20)
  
  pdf("results/figures/ML_06_xgboost_feature_importance.pdf", 
      width = 8, height = 6)
  
  par(mar = c(5, 12, 4, 2) + 0.1)  # Adjust margins for long feature names
  
  barplot(rev(top_features$Gain),
          horiz = TRUE,
          names.arg = rev(top_features$Feature),
          las = 1,
          col = "steelblue",
          main = "Top 20 Features by Importance (XGBoost)",
          xlab = "Gain",
          cex.names = 0.7)
  
  dev.off()
  cat("✓ Feature importance plot saved: results/figures/ML_06_xgboost_feature_importance.pdf\n")
}

# 6.3 Prediction Distribution
pdf("results/figures/ML_07_xgboost_prediction_distribution.pdf", 
    width = 8, height = 6)

par(mfrow = c(1, 2))

# Distribution by true class
hist(test_predictions[y_test_numeric == 0], 
     col = rgb(1, 0, 0, 0.5), 
     breaks = 20, 
     main = "Prediction Distribution\nby True Class",
     xlab = "Predicted Probability",
     xlim = c(0, 1),
     ylim = c(0, max(table(cut(test_predictions, breaks = 20)))))

hist(test_predictions[y_test_numeric == 1], 
     col = rgb(0, 0, 1, 0.5), 
     breaks = 20, 
     add = TRUE)

legend("topright", 
       legend = c("Class 0", "Class 1"), 
       fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))

# Calibration plot
calibration_data <- data.frame(
  Predicted = test_predictions,
  Actual = y_test_numeric
)

calibration_data$Bin <- cut(calibration_data$Predicted, 
                            breaks = seq(0, 1, by = 0.1),
                            include.lowest = TRUE)

bin_stats <- aggregate(Actual ~ Bin, data = calibration_data, FUN = mean)

plot(bin_stats$Actual ~ seq(0.05, 0.95, by = 0.1),
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Mean Predicted Probability",
     ylab = "Observed Proportion",
     main = "Calibration Plot",
     pch = 19, col = "blue")
abline(a = 0, b = 1, lty = 2, col = "gray")

par(mfrow = c(1, 1))
dev.off()
cat("✓ Prediction distribution plots saved: results/figures/ML_07_xgboost_prediction_distribution.pdf\n")

# ----------------------------
# 7. Save Model and Results
# ----------------------------
cat("\n", strrep("-", 40), "\n", sep = "")
cat("Step 6: Saving Model and Analysis Objects\n")
cat(strrep("-", 40), "\n")

# Save the trained model
xgb.save(final_xgb_model, "results/ML_xgboost_model.model")
cat("✓ XGBoost model saved: results/ML_xgboost_model.model\n")

# Save complete analysis objects
save(
  final_xgb_model,
  feature_importance,
  selected_features,
  performance_metrics,
  roc_curve,
  file = "results/ML_xgboost_analysis_results.RData"
)
cat("✓ Complete analysis objects saved: results/ML_xgboost_analysis_results.RData\n")

# ----------------------------
# 8. Analysis Summary
# ----------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat("XGBOOST FEATURE SELECTION - ANALYSIS SUMMARY\n")
cat(strrep("=", 60), "\n\n")

cat("Key Results:\n")
cat("  1. Optimal parameters identified through grid search\n")
cat("  2. Validation AUC:", round(best_auc, 4), "\n")
cat("  3. Test set AUC:", round(auc_value, 4), "\n")
cat("  4. Features selected:", length(selected_features), "\n\n")

cat("Output Files Generated:\n")
cat("  results/ML_xgboost_hyperparameter_tuning.csv\n")
cat("  results/ML_xgboost_feature_importance.csv\n")
cat("  results/ML_xgboost_selected_features.csv\n")
cat("  results/ML_xgboost_test_performance.csv\n")
cat("  results/ML_xgboost_model.model\n")
cat("  results/ML_xgboost_analysis_results.RData\n")
cat("  results/figures/ML_05_xgboost_roc_curve.pdf\n")
cat("  results/figures/ML_06_xgboost_feature_importance.pdf\n")
cat("  results/figures/ML_07_xgboost_prediction_distribution.pdf\n\n")

cat("Selected Features (Top", n_top_features, "):\n")
for (i in 1:min(n_top_features, length(selected_features))) {
  cat("  ", i, ". ", selected_features[i], "\n", sep = "")
}

cat("\nAnalysis complete. These features can be integrated with Lasso and SVM-RFE results.\n")
cat(strrep("=", 60), "\n")