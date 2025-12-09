# ==============================================================================
# Script: 02_mendelian_randomization.R
# Project: TRIM29_UC_Mechanism
# Purpose: Mendelian Randomization Analysis for UC and Lipid Metabolism
# Steps: 1. Load and prepare exposure/outcome IDs
#        2. Perform MR analysis for all exposure-outcome pairs
#        3. Sensitivity analyses (heterogeneity, pleiotropy, MR-PRESSO)
#        4. Export results to structured Excel files
# Input:  - Exposure IDs (UC-related GWAS)
#         - Outcome IDs (Lipid metabolism GWAS)
# Output: - MR results with odds ratios
#         - Sensitivity analysis results
#         - Harmonised data for each analysis
# ==============================================================================

# ----------------------------
# 0. Environment Setup
# ----------------------------
# Load required packages
library(TwoSampleMR)  # Core MR analysis functions
library(MRPRESSO)     # MR-PRESSO for outlier detection
library(ieugwasr)     # Interface to GWAS database
library(openxlsx)     # Excel file export

# Set seed for reproducibility
set.seed(123456)

# ----------------------------
# 1. Define Exposure and Outcome IDs
# ----------------------------
# UC-related exposure IDs from various consortia
exposure_ids <- c(
  "ieu-a-970", "ukb-b-19386", "ukb-b-7584", "ukb-a-104",
  "ebi-a-GCST90018933", "ebi-a-GCST90038684", "ebi-a-GCST90018713",
  "ukb-a-553", "ukb-d-ULCERNAS", "finn-b-ULCERNAS",
  "finn-b-K11_UC_STRICT_PSC", "finn-b-K11_UC_STRICT",
  "finn-b-K11_UC_STRICT2", "finn-b-K11_UC_NOCD",
  "finn-b-K11_ULCER", "finn-b-ULCEROTH", "ieu-a-972",
  "ebi-a-GCST90020072", "ieu-a-971", "ieu-a-973",
  "ebi-a-GCST000964", "ieu-a-968", "ieu-a-32",
  "ebi-a-GCST003045", "ebi-a-GCST004133", "ukb-a-553"
)

# Lipid metabolism outcome IDs
outcome_ids <- c(
  'ebi-a-GCST90104007', 'ebi-a-GCST90104006',
  'ebi-a-GCST90104005', 'ebi-a-GCST90104004',
  'ebi-a-GCST90104003', "ebi-a-GCST90090994"
)

# Report analysis plan
cat("Mendelian Randomization Analysis Plan:\n")
cat("  Exposures:", length(exposure_ids), "UC-related GWAS\n")
cat("  Outcomes: ", length(outcome_ids), "lipid metabolism GWAS\n")
cat("  Total pairs to analyze:", length(exposure_ids) * length(outcome_ids), "\n\n")

# Methods of interest for MR analysis
methods_of_interest <- c(
  "Inverse variance weighted (fixed effects)",
  "Inverse variance weighted (multiplicative random effects)",
  "Inverse variance weighted",
  "MR Egger",
  "Weighted median"
)

# ----------------------------
# 2. Initialize Results Storage
# ----------------------------
all_results <- list(
  harmonised_data = list(),
  mr_results = list(),
  or_results = list(),
  heterogeneity_results = list(),
  pleiotropy_results = list(),
  presso_results = list(),
  exposure_nsnp = list()  # Number of SNPs per exposure
)

# ----------------------------
# 3. Main Analysis Loop
# ----------------------------
cat(strrep("=", 60), "\n")
cat("STARTING MENDELIAN RANDOMIZATION ANALYSES\n")
cat(strrep("=", 60), "\n\n")

# Outer loop: Process each outcome
for (outcome_id in outcome_ids) {
  cat("\n", strrep("-", 50), "\n", sep = "")
  cat("PROCESSING OUTCOME:", outcome_id, "\n")
  cat(strrep("-", 50), "\n")
  
  # Initialize storage for current outcome
  all_results$harmonised_data[[outcome_id]] <- list()
  all_results$mr_results[[outcome_id]] <- list()
  all_results$or_results[[outcome_id]] <- list()
  all_results$heterogeneity_results[[outcome_id]] <- list()
  all_results$pleiotropy_results[[outcome_id]] <- list()
  all_results$presso_results[[outcome_id]] <- list()
  all_results$exposure_nsnp[[outcome_id]] <- list()
  
  # Inner loop: Process each exposure against current outcome
  for (exposure_id in exposure_ids) {
    cat("\n  Analysis: ", exposure_id, " → ", outcome_id, "\n", sep = "")
    
    # 3.1 Extract instrument variables (exposure data)
    exposure_data <- tryCatch({
      extract_instruments(
        outcomes = exposure_id,
        p1 = 5e-08,      # GWAS significance threshold
        clump = TRUE,    # LD clumping
        r2 = 0.001,      # LD r² threshold
        kb = 10000       # Clumping window
      )
    }, error = function(e) {
      cat("    WARNING: Failed to extract instruments for", exposure_id, "\n")
      cat("    Error:", e$message, "\n")
      return(NULL)
    })
    
    # Skip if no instruments found
    if (is.null(exposure_data) || nrow(exposure_data) == 0) {
      cat("    SKIPPING: No valid instruments found\n")
      next
    }
    
    cat("    Found", nrow(exposure_data), "instrumental variables\n")
    
    # 3.2 Extract outcome data for these SNPs
    outcome_data <- tryCatch({
      extract_outcome_data(
        snps = exposure_data$SNP,
        outcomes = outcome_id,
        proxies = FALSE,
        maf_threshold = 0.01
      )
    }, error = function(e) {
      cat("    WARNING: Failed to extract outcome data\n")
      cat("    Error:", e$message, "\n")
      return(NULL)
    })
    
    # Skip if no overlapping SNPs
    if (is.null(outcome_data) || nrow(outcome_data) == 0) {
      cat("    SKIPPING: No overlapping SNPs found\n")
      next
    }
    
    # 3.3 Harmonize exposure and outcome data
    harmonised_data <- tryCatch({
      harmonise_data(
        exposure_dat = exposure_data,
        outcome_dat = outcome_data
      )
    }, error = function(e) {
      cat("    WARNING: Failed to harmonize data\n")
      cat("    Error:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(harmonised_data) || nrow(harmonised_data) == 0) {
      cat("    SKIPPING: Harmonization failed\n")
      next
    }
    
    cat("    Harmonised", nrow(harmonised_data), "SNP pairs\n")
    
    # Store harmonised data
    all_results$harmonised_data[[outcome_id]][[exposure_id]] <- harmonised_data
    all_results$exposure_nsnp[[outcome_id]][[exposure_id]] <- nrow(exposure_data)
    
    # 3.4 Perform MR analysis
    mr_result <- tryCatch({
      mr(
        harmonised_data,
        method_list = c(
          'mr_ivw_fe',
          'mr_ivw_mre',
          'mr_ivw',
          'mr_egger_regression',
          'mr_weighted_median'
        )
      )
    }, error = function(e) {
      cat("    WARNING: MR analysis failed\n")
      cat("    Error:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(mr_result) || nrow(mr_result) == 0) {
      cat("    SKIPPING: MR analysis produced no results\n")
      all_results$mr_results[[outcome_id]][[exposure_id]] <- NULL
      next
    }
    
    # 3.5 Calculate odds ratios
    mr_result_or <- tryCatch({
      generate_odds_ratios(mr_result)
    }, error = function(e) {
      cat("    WARNING: OR calculation failed, using raw MR results\n")
      cat("    Error:", e$message, "\n")
      return(mr_result)
    })
    
    # Store MR results
    all_results$mr_results[[outcome_id]][[exposure_id]] <- mr_result_or
    all_results$or_results[[outcome_id]][[exposure_id]] <- mr_result_or
    
    # 3.6 Sensitivity analyses
    cat("    Running sensitivity analyses...\n")
    
    # Heterogeneity test
    heterogeneity_result <- tryCatch({
      mr_heterogeneity(
        harmonised_data,
        method_list = c("mr_egger_regression", "mr_ivw")
      )
    }, error = function(e) {
      cat("    WARNING: Heterogeneity test failed\n")
      return(NULL)
    })
    all_results$heterogeneity_results[[outcome_id]][[exposure_id]] <- heterogeneity_result
    
    # Pleiotropy test
    pleiotropy_result <- tryCatch({
      mr_pleiotropy_test(harmonised_data)
    }, error = function(e) {
      cat("    WARNING: Pleiotropy test failed\n")
      return(NULL)
    })
    all_results$pleiotropy_results[[outcome_id]][[exposure_id]] <- pleiotropy_result
    
    # MR-PRESSO outlier detection
    presso_result <- tryCatch({
      mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = harmonised_data,
        NbDistribution = 1000,
        SignifThreshold = 0.05,
        seed = 123456
      )
    }, error = function(e) {
      cat("    WARNING: MR-PRESSO analysis failed\n")
      return(NULL)
    })
    all_results$presso_results[[outcome_id]][[exposure_id]] <- presso_result
    
    cat("    ✓ Analysis completed successfully\n")
    
  }  # End inner loop (exposures)
  
  # ----------------------------
  # 4. Compile Results for Current Outcome
  # ----------------------------
  cat("\n  Compiling results for outcome:", outcome_id, "\n")
  
  # 4.1 Compile main MR results
  compiled_mr_results <- NULL
  
  if (length(all_results$mr_results[[outcome_id]]) > 0) {
    result_list <- list()
    
    for (exp_id in names(all_results$mr_results[[outcome_id]])) {
      res <- all_results$mr_results[[outcome_id]][[exp_id]]
      
      if (!is.null(res) && nrow(res) > 0) {
        # Filter for methods of interest
        selected_res <- res[res$method %in% methods_of_interest, , drop = FALSE]
        
        if (nrow(selected_res) > 0) {
          selected_res$exposure <- exp_id
          selected_res$outcome <- outcome_id
          result_list[[exp_id]] <- selected_res
        }
      }
    }
    
    if (length(result_list) > 0) {
      compiled_mr_results <- do.call(rbind, result_list)
      rownames(compiled_mr_results) <- NULL
      cat("    Compiled", nrow(compiled_mr_results), "MR results\n")
    }
  }
  
  # 4.2 Export results if available
  if (!is.null(compiled_mr_results) && nrow(compiled_mr_results) > 0) {
    
    # Create safe filename
    safe_outcome_name <- gsub("[^a-zA-Z0-9._-]", "_", outcome_id)
    
    # 4.2.1 Export all MR results
    all_results_file <- paste0("results/MR_Results_All_vs_", safe_outcome_name, ".xlsx")
    
    cat("    Exporting all MR results to:", all_results_file, "\n")
    tryCatch({
      writexl::write_xlsx(compiled_mr_results, all_results_file)
      cat("    ✓ All results exported successfully\n")
    }, error = function(e) {
      cat("    ERROR exporting all results:", e$message, "\n")
    })
    
    # 4.2.2 Export significant results (p < 0.05)
    significant_rows <- which(compiled_mr_results$pval < 0.05 & 
                                !is.na(compiled_mr_results$pval))
    
    if (length(significant_rows) > 0) {
      significant_results <- compiled_mr_results[significant_rows, , drop = FALSE]
      rownames(significant_results) <- NULL
      
      sig_results_file <- paste0("results/MR_Results_Significant_vs_", 
                                 safe_outcome_name, ".xlsx")
      
      cat("    Found", nrow(significant_results), "significant results (p < 0.05)\n")
      cat("    Exporting to:", sig_results_file, "\n")
      
      tryCatch({
        writexl::write_xlsx(significant_results, sig_results_file)
        cat("    ✓ Significant results exported successfully\n")
      }, error = function(e) {
        cat("    ERROR exporting significant results:", e$message, "\n")
      })
    } else {
      cat("    No significant results found (p < 0.05)\n")
    }
  }
  
  # ----------------------------
  # 5. Export Sensitivity Analysis Results
  # ----------------------------
  cat("\n  Exporting sensitivity analysis results...\n")
  
  # 5.1 Compile sensitivity results
  sensitivity_wb <- createWorkbook()
  sheets_added <- 0
  
  # Heterogeneity results
  if (length(all_results$heterogeneity_results[[outcome_id]]) > 0) {
    het_list <- list()
    
    for (exp_id in names(all_results$heterogeneity_results[[outcome_id]])) {
      het_res <- all_results$heterogeneity_results[[outcome_id]][[exp_id]]
      
      if (!is.null(het_res) && nrow(het_res) > 0) {
        het_res$exposure <- exp_id
        het_res$outcome <- outcome_id
        het_list[[exp_id]] <- het_res
      }
    }
    
    if (length(het_list) > 0) {
      het_combined <- do.call(rbind, het_list)
      addWorksheet(sensitivity_wb, "Heterogeneity")
      writeData(sensitivity_wb, "Heterogeneity", het_combined)
      sheets_added <- sheets_added + 1
      cat("    ✓ Heterogeneity results added\n")
    }
  }
  
  # Pleiotropy results
  if (length(all_results$pleiotropy_results[[outcome_id]]) > 0) {
    pleio_list <- list()
    
    for (exp_id in names(all_results$pleiotropy_results[[outcome_id]])) {
      pleio_res <- all_results$pleiotropy_results[[outcome_id]][[exp_id]]
      
      if (!is.null(pleio_res) && nrow(pleio_res) > 0) {
        pleio_res$exposure <- exp_id
        pleio_res$outcome <- outcome_id
        pleio_list[[exp_id]] <- pleio_res
      }
    }
    
    if (length(pleio_list) > 0) {
      pleio_combined <- do.call(rbind, pleio_list)
      addWorksheet(sensitivity_wb, "Pleiotropy")
      writeData(sensitivity_wb, "Pleiotropy", pleio_combined)
      sheets_added <- sheets_added + 1
      cat("    ✓ Pleiotropy results added\n")
    }
  }
  
  # MR-PRESSO results
  if (length(all_results$presso_results[[outcome_id]]) > 0) {
    presso_global_list <- list()
    presso_outlier_list <- list()
    
    for (exp_id in names(all_results$presso_results[[outcome_id]])) {
      presso_res <- all_results$presso_results[[outcome_id]][[exp_id]]
      
      if (!is.null(presso_res)) {
        # Global test
        if (!is.null(presso_res$`Global Test`) && 
            nrow(presso_res$`Global Test`) > 0) {
          global_df <- presso_res$`Global Test`
          global_df$exposure <- exp_id
          global_df$outcome <- outcome_id
          presso_global_list[[exp_id]] <- global_df
        }
        
        # Outlier test
        if (!is.null(presso_res$`Outlier Test`) && 
            nrow(presso_res$`Outlier Test`) > 0) {
          outlier_df <- presso_res$`Outlier Test`
          outlier_df$exposure <- exp_id
          outlier_df$outcome <- outcome_id
          presso_outlier_list[[exp_id]] <- outlier_df
        }
      }
    }
    
    # Add global test sheet
    if (length(presso_global_list) > 0) {
      presso_global_combined <- do.call(rbind, presso_global_list)
      addWorksheet(sensitivity_wb, "MR_PRESSO_Global")
      writeData(sensitivity_wb, "MR_PRESSO_Global", presso_global_combined)
      sheets_added <- sheets_added + 1
      cat("    ✓ MR-PRESSO Global test results added\n")
    }
    
    # Add outlier test sheet
    if (length(presso_outlier_list) > 0) {
      presso_outlier_combined <- do.call(rbind, presso_outlier_list)
      addWorksheet(sensitivity_wb, "MR_PRESSO_Outlier")
      writeData(sensitivity_wb, "MR_PRESSO_Outlier", presso_outlier_combined)
      sheets_added <- sheets_added + 1
      cat("    ✓ MR-PRESSO Outlier test results added\n")
    }
  }
  
  # 5.2 Save sensitivity workbook if sheets were added
  if (sheets_added > 0) {
    sensitivity_file <- paste0("results/MR_Sensitivity_", safe_outcome_name, ".xlsx")
    
    tryCatch({
      saveWorkbook(sensitivity_wb, sensitivity_file, overwrite = TRUE)
      cat("    ✓ Sensitivity results saved to:", sensitivity_file, "\n")
    }, error = function(e) {
      cat("    ERROR saving sensitivity results:", e$message, "\n")
    })
  } else {
    cat("    No sensitivity results to export\n")
  }
  
  cat("\n  ✓ Outcome", outcome_id, "processing complete\n")
  
}  # End outer loop (outcomes)

# ----------------------------
# 6. Save Complete Analysis Objects
# ----------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat("SAVING ANALYSIS OBJECTS\n")
cat(strrep("=", 60), "\n")

# Save key results for potential reuse
save(
  all_results,
  exposure_ids,
  outcome_ids,
  file = "results/MR_analysis_complete.RData"
)

cat("✓ Analysis objects saved to: results/MR_analysis_complete.RData\n")

# ----------------------------
# 7. Generate Analysis Summary
# ----------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat("ANALYSIS SUMMARY\n")
cat(strrep("=", 60), "\n")

# Count successful analyses
successful_analyses <- 0
for (outcome in names(all_results$mr_results)) {
  successful_analyses <- successful_analyses + 
    length(all_results$mr_results[[outcome]])
}

cat("Total exposure-outcome pairs analyzed:", successful_analyses, "\n")
cat("Exposure IDs processed:", length(exposure_ids), "\n")
cat("Outcome IDs processed:", length(outcome_ids), "\n\n")

cat("Output Files Generated:\n")
cat("  Main MR results: results/MR_Results_All_vs_[OUTCOME].xlsx\n")
cat("  Significant results: results/MR_Results_Significant_vs_[OUTCOME].xlsx\n")
cat("  Sensitivity analyses: results/MR_Sensitivity_[OUTCOME].xlsx\n")
cat("  Complete R objects: results/MR_analysis_complete.RData\n\n")

cat("MENDELIAN RANDOMIZATION ANALYSIS COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 60), "\n")