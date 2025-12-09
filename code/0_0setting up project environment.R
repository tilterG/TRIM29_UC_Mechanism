# ==============================================================================
# 0. PROJECT ENVIRONMENT SETUP
# ==============================================================================
cat("[0] SETTING UP PROJECT ENVIRONMENT\n")
cat(rep("-", 40), "\n")

# 0.1 Create project directory structure
cat("  Creating directory structure...\n")
required_dirs <- c(
  "results",                 # Main results directory
  "results/figures",         # For all plot outputs (PDF, PNG)
  "results/tables",          # For all data table outputs (CSV, XLSX)
  "results/intermediate",    # For intermediate R objects
  "data/raw",                # For raw input data (if any)
  "data/processed",          # For processed/cleaned data
  "code",                    # For analysis scripts
  "code/utils",              # For utility functions
  "logs"                     # For log files
)

# Create all required directories
dirs_created <- 0
for (dir_path in required_dirs) {
  if (!dir.exists(dir_path)) {
    success <- tryCatch({
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      TRUE
    }, error = function(e) FALSE)
    
    if (success) {
      cat("    ✓ Created: ", dir_path, "\n", sep = "")
      dirs_created <- dirs_created + 1
    } else {
      cat("    ✗ Failed to create: ", dir_path, "\n", sep = "")
    }
  } else {
    cat("    ✓ Exists: ", dir_path, "\n", sep = "")
  }
}

cat("  Directory setup complete: ", dirs_created, " directories checked/created\n\n")

# 0.2 Verify required analysis scripts exist
cat("  Verifying analysis scripts...\n")
required_scripts <- c(
  "code/01_WGCNA_analysis.R",
  "code/02_mendelian_randomization.R",
  "code/03_machine_learning_feature_selection.R",
  "code/04_xgboost_feature_selection.R"
)

script_status <- sapply(required_scripts, function(script_path) {
  if (file.exists(script_path)) {
    file_size <- file.info(script_path)$size
    list(exists = TRUE, size = file_size)
  } else {
    list(exists = FALSE, size = 0)
  }
})

missing_scripts <- names(script_status)[!sapply(script_status, function(x) x$exists)]

if (length(missing_scripts) > 0) {
  cat("  ✗ ERROR: Missing required script(s):\n")
  for (script in missing_scripts) {
    cat("      - ", script, "\n", sep = "")
  }
  cat("\n  Please ensure all analysis scripts are in the 'code/' directory.\n")
  stop("Critical scripts missing. Execution halted.")
} else {
  cat("  ✓ All required scripts found:\n")
  for (script_name in names(script_status)) {
    size_kb <- round(script_status[[script_name]]$size / 1024, 1)
    cat("      - ", basename(script_name), " (", size_kb, " KB)\n", sep = "")
  }
}
cat("  Script verification complete.\n\n")

# 0.3 Check for critical input data
cat("  Checking for input data...\n")
critical_data_files <- c(
  "data/processed/ml_input_data.csv",  # For machine learning
  "data/WGCNA.fpkm",                   # For WGCNA analysis (if saved as R object)
  "data/allTraits"                     # For WGCNA traits
)

data_status <- sapply(critical_data_files, file.exists)

if (!all(data_status)) {
  missing_data <- names(data_status)[!data_status]
  cat("  ⚠ WARNING: Some input data files not found:\n")
  for (data_file in missing_data) {
    cat("      - ", data_file, "\n", sep = "")
  }
  cat("\n  Please ensure required data files are in the correct locations.\n")
  cat("  The scripts may fail if data is not available.\n")
} else {
  cat("  ✓ All critical data files found.\n")
}
cat("  Data check complete.\n\n")

# 0.4 Initialize log file
log_file <- paste0("logs/analysis_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
if (!dir.exists("logs")) dir.create("logs", showWarnings = FALSE)

log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- paste("[", timestamp, "] ", msg, sep = "")
  cat(log_line, "\n", file = log_file, append = TRUE)
}

# Create initial log entry
log_message("TRIM29_UC_Mechanism Analysis Pipeline Started")
log_message(paste("Project root:", getwd()))
log_message(paste("R version:", R.version.string))
log_message(paste("Platform:", R.version$platform))

cat("  ✓ Log file initialized: ", basename(log_file), "\n\n")

# 0.5 Record system information
cat("  Recording system information...\n")
sys_info <- list(
  start_time = Sys.time(),
  r_version = R.version.string,
  platform = R.version$platform,
  os = Sys.info()["sysname"],
  user = Sys.info()["user"],
  working_dir = getwd(),
  total_scripts = length(required_scripts),
  scripts_found = length(required_scripts) - length(missing_scripts)
)

cat("    Start time:   ", format(sys_info$start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("    R version:    ", sys_info$r_version, "\n")
cat("    Working dir:  ", sys_info$working_dir, "\n")
cat("    Scripts:      ", sys_info$scripts_found, "/", sys_info$total_scripts, "\n")

# Log system info
log_message(paste("System info recorded - R", sys_info$r_version))
log_message(paste("OS:", sys_info$os))
log_message(paste("User:", sys_info$user))

# 0.6 Check for renv environment
cat("\n  Checking R environment...\n")
if (file.exists("renv.lock")) {
  if (!requireNamespace("renv", quietly = TRUE)) {
    cat("    ⚠ 'renv' package not installed. Installing...\n")
    install.packages("renv")
  }
  
  # Check if renv is activated
  if (is.null(renv::project())) {
    cat("    ℹ 'renv' environment found but not activated.\n")
    cat("    Consider running 'renv::restore()' to ensure package consistency.\n")
    log_message("renv.lock found but environment not activated")
  } else {
    cat("    ✓ 'renv' environment is active.\n")
    log_message("renv environment is active")
  }
} else {
  cat("    ⚠ No 'renv.lock' file found. Package versions are not locked.\n")
  log_message("No renv.lock file found - package versions not controlled")
}

# 0.7 Memory check (optional, for large datasets)
cat("\n  System resources:\n")
cat("    Memory available: ", round(system("free -h | awk 'NR==2{print $7}'", intern = TRUE), 1), " MB\n", sep = "")
cat("    CPU cores: ", parallel::detectCores(), "\n", sep = "")

log_message(paste("CPU cores:", parallel::detectCores()))

# 0.8 Summary
cat("\n", rep("-", 40), "\n", sep = "")
cat("ENVIRONMENT SETUP COMPLETE\n")
cat(rep("-", 40), "\n\n")

cat("Project is ready for execution.\n")
cat("Next: Starting analysis pipeline...\n\n")

# Log completion of setup
log_message("Environment setup completed successfully")
log_message(rep("-", 50))

# Pause briefly (optional, for user to read output)
if (interactive()) {
  Sys.sleep(1)
}