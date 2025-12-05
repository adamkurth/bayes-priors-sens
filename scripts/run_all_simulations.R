###--------------------------------------------------------------------------###
### Master Script: Run All Simulations and Analysis
###
### This script orchestrates the full analysis pipeline:
### 1. Generate data batches (if not already present)
### 2. Fit Stan models across all scenarios
### 3. Generate results for the write-up
###--------------------------------------------------------------------------###

# Set the working directory to the project root
setwd("/Users/adamkurth/Documents/RStudio/research/bayes-priors-sens")

# Configuration
REGENERATE_DATA <- FALSE  # Set TRUE to regenerate data batches
RUN_FULL_ANALYSIS <- TRUE # Set TRUE to run the complete analysis

# 1. Generate data batches if needed
if (REGENERATE_DATA || !dir.exists("data_batches")) {
  print("=== Generating Data Batches ===")
  source("R/02_generate.r")
} else {
  print("=== Data batches already exist. Skipping generation. ===")
  print("=== Set REGENERATE_DATA <- TRUE to regenerate. ===")
}

# 2. Run the main analysis
if (RUN_FULL_ANALYSIS) {
  print("=== Running Full Analysis ===")
  print("=== This may take several hours for all 1000 simulations ===")
  print("=== Reduce N_SIMS_TO_FIT in 01_test.r for faster testing ===")
  source("R/01_test.r")
}

print("=== Pipeline Complete ===")
print("Check /results folder for outputs.")