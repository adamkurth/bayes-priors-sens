###--------------------------------------------------------------------------###
### Project: Quantifying Uncertainty in Causal Effects
### Step: Batch Data Generation (Factorial Design)
### 
### Objective: 
### Generate batches of data across a GRID of simulation parameters.
### This ensures we test the model under varying degrees of "Confoundedness"
### AND varying degrees of "Sparsity" (Sample Size).
###--------------------------------------------------------------------------###

library(tidyverse)

# Load the updated flexible function
source("R/00_data.r") 

# setwd("/Users/adamkurth/Documents/RStudio/research/bayes-priors-sens")

if(!dir.exists("data_batches")) dir.create("data_batches")

###--------------------------------------------------------------------------###
### 1. Simulation Parameters
###--------------------------------------------------------------------------###
S <- 1000       # Number of simulations per scenario
SEED.BASE <- 2025

# --- FACTOR 1: SAMPLE SIZE (Crucial for Part 1 - Priors) ---
# N=200: High Sparsity. Some subgroups will have <10 patients. 
#        (Priors should matter A LOT here).
# N=1000: Moderate/Rich Data. 
#         (Priors should matter less, Data dominates).
sample.sizes <- c(200, 1000)

# --- FACTOR 2: CONFOUNDING (Crucial for Part 2 - Sensitivity) ---
# gamma: Selection Bias (U -> A)
# beta: Outcome Bias (U -> Y)
confounding.grid <- expand.grid(
  gamma.u = c(0.5, 2.0), # Weak vs Strong Selection
  beta.u  = c(2, 4, 8)   # Minor vs Moderate vs Severe Outcome Bias
)

print("--- Simulation Parameters Defined ---")
print(paste("Sample Sizes:", paste(sample.sizes, collapse=", ")))
print("Confounding Scenarios:")
print(confounding.grid)


###--------------------------------------------------------------------------###
### 2. Generate "Clean" Baselines (Randomized / No Unmeasured Confounding)
###--------------------------------------------------------------------------###
# We need these to prove the model works when assumptions hold.

for(n in sample.sizes) {
  
  batch.id <- paste0("BASELINE_Clean_N", n)
  print(paste(">>>> Generating Batch:", batch.id, "----------------"))
  
  sim.batch <- list()
  for(i in 1:S) {
    iter.seed <- SEED.BASE + n + i
    
    sim <- simulate.liver.data(
      N = n, 
      scenario = "clean", # Randomized
      seed = iter.seed
    )
    
    sim.batch[[i]] <- list(
      iter = i, seed = iter.seed,
      truth.data = sim$d, stan.data = sim$stan.data
    )
  }
  saveRDS(sim.batch, file = paste0("data_batches/batch_", batch.id, ".rds"))
}


###--------------------------------------------------------------------------###
### 3. Generate "Confounded" Scenarios (The Grid)
###--------------------------------------------------------------------------###

for(n in sample.sizes) {
  for(row in 1:nrow(confounding.grid)) {
    
    curr.gamma <- confounding.grid$gamma.u[row]
    curr.beta  <- confounding.grid$beta.u[row]
    
    batch.id <- paste0("Confounded_N", n, "_Gamma", curr.gamma, "_Beta", curr.beta)
    print(paste(">>>> Generating Batch:", batch.id, "----------------"))
    
    sim.batch <- list()
    for(i in 1:S) {
      # Unique seed logic to avoid overlap
      iter.seed <- SEED.BASE + (n * 10) + (row * 1000) + i
      
      # NOTE: Using underscores (gamma.u) to match 00.data.r
      sim <- simulate.liver.data(
        N = n, 
        scenario = "confounded", 
        gamma.u = curr.gamma, 
        beta.u = curr.beta,
        seed = iter.seed
      )
      
      sim.batch[[i]] <- list(
        iter = i,
        seed = iter.seed,
        params = list(N=n, gamma=curr.gamma, beta=curr.beta),
        truth.data = sim$d, 
        stan.data = sim$stan.data
      )
    }
    
    saveRDS(sim.batch, file = paste0("data_batches/batch_", batch.id, ".rds"))
    
    # --- Narrative Check (Validation) ---
    d.check <- sim.batch[[1]]$truth.data
    naive.diff <- mean(d.check$Y.obs[d.check$A==1]) - mean(d.check$Y.obs[d.check$A==0])
    true.ate   <- mean(d.check$true.cate)
    
    # Calculate Bias (Naive - Truth)
    bias <- naive.diff - true.ate
    
    print(paste("   N:", n, "| Gamma:", curr.gamma, "| Beta:", curr.beta,
                "| Bias:", round(bias, 2)))
  }
}

print("--- All Batches Complete ---")