###--------------------------------------------------------------------------###
### Quick Test Script: Verify Pipeline Works
###
### This script runs a minimal version of the analysis (1 simulation per batch)
### to verify the entire pipeline works before running the full analysis.
###--------------------------------------------------------------------------###

rm(list = ls())

# Set the working directory
setwd("/Users/adamkurth/Documents/RStudio/research/bayes-priors-sens")

library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)

# Load helpers
source("R/00_data.r")
source("R/utils/helpers.r")

# Create results directories
for (dir in c("results", "results/figures", "results/tables", "results/stan_diagnostics")) {
  if (!dir.exists(dir)) dir.create(dir)
}

# Parallel processing
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

###--------------------------------------------------------------------------###
### QUICK TEST: Run on ONE simulation from ONE batch
###--------------------------------------------------------------------------###

print("=== QUICK TEST: Single Simulation ===")

# Compile model
stan_model <- stan_model(file = "stan/partial_pool_continuous.stan")

# Load one batch
batch_file <- "data_batches/batch_BASELINE_Clean_N1000.rds"
batch_data <- readRDS(batch_file)
sim <- batch_data[[1]]

# Configure weak pooling prior
sim$stan.data$prior_tau_sigma <- 5
sim$stan.data$prior_mean_mu <- 0

# Fit model
print("Fitting model...")
fit <- sampling(
  stan_model,
  data = sim$stan.data,
  iter = 1000,  # Reduced for quick test
  warmup = 500,
  chains = 2,   # Reduced for quick test
  seed = 2025
)

# Check diagnostics
print("=== Diagnostics ===")
diagnostics <- extract_stan_diagnostics(fit)
cat(format_diagnostics_report(diagnostics))

# Extract CATE summary
print("=== CATE Estimates vs Truth ===")
cate_summary <- summary(fit, pars = "CATE")$summary %>% 
  as.data.frame() %>%
  rownames_to_column("param")

subgroup_map <- sim$truth.data %>%
  dplyr::select(any_of(c("subgroup_id", "subgroup_label", "true.cate"))) %>%
  distinct() %>%
  arrange(subgroup_id)

results <- cbind(subgroup_map, cate_summary)
print(results %>% dplyr::select(any_of(c("subgroup_label", "true.cate", "mean", "2.5%", "97.5%"))))

# Check if truth is within 95% CI
results <- results %>%
  mutate(covered = (true.cate >= `2.5%`) & (true.cate <= `97.5%`))

print(paste("\n95% Coverage:", mean(results$covered) * 100, "%"))

# Quick PPC
print("=== Generating PPC Plot ===")
y_rep <- as.matrix(fit, pars = "y_rep")
y_obs <- sim$stan.data$Y

p_ppc <- ppc_dens_overlay(y_obs, y_rep[1:50, ]) +
  ggtitle("Quick Test: Posterior Predictive Check")

ggsave("results/figures/quick_test_ppc.png", p_ppc, width = 8, height = 5)
print("PPC plot saved to results/figures/quick_test_ppc.png")

# Forest plot
print("=== Generating Forest Plot ===")
results$Sex <- ifelse(grepl("^F", results$subgroup_label), "Female", "Male")

p_forest <- ggplot(results, aes(x = subgroup_label, y = mean, color = Sex)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  geom_point(aes(y = true.cate), color = "black", shape = 4, size = 4, stroke = 1.5) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "Quick Test: CATE Forest Plot",
    subtitle = "Black X = True CATE, Dots = Estimates, Error bars = 95% CrI",
    x = "Subgroup",
    y = "Treatment Effect (MELD points)"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8"))

ggsave("results/figures/quick_test_forest.png", p_forest, width = 10, height = 6)
print("Forest plot saved to results/figures/quick_test_forest.png")

print("\n=== QUICK TEST COMPLETE ===")
print("If no errors occurred, run the full analysis with:")
print("  source('R/01_test.r')")
