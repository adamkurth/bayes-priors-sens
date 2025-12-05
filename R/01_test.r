###--------------------------------------------------------------------------###
### Project: Quantifying Uncertainty in Causal Effects
### Step 1: Model Evaluation, Posterior Summaries & Sensitivity Analysis
### 
### Objective: 
### 1. Fit the Partial Pooling Stan model across ALL data batches.
### 2. Compare Weak vs Strong Pooling priors on tau.
### 3. Extract posterior summaries, diagnostics, and PPCs.
### 4. Conduct sensitivity analysis across (gamma_U, beta_U) grid.
### 5. Save all results to /results for LaTeX write-up.
###--------------------------------------------------------------------------###
rm(list=ls())

library(tidyverse)
library(latex2exp)
library(rstan)
library(shinystan)
library(bayesplot)
library(loo)        # For model comparison (PSIS-LOO)
library(xtable)     # For LaTeX tables

# Set working directory
setwd("/Users/adamkurth/Documents/RStudio/research/bayes-priors-sens")

# Load data generating function 
source("R/00_data.r")

# Create results directory if needed
if(!dir.exists("results")) dir.create("results")
if(!dir.exists("results/figures")) dir.create("results/figures")
if(!dir.exists("results/tables")) dir.create("results/tables")
if(!dir.exists("results/stan_diagnostics")) dir.create("results/stan_diagnostics")

# Optimize for parallel processing
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Color scheme for plots
color_scheme_set("blue")

###--------------------------------------------------------------------------###
### CONFIGURATION
###--------------------------------------------------------------------------###

# Number of simulations to use per batch (use subset for speed, full = 1000)
N_SIMS_TO_FIT <- 5  # Set to 1000 for full analysis, 5 for testing

# Prior configurations to compare
PRIOR_CONFIGS <- list(
  weak_pooling = list(
    name = "Weak Pooling",
    tau_sigma = 5,    # Wide prior on heterogeneity
    mu_mean = 0       # Weakly informative prior on ATE
  ),
  strong_pooling = list(
    name = "Strong Pooling", 
    tau_sigma = 0.5,  # Tight prior -> shrinks CATEs toward global mean
    mu_mean = 0
  )
)

# Stan settings
STAN_ITER <- 2000
STAN_WARMUP <- 1000
STAN_CHAINS <- 4
STAN_SEED <- 2025

###--------------------------------------------------------------------------###
### HELPER FUNCTIONS
###--------------------------------------------------------------------------###

#' Fit Stan model and extract results
#' 
#' @param stan.data List of data for Stan
#' @param prior_config List with tau_sigma and mu_mean
#' @param stan_model Compiled Stan model object
#' @param truth_data Full data frame with true CATEs
fit_and_extract <- function(stan.data, prior_config, stan_model, truth_data) {
  
  # Apply prior configuration
  stan.data$prior_tau_sigma <- prior_config$tau_sigma
  stan.data$prior_mean_mu <- prior_config$mu_mean
  
  # Fit model
  fit <- sampling(
    stan_model,
    data = stan.data,
    iter = STAN_ITER,
    warmup = STAN_WARMUP,
    chains = STAN_CHAINS,
    seed = STAN_SEED,
    refresh = 0,  # Suppress iteration output
    control = list(max_treedepth = 12, adapt_delta = 0.9)  # Reduce warnings
  )
  
  # Extract diagnostics
  diagnostics <- list(
    rhat = summary(fit)$summary[, "Rhat"],
    n_eff = summary(fit)$summary[, "n_eff"],
    divergences = sum(get_divergent_iterations(fit)),
    max_treedepth = sum(get_max_treedepth_iterations(fit))
  )
  
  # Extract CATE summaries
  cate_summary <- summary(fit, pars = "CATE")$summary %>% 
    as.data.frame() %>%
    rownames_to_column("param")
  
 # Add subgroup labels - use dplyr::select explicitly and handle column names
  subgroup_map <- truth_data %>% 
    dplyr::select(any_of(c("subgroup_id", "subgroup_label", "true.cate", "sex", "ascites"))) %>%
    distinct() %>%
    arrange(subgroup_id)
  
  cate_summary <- cbind(subgroup_map, cate_summary)
  
  # Extract global parameters
  global_params <- summary(fit, pars = c("mu", "tau", "sigma"))$summary %>%
    as.data.frame() %>%
    rownames_to_column("param")
  
  # Posterior samples for detailed analysis
  posterior_cates <- as.matrix(fit, pars = "CATE")
  colnames(posterior_cates) <- subgroup_map$subgroup_label
  
  # Posterior predictive samples
  y_rep <- as.matrix(fit, pars = "y_rep")
  
  return(list(
    fit = fit,
    cate_summary = cate_summary,
    global_params = global_params,
    posterior_cates = posterior_cates,
    y_rep = y_rep,
    diagnostics = diagnostics
  ))
}

#' Calculate bias and coverage metrics
calc_performance <- function(cate_summary) {
  cate_summary %>%
    mutate(
      bias = mean - true.cate,
      abs_bias = abs(bias),
      coverage_95 = (true.cate >= `2.5%`) & (true.cate <= `97.5%`),
      ci_width = `97.5%` - `2.5%`,
      rmse = sqrt(bias^2)
    ) %>%
    summarise(
      mean_bias = mean(bias),
      mean_abs_bias = mean(abs_bias),
      coverage_rate = mean(coverage_95),
      mean_ci_width = mean(ci_width),
      rmse = sqrt(mean(bias^2))
    )
}

#' Calculate probability of benefit (CATE < 0)
calc_prob_benefit <- function(posterior_cates) {
  apply(posterior_cates, 2, function(x) mean(x < 0))
}

###--------------------------------------------------------------------------###
### PART 1: COMPILE STAN MODEL
###--------------------------------------------------------------------------###

print("=== Compiling Stan Model ===")
stan_model <- stan_model(file = "stan/partial_pool_continuous.stan")

###--------------------------------------------------------------------------###
### PART 2: FIT BASELINE (CLEAN) SCENARIOS
###--------------------------------------------------------------------------###

print("=== PART 2: Fitting Baseline (Clean) Scenarios ===")

baseline_results <- list()

for (n_size in c(200, 1000)) {
  
  batch_file <- paste0("data_batches/batch_BASELINE_Clean_N", n_size, ".rds")
  print(paste("Loading:", batch_file))
  batch_data <- readRDS(batch_file)
  
  for (prior_name in names(PRIOR_CONFIGS)) {
    
    prior_config <- PRIOR_CONFIGS[[prior_name]]
    result_key <- paste0("N", n_size, "_", prior_name)
    print(paste("  Fitting:", result_key))
    
    # Fit on multiple simulations and aggregate
    sim_results <- list()
    
    for (i in 1:min(N_SIMS_TO_FIT, length(batch_data))) {
      
      sim <- batch_data[[i]]
      
      tryCatch({
        res <- fit_and_extract(
          stan.data = sim$stan.data,
          prior_config = prior_config,
          stan_model = stan_model,
          truth_data = sim$truth.data
        )
        
        sim_results[[i]] <- list(
          cate_summary = res$cate_summary,
          global_params = res$global_params,
          performance = calc_performance(res$cate_summary),
          prob_benefit = calc_prob_benefit(res$posterior_cates),
          diagnostics = res$diagnostics
        )
        
        # Save first simulation fit for detailed diagnostics
        if (i == 1) {
          baseline_results[[result_key]]$first_fit <- res$fit
          baseline_results[[result_key]]$first_y_rep <- res$y_rep
          baseline_results[[result_key]]$first_y_obs <- sim$stan.data$Y
          baseline_results[[result_key]]$first_truth <- sim$truth.data
        }
        
      }, error = function(e) {
        print(paste("    Error in sim", i, ":", e$message))
      })
    }
    
    baseline_results[[result_key]]$sims <- sim_results
    baseline_results[[result_key]]$prior_config <- prior_config
    baseline_results[[result_key]]$n_size <- n_size
  }
}

###--------------------------------------------------------------------------###
### PART 3: FIT CONFOUNDED SCENARIOS (Sensitivity Analysis Grid)
###--------------------------------------------------------------------------###

print("=== PART 3: Fitting Confounded Scenarios (Sensitivity Grid) ===")

# Define the sensitivity grid
sensitivity_grid <- expand.grid(
  gamma = c(0.5, 2.0),
  beta = c(2, 4, 8),
  n_size = c(200, 1000)
)

confounded_results <- list()

for (row_idx in 1:nrow(sensitivity_grid)) {
  
  curr_gamma <- sensitivity_grid$gamma[row_idx]
  curr_beta <- sensitivity_grid$beta[row_idx]
  curr_n <- sensitivity_grid$n_size[row_idx]
  
  batch_file <- paste0("data_batches/batch_Confounded_N", curr_n, 
                       "_Gamma", curr_gamma, "_Beta", curr_beta, ".rds")
  
  if (!file.exists(batch_file)) {
    print(paste("Skipping (not found):", batch_file))
    next
  }
  
  print(paste("Loading:", batch_file))
  batch_data <- readRDS(batch_file)
  
  for (prior_name in names(PRIOR_CONFIGS)) {
    
    prior_config <- PRIOR_CONFIGS[[prior_name]]
    result_key <- paste0("N", curr_n, "_G", curr_gamma, "_B", curr_beta, "_", prior_name)
    print(paste("  Fitting:", result_key))
    
    sim_results <- list()
    
    for (i in 1:min(N_SIMS_TO_FIT, length(batch_data))) {
      
      sim <- batch_data[[i]]
      
      tryCatch({
        res <- fit_and_extract(
          stan.data = sim$stan.data,
          prior_config = prior_config,
          stan_model = stan_model,
          truth_data = sim$truth.data
        )
        
        # Calculate naive estimate for comparison
        d_truth <- sim$truth.data
        naive_ate <- mean(d_truth$Y.obs[d_truth$A==1]) - mean(d_truth$Y.obs[d_truth$A==0])
        true_ate <- mean(d_truth$true.cate)
        
        sim_results[[i]] <- list(
          cate_summary = res$cate_summary,
          global_params = res$global_params,
          performance = calc_performance(res$cate_summary),
          prob_benefit = calc_prob_benefit(res$posterior_cates),
          diagnostics = res$diagnostics,
          naive_ate = naive_ate,
          true_ate = true_ate,
          naive_bias = naive_ate - true_ate
        )
        
        # Save first fit for diagnostics
        if (i == 1) {
          confounded_results[[result_key]]$first_fit <- res$fit
          confounded_results[[result_key]]$first_posterior_cates <- res$posterior_cates
        }
        
      }, error = function(e) {
        print(paste("    Error in sim", i, ":", e$message))
      })
    }
    
    confounded_results[[result_key]]$sims <- sim_results
    confounded_results[[result_key]]$prior_config <- prior_config
    confounded_results[[result_key]]$gamma <- curr_gamma
    confounded_results[[result_key]]$beta <- curr_beta
    confounded_results[[result_key]]$n_size <- curr_n
  }
}

###--------------------------------------------------------------------------###
### PART 4: AGGREGATE RESULTS AND CREATE TABLES
###--------------------------------------------------------------------------###

print("=== PART 4: Aggregating Results ===")

# --- Table 1: Posterior Summaries for Baseline (Clean) Scenarios ---
baseline_summary <- do.call(rbind, lapply(names(baseline_results), function(key) {
  res <- baseline_results[[key]]
  if (is.null(res$sims) || length(res$sims) == 0) return(NULL)
  
  # Aggregate CATE estimates across simulations
  all_cate <- do.call(rbind, lapply(res$sims, function(s) s$cate_summary))
  
  agg <- all_cate %>%
    group_by(subgroup_label, sex, ascites) %>%
    summarise(
      true_cate = mean(true.cate),
      est_mean = mean(mean),
      est_sd = mean(sd),
      ci_lower = mean(`2.5%`),
      ci_upper = mean(`97.5%`),
      .groups = "drop"
    ) %>%
    mutate(
      prior = res$prior_config$name,
      n_size = res$n_size,
      scenario = "Clean (RCT)"
    )
  
  return(agg)
}))

# --- Table 2: Performance Metrics by Prior ---
performance_summary <- do.call(rbind, lapply(names(baseline_results), function(key) {
  res <- baseline_results[[key]]
  if (is.null(res$sims) || length(res$sims) == 0) return(NULL)
  
  perf <- do.call(rbind, lapply(res$sims, function(s) s$performance))
  
  data.frame(
    prior = res$prior_config$name,
    n_size = res$n_size,
    scenario = "Clean",
    mean_bias = mean(perf$mean_bias),
    rmse = mean(perf$rmse),
    coverage_95 = mean(perf$coverage_rate),
    ci_width = mean(perf$mean_ci_width)
  )
}))

# --- Table 3: Sensitivity Analysis Results ---
sensitivity_summary <- do.call(rbind, lapply(names(confounded_results), function(key) {
  res <- confounded_results[[key]]
  if (is.null(res$sims) || length(res$sims) == 0) return(NULL)
  
  # Aggregate performance
  perf <- do.call(rbind, lapply(res$sims, function(s) s$performance))
  
  # Aggregate probability of benefit by subgroup
  prob_ben <- do.call(rbind, lapply(res$sims, function(s) s$prob_benefit))
  avg_prob_ben <- colMeans(prob_ben, na.rm = TRUE)
  
  # Naive bias
  naive_bias <- mean(sapply(res$sims, function(s) s$naive_bias))
  
  data.frame(
    prior = res$prior_config$name,
    n_size = res$n_size,
    gamma = res$gamma,
    beta = res$beta,
    mean_bias = mean(perf$mean_bias),
    rmse = mean(perf$rmse),
    coverage_95 = mean(perf$coverage_rate),
    naive_bias = naive_bias,
    # Probability benefit remains for Female subgroups
    prob_benefit_female = mean(avg_prob_ben[grepl("^F_", names(avg_prob_ben))]),
    prob_benefit_male = mean(avg_prob_ben[grepl("^M_", names(avg_prob_ben))])
  )
}))

###--------------------------------------------------------------------------###
### PART 5: CREATE VALIDITY REGION ANALYSIS
###--------------------------------------------------------------------------###

print("=== PART 5: Validity Region Analysis ===")

# A subgroup's benefit is "robust" if P(CATE < 0 | data) > 0.95
# across all confounding scenarios

validity_analysis <- sensitivity_summary %>%
  group_by(gamma, beta, prior) %>%
  summarise(
    female_robust = mean(prob_benefit_female > 0.95),
    male_robust = mean(prob_benefit_male > 0.95),
    any_robust = mean(prob_benefit_female > 0.95 | prob_benefit_male > 0.95),
    .groups = "drop"
  )

# Identify regions where confounding overturns benefit
overturn_region <- sensitivity_summary %>%
  filter(prob_benefit_female < 0.5 | prob_benefit_male < 0.5) %>%
  select(gamma, beta, n_size, prior, prob_benefit_female, prob_benefit_male)

###--------------------------------------------------------------------------###
### PART 6: GENERATE FIGURES
###--------------------------------------------------------------------------###

print("=== PART 6: Generating Figures ===")

# --- Figure 1: Forest Plot (Weak vs Strong Pooling) ---
if (!is.null(baseline_summary) && nrow(baseline_summary) > 0) {
  
  p_forest_compare <- baseline_summary %>%
    filter(n_size == 1000) %>%
    mutate(Sex = ifelse(sex == 1, "Female", "Male")) %>%
    ggplot(aes(x = subgroup_label, y = est_mean, color = prior)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  position = position_dodge(width = 0.5), width = 0.2) +
    geom_point(aes(y = true_cate), color = "black", shape = 4, size = 4, stroke = 1.5) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(
      title = "CATE Estimates: Weak vs Strong Pooling Priors",
      subtitle = "Black X = True CATE. Error bars = 95% CrI. (N = 1000, Clean RCT)",
      x = "Subgroup",
      y = "Treatment Effect (Change in MELD Score)",
      color = "Prior"
    ) +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")
  
  ggsave("results/figures/fig1_forest_prior_comparison.pdf", 
         p_forest_compare, width = 10, height = 6)
  ggsave("results/figures/fig1_forest_prior_comparison.png", 
         p_forest_compare, width = 10, height = 6, dpi = 300)
}

# --- Figure 2: Posterior Predictive Check ---
# Use first baseline fit with N=1000, weak pooling
key_for_ppc <- "N1000_weak_pooling"
if (!is.null(baseline_results[[key_for_ppc]]$first_fit)) {
  
  y_rep <- baseline_results[[key_for_ppc]]$first_y_rep
  y_obs <- baseline_results[[key_for_ppc]]$first_y_obs
  
  p_ppc <- ppc_dens_overlay(y_obs, y_rep[1:50, ]) +
    ggtitle("Posterior Predictive Check",
            subtitle = "Black = Observed MELD Scores, Blue = Model Simulations (N=1000, Clean)")
  
  ggsave("results/figures/fig2_ppc_check.pdf", p_ppc, width = 8, height = 5)
  ggsave("results/figures/fig2_ppc_check.png", p_ppc, width = 8, height = 5, dpi = 300)
}

# --- Figure 3: Sensitivity Grid Heatmap ---
if (!is.null(sensitivity_summary) && nrow(sensitivity_summary) > 0) {
  
  # Heatmap: Probability of Female Benefit across (gamma, beta)
  p_sens_female <- sensitivity_summary %>%
    filter(prior == "Weak Pooling", n_size == 1000) %>%
    ggplot(aes(x = factor(gamma), y = factor(beta), fill = prob_benefit_female)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", prob_benefit_female)), color = "black", size = 4) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                         midpoint = 0.5, limits = c(0, 1),
                         name = "P(Benefit)") +
    labs(
      title = "Sensitivity Analysis: Female Subgroups",
      subtitle = "P(CATE < 0 | data) across confounding scenarios (N=1000, Weak Pooling)",
      x = expression(gamma[U] ~ "(Selection Bias)"),
      y = expression(beta[U] ~ "(Outcome Bias)")
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  ggsave("results/figures/fig3_sensitivity_female.pdf", p_sens_female, width = 7, height = 5)
  ggsave("results/figures/fig3_sensitivity_female.png", p_sens_female, width = 7, height = 5, dpi = 300)
  
  # Heatmap: Probability of Male Benefit
  p_sens_male <- sensitivity_summary %>%
    filter(prior == "Weak Pooling", n_size == 1000) %>%
    ggplot(aes(x = factor(gamma), y = factor(beta), fill = prob_benefit_male)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", prob_benefit_male)), color = "black", size = 4) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                         midpoint = 0.5, limits = c(0, 1),
                         name = "P(Benefit)") +
    labs(
      title = "Sensitivity Analysis: Male Subgroups",
      subtitle = "P(CATE < 0 | data) across confounding scenarios (N=1000, Weak Pooling)",
      x = expression(gamma[U] ~ "(Selection Bias)"),
      y = expression(beta[U] ~ "(Outcome Bias)")
    ) +
    theme_minimal()
  
  ggsave("results/figures/fig3_sensitivity_male.pdf", p_sens_male, width = 7, height = 5)
  ggsave("results/figures/fig3_sensitivity_male.png", p_sens_male, width = 7, height = 5, dpi = 300)
}

# --- Figure 4: Bias vs Confounding Strength ---
if (!is.null(sensitivity_summary) && nrow(sensitivity_summary) > 0) {
  
  p_bias_conf <- sensitivity_summary %>%
    filter(n_size == 1000) %>%
    mutate(confound_strength = gamma * beta) %>%
    ggplot(aes(x = confound_strength, y = naive_bias, color = prior, shape = prior)) +
    geom_point(size = 3) +
    geom_line(aes(group = prior)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = "Naive Bias vs Confounding Strength",
      subtitle = "Bias = Naive ATE - True ATE (N = 1000)",
      x = expression(gamma[U] %*% beta[U] ~ "(Confounding Product)"),
      y = "Bias in Naive Estimate (MELD points)"
    ) +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")
  
  ggsave("results/figures/fig4_bias_vs_confounding.pdf", p_bias_conf, width = 8, height = 5)
  ggsave("results/figures/fig4_bias_vs_confounding.png", p_bias_conf, width = 8, height = 5, dpi = 300)
}

# --- Figure 5: Sex Disparity Posterior Densities ---
if (!is.null(baseline_results[["N1000_weak_pooling"]]$first_fit)) {
  
  posterior_cates <- as.matrix(baseline_results[["N1000_weak_pooling"]]$first_fit, pars = "CATE")
  truth_data <- baseline_results[["N1000_weak_pooling"]]$first_truth
  
  subgroup_labels <- truth_data %>%
    dplyr::select(any_of(c("subgroup_id", "subgroup_label"))) %>%
    distinct() %>%
    arrange(subgroup_id)
  
  colnames(posterior_cates) <- subgroup_labels$subgroup_label
  
  d_post_long <- as.data.frame(posterior_cates) %>%
    pivot_longer(everything(), names_to = "Subgroup", values_to = "CATE") %>%
    mutate(Sex = ifelse(grepl("^F", Subgroup), "Female", "Male"))
  
  p_disparity <- ggplot(d_post_long, aes(x = CATE, fill = Sex)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = c(-5, -3), linetype = "dashed", color = c("red", "blue")) +
    annotate("text", x = -4.7, y = 0.3, label = "Female Truth (-5)", angle = 90, color = "red") +
    annotate("text", x = -2.7, y = 0.3, label = "Male Truth (-3)", angle = 90, color = "blue") +
    labs(
      title = "Sex Disparity in Treatment Effect",
      subtitle = "Posterior densities show Females benefit more than Males (N=1000, Clean)",
      x = "Estimated Change in MELD Score (CATE)",
      y = "Posterior Density"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8"))
  
  ggsave("results/figures/fig5_sex_disparity.pdf", p_disparity, width = 8, height = 5)
  ggsave("results/figures/fig5_sex_disparity.png", p_disparity, width = 8, height = 5, dpi = 300)
}

###--------------------------------------------------------------------------###
### PART 7: GENERATE LATEX TABLES
###--------------------------------------------------------------------------###

print("=== PART 7: Generating LaTeX Tables ===")

# --- Table 1: CATE Posterior Summaries ---
if (!is.null(baseline_summary) && nrow(baseline_summary) > 0) {
  
  table1 <- baseline_summary %>%
    filter(n_size == 1000) %>%
    select(Subgroup = subgroup_label, Prior = prior, 
           `True CATE` = true_cate, 
           Mean = est_mean, 
           `95% CI Lower` = ci_lower,
           `95% CI Upper` = ci_upper) %>%
    mutate(across(where(is.numeric), ~round(., 2)))
  
  print(xtable(table1, 
               caption = "Posterior summaries of subgroup-specific CATEs under weak and strong pooling priors (N=1000, Clean RCT scenario).",
               label = "tab:cate_summaries"),
        file = "results/tables/table1_cate_summaries.tex",
        include.rownames = FALSE)
}

# --- Table 2: Performance Comparison ---
if (!is.null(performance_summary) && nrow(performance_summary) > 0) {
  
  table2 <- performance_summary %>%
    select(Prior = prior, N = n_size, 
           `Mean Bias` = mean_bias, 
           RMSE = rmse,
           `95% Coverage` = coverage_95,
           `CI Width` = ci_width) %>%
    mutate(across(where(is.numeric), ~round(., 3)))
  
  print(xtable(table2,
               caption = "Performance metrics for CATE estimation under different priors and sample sizes.",
               label = "tab:performance"),
        file = "results/tables/table2_performance.tex",
        include.rownames = FALSE)
}

# --- Table 3: Sensitivity Analysis Grid ---
if (!is.null(sensitivity_summary) && nrow(sensitivity_summary) > 0) {
  
  table3 <- sensitivity_summary %>%
    filter(prior == "Weak Pooling") %>%
    select(`N` = n_size,
           `$\\gamma_U$` = gamma, 
           `$\\beta_U$` = beta,
           `Naive Bias` = naive_bias,
           `Model Bias` = mean_bias,
           `P(Benefit|Female)` = prob_benefit_female,
           `P(Benefit|Male)` = prob_benefit_male) %>%
    mutate(across(where(is.numeric), ~round(., 2)))
  
  print(xtable(table3,
               caption = "Sensitivity analysis: Posterior probability of benefit across confounding scenarios (Weak Pooling prior).",
               label = "tab:sensitivity"),
        file = "results/tables/table3_sensitivity.tex",
        include.rownames = FALSE,
        sanitize.colnames.function = identity)
}

# --- Table 4: Validity Regions ---
if (!is.null(validity_analysis) && nrow(validity_analysis) > 0) {
  
  table4 <- validity_analysis %>%
    select(`$\\gamma_U$` = gamma,
           `$\\beta_U$` = beta,
           Prior = prior,
           `Female Robust` = female_robust,
           `Male Robust` = male_robust) %>%
    mutate(across(where(is.numeric), ~round(., 2)))
  
  print(xtable(table4,
               caption = "Validity regions: Proportion of scenarios where P(CATE < 0) > 0.95 (robust benefit).",
               label = "tab:validity"),
        file = "results/tables/table4_validity.tex",
        include.rownames = FALSE,
        sanitize.colnames.function = identity)
}

###--------------------------------------------------------------------------###
### PART 8: SAVE STAN DIAGNOSTICS
###--------------------------------------------------------------------------###

print("=== PART 8: Saving Stan Diagnostics ===")

# Diagnostics summary
diag_summary <- list()

for (key in names(baseline_results)) {
  if (!is.null(baseline_results[[key]]$sims)) {
    diags <- lapply(baseline_results[[key]]$sims, function(s) s$diagnostics)
    
    diag_summary[[key]] <- list(
      max_rhat = max(sapply(diags, function(d) max(d$rhat, na.rm = TRUE))),
      min_neff = min(sapply(diags, function(d) min(d$n_eff, na.rm = TRUE))),
      total_divergences = sum(sapply(diags, function(d) d$divergences)),
      total_treedepth = sum(sapply(diags, function(d) d$max_treedepth))
    )
  }
}

saveRDS(diag_summary, "results/stan_diagnostics/diagnostics_summary.rds")

# Write diagnostics report
sink("results/stan_diagnostics/diagnostics_report.txt")
cat("=== Stan Diagnostics Report ===\n\n")
cat("Generated:", as.character(Sys.time()), "\n\n")

for (key in names(diag_summary)) {
  cat(paste0("--- ", key, " ---\n"))
  cat(paste0("  Max R-hat: ", round(diag_summary[[key]]$max_rhat, 4), "\n"))
  cat(paste0("  Min n_eff: ", round(diag_summary[[key]]$min_neff, 0), "\n"))
  cat(paste0("  Total divergences: ", diag_summary[[key]]$total_divergences, "\n"))
  cat(paste0("  Treedepth warnings: ", diag_summary[[key]]$total_treedepth, "\n\n"))
}
sink()

###--------------------------------------------------------------------------###
### PART 9: SAVE ALL RESULTS
###--------------------------------------------------------------------------###

print("=== PART 9: Saving All Results ===")

# Save aggregated results
saveRDS(baseline_summary, "results/baseline_summary.rds")
saveRDS(performance_summary, "results/performance_summary.rds")
saveRDS(sensitivity_summary, "results/sensitivity_summary.rds")
saveRDS(validity_analysis, "results/validity_analysis.rds")

# Save raw results (warning: can be large)
saveRDS(baseline_results, "results/baseline_results_raw.rds")
saveRDS(confounded_results, "results/confounded_results_raw.rds")

print("=== ANALYSIS COMPLETE ===")
print(paste("Results saved to:", normalizePath("results")))
print("Files generated:")
print(list.files("results", recursive = TRUE))

###--------------------------------------------------------------------------###
### OPTIONAL: Launch Shiny Stan for Interactive Exploration
###--------------------------------------------------------------------------###

# Uncomment to launch interactive diagnostics for a specific fit
# launch_shinystan(baseline_results[["N1000_weak_pooling"]]$first_fit)
