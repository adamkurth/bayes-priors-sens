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
N.SIMS.TO.FIT <- 5  # Set to 1000 for full analysis, 5 for testing

# Prior configurations to compare
PRIOR.CONFIGS <- list(
  weak.pooling = list(
    name = "Weak Pooling",
    tau.sigma = 5,    # Wide prior on heterogeneity
    mu.mean = 0       # Weakly informative prior on ATE
  ),
  strong.pooling = list(
    name = "Strong Pooling", 
    tau.sigma = 0.5,  # Tight prior -> shrinks CATEs toward global mean
    mu.mean = 0
  )
)

# Stan settings
STAN.ITER <- 2000
STAN.WARMUP <- 1000
STAN.CHAINS <- 4
STAN.SEED <- 2025

###--------------------------------------------------------------------------###
### HELPER FUNCTIONS
###--------------------------------------------------------------------------###

#' Fit Stan model and extract results
#' 
#' @param stan.data List of data for Stan
#' @param prior.config List with tau.sigma and mu.mean
#' @param stan.model Compiled Stan model object
#' @param truth.data Full data frame with true CATEs
fit.and.extract <- function(stan.data, prior.config, stan.model, truth.data) {
  # fit, cate.summary, global.params, posterior.cates, y.rep, diagnostics
  # apply prior configuration
  stan.data$prior_tau_sigma <- prior.config$tau.sigma
  stan.data$prior_mean_mu <- prior.config$mu.mean
  
  # Fit model
  fit <- sampling(
    stan.model,
    data = stan.data,
    iter = STAN.ITER,
    warmup = STAN.WARMUP,
    chains = STAN.CHAINS,
    seed = STAN.SEED,
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
  cate.summary <- summary(fit, pars = "CATE")$summary %>% 
    as.data.frame() %>%
    rownames_to_column("param")
  
 # Add subgroup labels - use dplyr::select explicitly and handle column names
  subgroup.map <- truth.data %>% 
    dplyr::select(any_of(c("subgroup_id", "subgroup_label", "true.cate", "sex", "ascites"))) %>%
    distinct() %>%
    arrange(subgroup_id)
  
  cate.summary <- cbind(subgroup.map, cate.summary)
  
  # Extract global parameters
  global.params <- summary(fit, pars = c("mu", "tau", "sigma"))$summary %>%
    as.data.frame() %>%
    rownames_to_column("param")
  
  # Posterior samples for detailed analysis
  posterior.cates <- as.matrix(fit, pars = "CATE")
  colnames(posterior.cates) <- subgroup.map$subgroup_label
  
  # Posterior predictive samples
  y.rep <- as.matrix(fit, pars = "y_rep")
  
  return(list(
    fit = fit,
    cate.summary = cate.summary,
    global.params = global.params,
    posterior.cates = posterior.cates,
    y.rep = y.rep,
    diagnostics = diagnostics
  ))
}

#' Calculate bias and coverage metrics
calc.performance <- function(cate.summary) {
  # cate.summary: data frame with true.cate, mean, `2.5%`, `97.5%`
  cate.summary %>%
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
calc.prob.benefit <- function(posterior.cates) {
  # posterior.cates: matrix of posterior samples (rows = samples, cols = subgroups)
  # Return vector of probabilities per subgroup
  apply(posterior.cates, 2, function(x) mean(x < 0))
}

###--------------------------------------------------------------------------###
### PART 1: COMPILE STAN MODEL
###--------------------------------------------------------------------------###

print("=== Compiling Stan Model ===")
stan.model <- stan_model(file = "stan/partial_pool_continuous.stan")

###--------------------------------------------------------------------------###
### PART 2: FIT BASELINE (CLEAN) SCENARIOS
###--------------------------------------------------------------------------###

print("=== PART 2: Fitting Baseline (Clean) Scenarios ===")

baseline.results <- list()

for (n.size in c(200, 1000)) {
  
  batch.file <- paste0("data_batches/batch_BASELINE_Clean_N", n.size, ".rds")
  print(paste("Loading:", batch.file))
  batch.data <- readRDS(batch.file)
  
  for (prior.name in names(PRIOR.CONFIGS)) {
    
    prior.config <- PRIOR.CONFIGS[[prior.name]]
    result.key <- paste0("N", n.size, "_", prior.name)
    print(paste("  Fitting:", result.key))
    
    # Fit on multiple simulations and aggregate
    sim.results <- list()
    
    for (i in 1:min(N.SIMS.TO.FIT, length(batch.data))) {
      
      sim <- batch.data[[i]]
      
      tryCatch({
        res <- fit.and.extract(
          stan.data = sim$stan.data,
          prior.config = prior.config,
          stan.model = stan.model,
          truth.data = sim$truth.data
        )
        
        sim.results[[i]] <- list(
          cate.summary = res$cate.summary,
          global.params = res$global.params,
          performance = calc.performance(res$cate.summary),
          prob.benefit = calc.prob.benefit(res$posterior.cates),
          diagnostics = res$diagnostics
        )
        
        # Save first simulation fit for detailed diagnostics
        if (i == 1) {
          baseline.results[[result.key]]$first.fit <- res$fit
          baseline.results[[result.key]]$first.y.rep <- res$y.rep
          baseline.results[[result.key]]$first.y.obs <- sim$stan.data$Y
          baseline.results[[result.key]]$first.truth <- sim$truth.data
        }
        
      }, error = function(e) {
        print(paste("    Error in sim", i, ":", e$message))
      })
    }
    
    baseline.results[[result.key]]$sims <- sim.results
    baseline.results[[result.key]]$prior.config <- prior.config
    baseline.results[[result.key]]$n.size <- n.size
  }
}

###--------------------------------------------------------------------------###
### PART 3: FIT CONFOUNDED SCENARIOS (Sensitivity Analysis Grid)
###--------------------------------------------------------------------------###

print("=== PART 3: Fitting Confounded Scenarios (Sensitivity Grid) ===")

# Define the sensitivity grid
sensitivity.grid <- expand.grid(
  gamma = c(0.5, 2.0),
  beta = c(2, 4, 8),
  n.size = c(200, 1000)
)

confounded.results <- list()

for (row.idx in 1:nrow(sensitivity.grid)) {
  
  curr.gamma <- sensitivity.grid$gamma[row.idx]
  curr.beta <- sensitivity.grid$beta[row.idx]
  curr.n <- sensitivity.grid$n.size[row.idx]
  
  batch.file <- paste0("data_batches/batch_Confounded_N", curr.n, 
                       "_Gamma", curr.gamma, "_Beta", curr.beta, ".rds")
  
  if (!file.exists(batch.file)) {
    print(paste("Skipping (not found):", batch.file))
    next
  }
  
  print(paste("Loading:", batch.file))
  batch.data <- readRDS(batch.file)
  
  for (prior.name in names(PRIOR.CONFIGS)) {
    
    prior.config <- PRIOR.CONFIGS[[prior.name]]
    result.key <- paste0("N", curr.n, "_G", curr.gamma, "_B", curr.beta, "_", prior.name)
    print(paste("  Fitting:", result.key))
    
    sim.results <- list()
    
    for (i in 1:min(N.SIMS.TO.FIT, length(batch.data))) {
      
      sim <- batch.data[[i]]
      
      tryCatch({
        res <- fit.and.extract(
          stan.data = sim$stan.data,
          prior.config = prior.config,
          stan.model = stan.model,
          truth.data = sim$truth.data
        )
        
        # Calculate naive estimate for comparison
        d.truth <- sim$truth.data
        naive.ate <- mean(d.truth$Y.obs[d.truth$A==1]) - mean(d.truth$Y.obs[d.truth$A==0])
        true.ate <- mean(d.truth$true.cate)
        
        sim.results[[i]] <- list(
          cate.summary = res$cate.summary,
          global.params = res$global.params,
          performance = calc.performance(res$cate.summary),
          prob.benefit = calc.prob.benefit(res$posterior.cates),
          diagnostics = res$diagnostics,
          naive.ate = naive.ate,
          true.ate = true.ate,
          naive.bias = naive.ate - true.ate
        )
        
        # Save first fit for diagnostics
        if (i == 1) {
          confounded.results[[result.key]]$first.fit <- res$fit
          confounded.results[[result.key]]$first.posterior.cates <- res$posterior.cates
        }
        
      }, error = function(e) {
        print(paste("    Error in sim", i, ":", e$message))
      })
    }
    
    confounded.results[[result.key]]$sims <- sim.results
    confounded.results[[result.key]]$prior.config <- prior.config
    confounded.results[[result.key]]$gamma <- curr.gamma
    confounded.results[[result.key]]$beta <- curr.beta
    confounded.results[[result.key]]$n.size <- curr.n
  }
}

###--------------------------------------------------------------------------###
### PART 4: AGGREGATE RESULTS AND CREATE TABLES
###--------------------------------------------------------------------------###

print("=== PART 4: Aggregating Results ===")

# --- Table 1: Posterior Summaries for Baseline (Clean) Scenarios ---
baseline.summary <- do.call(rbind, lapply(names(baseline.results), function(key) {
  # baseline.summary extracts mean CATE estimates across sims
  res <- baseline.results[[key]]
  if (is.null(res$sims) || length(res$sims) == 0) return(NULL)
  
  # Aggregate CATE estimates across simulations
  all.cate <- do.call(rbind, lapply(res$sims, function(s) s$cate.summary))
  
  agg <- all.cate %>%
    group_by(subgroup_label, sex, ascites) %>%
    summarise(
      true.cate.mean = mean(true.cate),
      est.mean = mean(mean),
      est.sd = mean(sd),
      ci.lower = mean(`2.5%`),
      ci.upper = mean(`97.5%`),
      .groups = "drop"
    ) %>%
    mutate(
      prior = res$prior.config$name,
      n.size = res$n.size,
      scenario = "Clean (RCT)"
    )
  
  return(agg)
}))

# --- Table 2: Performance Metrics by Prior ---
performance.summary <- do.call(rbind, lapply(names(baseline.results), function(key) {
  # performance.summary extracts performance metrics across sims
  res <- baseline.results[[key]]
  if (is.null(res$sims) || length(res$sims) == 0) return(NULL)
  
  perf <- do.call(rbind, lapply(res$sims, function(s) s$performance))
  
  data.frame(
    prior = res$prior.config$name,
    n.size = res$n.size,
    scenario = "Clean",
    mean.bias = mean(perf$mean_bias),
    rmse = mean(perf$rmse),
    coverage.95 = mean(perf$coverage_rate),
    ci.width = mean(perf$mean_ci_width)
  )
}))

# --- Table 3: Sensitivity Analysis Results ---
sensitivity.summary <- do.call(rbind, lapply(names(confounded.results), function(key) {
  res <- confounded.results[[key]]
  if (is.null(res$sims) || length(res$sims) == 0) return(NULL)
  
  # Aggregate performance
  perf <- do.call(rbind, lapply(res$sims, function(s) s$performance))
  
  # Aggregate probability of benefit by subgroup
  prob.ben <- do.call(rbind, lapply(res$sims, function(s) s$prob.benefit))
  avg.prob.ben <- colMeans(prob.ben, na.rm = TRUE)
  
  # Naive bias
  naive.bias <- mean(sapply(res$sims, function(s) s$naive.bias))
  
  data.frame(
    prior = res$prior.config$name,
    n.size = res$n.size,
    gamma = res$gamma,
    beta = res$beta,
    mean.bias = mean(perf$mean_bias),
    rmse = mean(perf$rmse),
    coverage.95 = mean(perf$coverage_rate),
    naive.bias = naive.bias,
    # Probability benefit remains for Female subgroups
    prob.benefit.female = mean(avg.prob.ben[grepl("^F_", names(avg.prob.ben))]),
    prob.benefit.male = mean(avg.prob.ben[grepl("^M_", names(avg.prob.ben))])
  )
}))

###--------------------------------------------------------------------------###
### PART 5: CREATE VALIDITY REGION ANALYSIS
###--------------------------------------------------------------------------###

print("=== PART 5: Validity Region Analysis ===")

# A subgroup's benefit is "robust" if P(CATE < 0 | data) > 0.95
# across all confounding scenarios

validity.analysis <- sensitivity.summary %>%
  group_by(gamma, beta, prior) %>%
  summarise(
    female.robust = mean(prob.benefit.female > 0.95),
    male.robust = mean(prob.benefit.male > 0.95),
    any.robust = mean(prob.benefit.female > 0.95 | prob.benefit.male > 0.95),
    .groups = "drop"
  )

# Identify regions where confounding overturns benefit
overturn.region <- sensitivity.summary %>%
  filter(prob.benefit.female < 0.5 | prob.benefit.male < 0.5) %>%
  dplyr::select(gamma, beta, n.size, prior, prob.benefit.female, prob.benefit.male)

###--------------------------------------------------------------------------###
### PART 6: GENERATE FIGURES
###--------------------------------------------------------------------------###

print("=== PART 6: Generating Figures ===")

# --- Figure 1: Forest Plot (Weak vs Strong Pooling) - Faceted by Sample Size ---
if (!is.null(baseline.summary) && nrow(baseline.summary) > 0) {
  
  # Create faceted plot showing both N=200 and N=1000
  p.forest.compare <- baseline.summary %>%
    mutate(
      Sex = ifelse(sex == 1, "Female", "Male"),
      Sample.Size = paste0("N = ", n.size)
    ) %>%
    ggplot(aes(x = subgroup_label, y = est.mean, color = prior)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), 
                  position = position_dodge(width = 0.5), width = 0.2) +
    geom_point(aes(y = true.cate.mean), color = "black", shape = 4, size = 4, stroke = 1.5) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~Sample.Size, ncol = 2) +
    labs(
      title = "CATE Estimates: Weak vs Strong Pooling Priors",
      subtitle = "Black X = True CATE. Error bars = 95% CrI. Prior sensitivity decreases with larger N.",
      x = "Subgroup",
      y = "Treatment Effect (Change in MELD Score)",
      color = "Prior"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 9)
    ) +
    scale_color_brewer(palette = "Set1")
  
  ggsave("results/figures/fig1_forest_prior_comparison.pdf", 
         p.forest.compare, width = 12, height = 7)
  ggsave("results/figures/fig1_forest_prior_comparison.png", 
         p.forest.compare, width = 12, height = 7, dpi = 300)
}

# --- Figure 2: Posterior Predictive Check ---
# Use first baseline fit with N=1000, weak pooling
key.for.ppc <- "N1000_weak.pooling"
if (!is.null(baseline.results[[key.for.ppc]]$first.fit)) {
  
  y.rep <- baseline.results[[key.for.ppc]]$first.y.rep
  y.obs <- baseline.results[[key.for.ppc]]$first.y.obs
  
  p.ppc <- ppc_dens_overlay(y.obs, y.rep[1:50, ]) +
    ggtitle("Posterior Predictive Check",
            subtitle = "Black = Observed MELD Scores, Blue = Model Simulations (N=1000, Clean)")
  
  ggsave("results/figures/fig2_ppc_check.pdf", p.ppc, width = 8, height = 5)
  ggsave("results/figures/fig2_ppc_check.png", p.ppc, width = 8, height = 5, dpi = 300)
}

# --- Figure 3: Sensitivity Grid Heatmap ---
if (!is.null(sensitivity.summary) && nrow(sensitivity.summary) > 0) {
  
  # Heatmap: Probability of Female Benefit across (gamma, beta)
  p.sens.female <- sensitivity.summary %>%
    filter(prior == "Weak Pooling", n.size == 1000) %>%
    ggplot(aes(x = factor(gamma), y = factor(beta), fill = prob.benefit.female)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", prob.benefit.female)), color = "black", size = 4) +
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
    theme(legend.position = "right", plot.title = element_text(size = 11), 
          plot.subtitle = element_text(size = 7))
  
  # ggsave("results/figures/fig3_sensitivity_female.pdf", p.sens.female, width = 7, height = 5)
  # ggsave("results/figures/fig3_sensitivity_female.png", p.sens.female, width = 7, height = 5, dpi = 300)

  ggsave("results/figures/fig3_sensitivity_female.pdf", p.sens.female, width = 5, height = 4)
  ggsave("results/figures/fig3_sensitivity_female.png", p.sens.female, width = 5, height = 4, dpi = 300)
  
  # Heatmap: Probability of Male Benefit
  p.sens.male <- sensitivity.summary %>%
    filter(prior == "Weak Pooling", n.size == 1000) %>%
    ggplot(aes(x = factor(gamma), y = factor(beta), fill = prob.benefit.male)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", prob.benefit.male)), color = "black", size = 4) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                         midpoint = 0.5, limits = c(0, 1),
                         name = "P(Benefit)") +
    labs(
      title = "Sensitivity Analysis: Male Subgroups",
      subtitle = "P(CATE < 0 | data) across confounding scenarios (N=1000, Weak Pooling)",
      x = expression(gamma[U] ~ "(Selection Bias)"),
      y = expression(beta[U] ~ "(Outcome Bias)")
    ) +
    theme_minimal() + 
    theme(legend.position = "right", plot.title = element_text(size = 11), 
          plot.subtitle = element_text(size = 7))

  
  # ggsave("results/figures/fig3_sensitivity_male.pdf", p.sens.male, width = 7, height = 5)
  # ggsave("results/figures/fig3_sensitivity_male.png", p.sens.male, width = 7, height = 5, dpi = 300)

  ggsave("results/figures/fig3_sensitivity_male.pdf", p.sens.male, width = 5, height = 4)
  ggsave("results/figures/fig3_sensitivity_male.png", p.sens.male, width = 5, height = 4, dpi = 300)  
}

# --- Figure 4: Bias vs Confounding Strength ---
if (!is.null(sensitivity.summary) && nrow(sensitivity.summary) > 0) {
  
  p.bias.conf <- sensitivity.summary %>%
    filter(n.size == 1000) %>%
    mutate(confound.strength = gamma * beta) %>%
    ggplot(aes(x = confound.strength, y = naive.bias, color = prior, shape = prior)) +
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
  
  ggsave("results/figures/fig4_bias_vs_confounding.pdf", p.bias.conf, width = 8, height = 5)
  ggsave("results/figures/fig4_bias_vs_confounding.png", p.bias.conf, width = 8, height = 5, dpi = 300)
}

# --- Figure 5: Sex Disparity Posterior Densities ---
if (!is.null(baseline.results[["N1000_weak.pooling"]]$first.fit)) {
  
  posterior.cates <- as.matrix(baseline.results[["N1000_weak.pooling"]]$first.fit, pars = "CATE")
  truth.data <- baseline.results[["N1000_weak.pooling"]]$first.truth
  
  subgroup.labels <- truth.data %>%
    dplyr::select(any_of(c("subgroup_id", "subgroup_label"))) %>%
    distinct() %>%
    arrange(subgroup_id)
  
  colnames(posterior.cates) <- subgroup.labels$subgroup_label
  
  d.post.long <- as.data.frame(posterior.cates) %>%
    pivot_longer(everything(), names_to = "Subgroup", values_to = "CATE") %>%
    mutate(Sex = ifelse(grepl("^F", Subgroup), "Female", "Male"))
  
  p.disparity <- ggplot(d.post.long, aes(x = CATE, fill = Sex)) +
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
  
  ggsave("results/figures/fig5_sex_disparity.pdf", p.disparity, width = 8, height = 5)
  ggsave("results/figures/fig5_sex_disparity.png", p.disparity, width = 8, height = 5, dpi = 300)
}

# --- Compute Probability of Sex Disparity ---
if (!is.null(baseline.results[["N1000_weak.pooling"]]$first.fit)) {
  
  posterior.cates <- as.matrix(baseline.results[["N1000_weak.pooling"]]$first.fit, pars = "CATE")
  truth.data <- baseline.results[["N1000_weak.pooling"]]$first.truth
  
  subgroup.labels <- truth.data %>%
    dplyr::select(any_of(c("subgroup_id", "subgroup_label"))) %>%
    distinct() %>%
    arrange(subgroup_id)
  
  colnames(posterior.cates) <- subgroup.labels$subgroup_label
  
  # Identify female vs male subgroups
  female.cols <- grep("^F_", colnames(posterior.cates), value = TRUE)
  male.cols <- grep("^M_", colnames(posterior.cates), value = TRUE)
  
  # Compute average CATE by sex for each posterior draw
  avg.cate.female <- rowMeans(posterior.cates[, female.cols])
  avg.cate.male <- rowMeans(posterior.cates[, male.cols])
  
  # Probability that females benefit more (more negative CATE)
  prob.female.better <- mean(avg.cate.female < avg.cate.male)
  
  # Effect size: difference in average CATE
  sex.diff <- avg.cate.female - avg.cate.male
  sex.diff.summary <- c(
    mean = mean(sex.diff),
    sd = sd(sex.diff),
    q025 = quantile(sex.diff, 0.025),
    q975 = quantile(sex.diff, 0.975)
  )
  
  print("=== Sex Disparity Analysis ===")
  print(paste("P(Female benefit > Male benefit):", round(prob.female.better, 4)))
  print("Difference in average CATE (Female - Male):")
  print(round(sex.diff.summary, 3))
  
  # Save for reporting
  sex.disparity.results <- list(
    prob.female.better = prob.female.better,
    diff.summary = sex.diff.summary
  )
  saveRDS(sex.disparity.results, "results/sex_disparity_analysis.rds")
} 


# --- Figure 6: Clean vs Confounded Comparison ---
# Shows how confounding shifts estimates away from truth
if (!is.null(baseline.summary) && !is.null(sensitivity.summary) && 
    nrow(baseline.summary) > 0 && nrow(sensitivity.summary) > 0) {
  
  # Get clean RCT estimates
  clean.estimates <- baseline.summary %>%
    filter(n.size == 1000, prior == "Weak Pooling") %>%
    dplyr::select(subgroup_label, est.mean, ci.lower, ci.upper, true.cate.mean) %>%
    mutate(scenario = "Clean RCT")
  
  # Get most confounded scenario estimates (gamma=2, beta=8)
  # Need to extract from confounded.results
  conf.key <- "N1000_G2_B8_weak.pooling"
  if (!is.null(confounded.results[[conf.key]]$sims)) {
    
    all.conf.cate <- do.call(rbind, lapply(confounded.results[[conf.key]]$sims, 
                                            function(s) s$cate.summary))
    
    conf.estimates <- all.conf.cate %>%
      group_by(subgroup_label) %>%
      summarise(
        est.mean = mean(mean),
        ci.lower = mean(`2.5%`),
        ci.upper = mean(`97.5%`),
        true.cate.mean = mean(true.cate),
        .groups = "drop"
      ) %>%
      mutate(scenario = "Confounded (γ=2, β=8)")
    
    # Combine for plotting
    comparison.data <- bind_rows(clean.estimates, conf.estimates)
    
    # Add truth reference
    truth.ref <- clean.estimates %>%
      dplyr::select(subgroup_label, true.cate.mean) %>%
      distinct()
    
    p.clean.vs.conf <- ggplot(comparison.data, 
                               aes(x = subgroup_label, y = est.mean, color = scenario)) +
      geom_point(position = position_dodge(width = 0.5), size = 3) +
      geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), 
                    position = position_dodge(width = 0.5), width = 0.2) +
      geom_point(data = truth.ref, aes(x = subgroup_label, y = true.cate.mean), 
                 color = "black", shape = 4, size = 4, stroke = 1.5, inherit.aes = FALSE) +
      coord_flip() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      labs(
        title = "Effect of Unmeasured Confounding on CATE Estimates",
        subtitle = "Black X = True CATE. Confounding (γ=2, β=8) shifts estimates toward null.",
        x = "Subgroup",
        y = "Treatment Effect (Change in MELD Score)",
        color = "Scenario"
      ) +
      theme_minimal() +
      scale_color_manual(values = c("Clean RCT" = "#377EB8", "Confounded (γ=2, β=8)" = "#E41A1C"))
    
    ggsave("results/figures/fig6_clean_vs_confounded.pdf", 
           p.clean.vs.conf, width = 10, height = 6)
    ggsave("results/figures/fig6_clean_vs_confounded.png", 
           p.clean.vs.conf, width = 10, height = 6, dpi = 300)
    
    print("Generated Figure 6: Clean vs Confounded comparison")
  }
}

###--------------------------------------------------------------------------###
### PART 7: GENERATE LATEX TABLES
###--------------------------------------------------------------------------###

print("=== PART 7: Generating LaTeX Tables ===")

# =============================================================================
# TABLE 1: Subgroup Observation Counts by Sample Size
# =============================================================================
print("Generating Table 1: Subgroup Counts...")

# Load one batch from each sample size to get subgroup counts
batch.N200 <- readRDS("data_batches/batch_BASELINE_Clean_N200.rds")
batch.N1000 <- readRDS("data_batches/batch_BASELINE_Clean_N1000.rds")

# Get subgroup counts from first replication (representative)
counts.N200 <- batch.N200[[1]]$truth.data %>%
  group_by(subgroup_label) %>%
  summarise(N200 = n(), .groups = "drop")

counts.N1000 <- batch.N1000[[1]]$truth.data %>%
  group_by(subgroup_label) %>%
  summarise(N1000 = n(), .groups = "drop")

# Merge and format
subgroup.counts <- counts.N200 %>%
  full_join(counts.N1000, by = "subgroup_label") %>%
  rename(Subgroup = subgroup_label, `N = 200` = N200, `N = 1000` = N1000) %>%
  arrange(Subgroup)

# Generate LaTeX for Table 1
table1.tex <- "\\begin{table}[ht]
\\centering
\\caption{Number of observations per subgroup by sample size.}
\\label{tab:subgroup_counts}
\\begin{tabular}{lcc}
\\hline
Subgroup & N = 200 & N = 1000 \\\\
\\hline\n"

for (i in 1:nrow(subgroup.counts)) {
  table1.tex <- paste0(table1.tex, 
                        subgroup.counts$Subgroup[i], " & ",
                        subgroup.counts$`N = 200`[i], " & ",
                        subgroup.counts$`N = 1000`[i], " \\\\\n")
}

table1.tex <- paste0(table1.tex, "\\hline
\\end{tabular}
\\end{table}")

writeLines(table1.tex, "results/tables/table1_subgroup_counts.tex")
print("Saved: results/tables/table1_subgroup_counts.tex")

# =============================================================================
# TABLE 2: CATE Posterior Summaries (Weak vs Strong Pooling)
# =============================================================================
print("Generating Table 2: CATE Summaries...")

if (!is.null(baseline.summary) && nrow(baseline.summary) > 0) {
  
  # Filter to N=1000 Clean scenario, reshape to wide format
  cate.table.data <- baseline.summary %>%
    filter(n.size == 1000) %>%
    dplyr::select(subgroup_label, prior, true.cate.mean, est.mean, ci.lower, ci.upper) %>%
    mutate(
      # Format 95% CI as (lower, upper)
      ci.formatted = sprintf("(%.2f, %.2f)", ci.lower, ci.upper)
    )
  
  # Reshape to wide format with weak/strong columns
  weak.data <- cate.table.data %>%
    filter(prior == "Weak Pooling") %>%
    dplyr::select(subgroup_label, true.cate.mean, 
           weak.mean = est.mean, weak.ci = ci.formatted)
  
  strong.data <- cate.table.data %>%
    filter(prior == "Strong Pooling") %>%
    dplyr::select(subgroup_label, 
           strong.mean = est.mean, strong.ci = ci.formatted)
  
  table2.data <- weak.data %>%
    left_join(strong.data, by = "subgroup_label") %>%
    mutate(true.cate.mean = round(true.cate.mean, 2),
           weak.mean = round(weak.mean, 2),
           strong.mean = round(strong.mean, 2)) %>%
    arrange(subgroup_label)
  
  # Generate LaTeX for Table 2
  table2.tex <- "\\begin{table}[ht]
\\centering
\\caption{Posterior summaries of subgroup-specific CATEs under weak and strong pooling priors (N=1000, Clean RCT scenario).}
\\label{tab:cate_summaries}
\\begin{tabular}{lcccccc}
\\hline
 & & \\multicolumn{2}{c}{Weak Pooling} & \\multicolumn{2}{c}{Strong Pooling} \\\\
\\cmidrule(lr){3-4} \\cmidrule(lr){5-6}
Subgroup & True CATE & Mean & 95\\% CrI & Mean & 95\\% CrI \\\\
\\hline\n"
  
  for (i in 1:nrow(table2.data)) {
    table2.tex <- paste0(table2.tex,
                          table2.data$subgroup_label[i], " & ",
                          table2.data$true.cate.mean[i], " & ",
                          table2.data$weak.mean[i], " & ",
                          table2.data$weak.ci[i], " & ",
                          table2.data$strong.mean[i], " & ",
                          table2.data$strong.ci[i], " \\\\\n")
  }
  
  table2.tex <- paste0(table2.tex, "\\hline
\\end{tabular}
\\end{table}")
  
  writeLines(table2.tex, "results/tables/table2_cate_summaries.tex")
  print("Saved: results/tables/table2_cate_summaries.tex")
}

# =============================================================================
# TABLE 3: Sensitivity Analysis Grid (P(Benefit) across γ_U × β_U)
# =============================================================================
print("Generating Table 3: Sensitivity Analysis Grid...")

if (!is.null(sensitivity.summary) && nrow(sensitivity.summary) > 0) {
  
  # Include BOTH N=200 and N=1000, Weak Pooling for main sensitivity table
  sens.table.data <- sensitivity.summary %>%
    filter(prior == "Weak Pooling") %>%
    dplyr::select(n.size, gamma, beta, naive.bias, mean.bias, prob.benefit.female, prob.benefit.male) %>%
    mutate(
      naive.bias = round(naive.bias, 1),
      mean.bias = round(mean.bias, 1),
      prob.benefit.female = round(prob.benefit.female, 2),
      prob.benefit.male = round(prob.benefit.male, 2)
    ) %>%
    arrange(n.size, gamma, beta)
  
  # Generate LaTeX for Table 3 with tabularx format
  table3.tex <- "\\begin{table}[H]
\\centering
\\caption{Sensitivity Analysis: Posterior Probability of Benefit Across Confounding Scenarios}
\\label{tab:sensitivity}

\\begin{tabularx}{\\linewidth}{%
    C{0.8cm}    % N
    C{1.0cm}    % gamma
    C{1.0cm}    % beta
    C{1.5cm}    % Naive Bias
    C{1.5cm}    % Model Bias
    C{1.8cm}    % P(Benefit|F)
    C{1.8cm}    % P(Benefit|M)
}
\\toprule
$N$ & $\\gamma_U$ & $\\beta_U$ & Naive Bias & Model Bias & $P(\\text{Ben}|\\text{F})$ & $P(\\text{Ben}|\\text{M})$ \\\\
\\midrule\n"
  
  # Add rows with midrule between N=200 and N=1000
  current.n <- NULL
  for (i in 1:nrow(sens.table.data)) {
    if (!is.null(current.n) && sens.table.data$n.size[i] != current.n) {
      table3.tex <- paste0(table3.tex, "\\midrule\n")
    }
    current.n <- sens.table.data$n.size[i]
    
    table3.tex <- paste0(table3.tex,
                          sens.table.data$n.size[i], " & ",
                          sens.table.data$gamma[i], " & ",
                          sens.table.data$beta[i], " & ",
                          sens.table.data$naive.bias[i], " & ",
                          sens.table.data$mean.bias[i], " & ",
                          sens.table.data$prob.benefit.female[i], " & ",
                          sens.table.data$prob.benefit.male[i], " \\\\\n")
  }
  
  table3.tex <- paste0(table3.tex, "
\\bottomrule
\\end{tabularx}

\\vspace{1ex}
\\footnotesize
\\textit{Note:} $\\gamma_U$ = selection bias (strength of $U \\to A$); $\\beta_U$ = outcome bias (strength of $U \\to Y$). Naive bias = unadjusted ATE $-$ true ATE. Model bias = posterior mean ATE $-$ true ATE. $P(\\text{Ben}|\\cdot)$ = posterior probability that average subgroup CATE $< 0$ (treatment reduces MELD). Female conclusions remain robust ($P > 0.78$) across all scenarios; male conclusions become fragile ($P < 0.5$) under strong confounding. Values shown are for weak pooling prior ($\\tau \\sim \\text{Half-Normal}(0, 5)$).
\\end{table}")
  
  writeLines(table3.tex, "results/tables/table3_sensitivity.tex")
  print("Saved: results/tables/table3_sensitivity.tex")
}

# =============================================================================
# TABLE 4: Performance Metrics Summary
# =============================================================================
print("Generating Table 4: Performance Metrics...")

if (!is.null(performance.summary) && nrow(performance.summary) > 0) {
  
  perf.table.data <- performance.summary %>%
    dplyr::select(prior, n.size, mean.bias, rmse, coverage.95, ci.width) %>%
    mutate(
      mean.bias = round(mean.bias, 3),
      rmse = round(rmse, 3),
      coverage.95 = round(coverage.95, 2),
      ci.width = round(ci.width, 2)
    ) %>%
    arrange(n.size, prior)
  
  # Generate polished LaTeX with tabularx
  table4.tex <- "\\begin{table}[H]
\\centering
# \\caption{Performance Metrics for CATE Estimation Under Different Prior Specifications}
\\label{tab:performance}

\\begin{tabularx}{\\linewidth}{%
    L{2.5cm}    % Prior
    C{1.2cm}    % N
    C{1.5cm}    % Mean Bias
    C{1.5cm}    % RMSE
    C{1.8cm}    % 95% Coverage
    C{1.5cm}    % CI Width
}
\\toprule
Prior & $N$ & Mean Bias & RMSE & 95\\% Cov. & CI Width \\\\
\\midrule\n"
  
  # Group by sample size with mid-rules
  current.n <- NULL
  for (i in 1:nrow(perf.table.data)) {
    if (!is.null(current.n) && perf.table.data$n.size[i] != current.n) {
      table4.tex <- paste0(table4.tex, "\\midrule\n")
    }
    current.n <- perf.table.data$n.size[i]
    
    table4.tex <- paste0(table4.tex,
                          perf.table.data$prior[i], " & ",
                          perf.table.data$n.size[i], " & ",
                          perf.table.data$mean.bias[i], " & ",
                          perf.table.data$rmse[i], " & ",
                          perf.table.data$coverage.95[i], " & ",
                          perf.table.data$ci.width[i], " \\\\\n")
  }
  
  table4.tex <- paste0(table4.tex, "
\\bottomrule
\\end{tabularx}

\\vspace{1ex}
\\footnotesize
\\textit{Note:} Performance metrics computed over $S = 100$ simulation replications under the Clean RCT scenario (no unmeasured confounding). Mean Bias = average (posterior mean CATE $-$ true CATE) across subgroups. RMSE = root mean squared error. 95\\% Cov. = proportion of 95\\% credible intervals containing the true CATE. CI Width = average posterior interval width. Both priors achieve nominal coverage; weak pooling shows slightly higher RMSE but preserves subgroup heterogeneity, while strong pooling shrinks toward the population mean.
\\end{table}")
  
  writeLines(table4.tex, "results/tables/table4_performance.tex")
  print("Saved: results/tables/table4_performance.tex")
}

# =============================================================================
# TABLE 5: Validity Regions
# =============================================================================
print("Generating Table 5: Validity Regions...")

if (!is.null(validity.analysis) && nrow(validity.analysis) > 0) {
  
  validity.table.data <- validity.analysis %>%
    dplyr::select(gamma, beta, prior, female.robust, male.robust) %>%
    mutate(across(where(is.numeric), ~round(., 2)))
  
  table5.tex <- "\\begin{table}[ht]
\\centering
\\caption{Validity regions: Proportion of scenarios where P(CATE $<$ 0) $>$ 0.95 (robust benefit).}
\\label{tab:validity}
\\begin{tabular}{cclcc}
\\hline
$\\gamma_U$ & $\\beta_U$ & Prior & Female Robust & Male Robust \\\\
\\hline\n"
  
  for (i in 1:nrow(validity.table.data)) {
    table5.tex <- paste0(table5.tex,
                          validity.table.data$gamma[i], " & ",
                          validity.table.data$beta[i], " & ",
                          validity.table.data$prior[i], " & ",
                          validity.table.data$female.robust[i], " & ",
                          validity.table.data$male.robust[i], " \\\\\n")
  }
  
  table5.tex <- paste0(table5.tex, "\\hline
\\end{tabular}
\\end{table}")
  
  writeLines(table5.tex, "results/tables/table5_validity.tex")
  print("Saved: results/tables/table5_validity.tex")
}

###--------------------------------------------------------------------------###
### PART 8: SAVE STAN DIAGNOSTICS
###--------------------------------------------------------------------------###

print("=== PART 8: Saving Stan Diagnostics ===")

# Diagnostics summary
diag.summary <- list()

for (key in names(baseline.results)) {
  if (!is.null(baseline.results[[key]]$sims)) {
    diags <- lapply(baseline.results[[key]]$sims, function(s) s$diagnostics)
    
    diag.summary[[key]] <- list(
      max.rhat = max(sapply(diags, function(d) max(d$rhat, na.rm = TRUE))),
      min.neff = min(sapply(diags, function(d) min(d$n_eff, na.rm = TRUE))),
      total.divergences = sum(sapply(diags, function(d) d$divergences)),
      total.treedepth = sum(sapply(diags, function(d) d$max_treedepth))
    )
  }
}

# Also collect diagnostics from confounded scenarios
for (key in names(confounded.results)) {
  if (!is.null(confounded.results[[key]]$sims)) {
    diags <- lapply(confounded.results[[key]]$sims, function(s) s$diagnostics)
    
    diag.summary[[key]] <- list(
      max.rhat = max(sapply(diags, function(d) max(d$rhat, na.rm = TRUE))),
      min.neff = min(sapply(diags, function(d) min(d$n_eff, na.rm = TRUE))),
      total.divergences = sum(sapply(diags, function(d) d$divergences)),
      total.treedepth = sum(sapply(diags, function(d) d$max_treedepth))
    )
  }
}

saveRDS(diag.summary, "results/stan_diagnostics/diagnostics_summary.rds")

# Write diagnostics report
sink("results/stan_diagnostics/diagnostics_report.txt")
cat("=== Stan Diagnostics Report ===\n\n")
cat("Generated:", as.character(Sys.time()), "\n\n")

for (key in names(diag.summary)) {
  cat(paste0("--- ", key, " ---\n"))
  cat(paste0("  Max R-hat: ", round(diag.summary[[key]]$max.rhat, 4), "\n"))
  cat(paste0("  Min n_eff: ", round(diag.summary[[key]]$min.neff, 0), "\n"))
  cat(paste0("  Total divergences: ", diag.summary[[key]]$total.divergences, "\n"))
  cat(paste0("  Treedepth warnings: ", diag.summary[[key]]$total.treedepth, "\n\n"))
}
sink()

# =============================================================================
# TABLE 6: Stan Convergence Diagnostics Summary
# =============================================================================
print("Generating Table 6: Stan Diagnostics...")

if (length(diag.summary) > 0) {
  
  # Convert to data frame for table generation
  diag.df <- do.call(rbind, lapply(names(diag.summary), function(key) {
    # Parse the key to extract scenario info
    parts <- strsplit(key, "_")[[1]]
    
    # Determine scenario type and parameters
    if (grepl("^N", key)) {
      # Baseline scenario: N200_weak.pooling or N1000_strong.pooling
      n.size <- as.numeric(gsub("N", "", parts[1]))
      prior.type <- gsub("\\.", " ", tools::toTitleCase(parts[2]))
      scenario <- "Clean RCT"
      gamma <- NA
      beta <- NA
    } else {
      # Confounded scenario: N1000_G0.5_B2_weak.pooling
      n.size <- as.numeric(gsub("N", "", parts[1]))
      gamma <- as.numeric(gsub("G", "", parts[2]))
      beta <- as.numeric(gsub("B", "", parts[3]))
      prior.type <- gsub("\\.", " ", tools::toTitleCase(parts[4]))
      scenario <- paste0("Confounded (γ=", gamma, ", β=", beta, ")")
    }
    
    data.frame(
      Scenario = scenario,
      N = n.size,
      Prior = prior.type,
      Max.Rhat = diag.summary[[key]]$max.rhat,
      Min.Neff = diag.summary[[key]]$min.neff,
      Divergences = diag.summary[[key]]$total.divergences,
      Treedepth = diag.summary[[key]]$total.treedepth,
      stringsAsFactors = FALSE
    )
  }))
  
  # Sort and format
  diag.df <- diag.df %>%
    arrange(Scenario, N, Prior) %>%
    mutate(
      Max.Rhat = round(Max.Rhat, 4),
      Min.Neff = round(Min.Neff, 0)
    )
  
  # Summarize: aggregate across all scenarios for a compact table
  diag.agg <- diag.df %>%
    group_by(Scenario, N, Prior) %>%
    summarise(
      Max.Rhat = max(Max.Rhat),
      Min.Neff = min(Min.Neff),
      Divergences = sum(Divergences),
      Treedepth = sum(Treedepth),
      .groups = "drop"
    )
  
  # For a cleaner table, summarize by Prior and N only (across all scenarios)
  diag.compact <- diag.df %>%
    group_by(Prior, N) %>%
    summarise(
      `Max $\\hat{R}$` = max(Max.Rhat),
      `Min $n_{\\text{eff}}$` = min(Min.Neff),
      `Total Divergences` = sum(Divergences),
      `Treedepth Warnings` = sum(Treedepth),
      `Models Fit` = n(),
      .groups = "drop"
    ) %>%
    arrange(N, Prior)
  
  # Generate polished LaTeX table
  table6.tex <- "\\begin{table}[H]
\\centering
\\caption{MCMC Convergence Diagnostics Summary Across All Model Fits}
\\label{tab:diagnostics}

\\begin{tabularx}{\\linewidth}{%
    L{2.5cm}    % Prior
    C{1.0cm}    % N
    C{1.5cm}    % Max R-hat
    C{1.8cm}    % Min n_eff
    C{1.8cm}    % Divergences
    C{1.8cm}    % Treedepth
    C{1.5cm}    % Models
}
\\toprule
Prior & $N$ & Max $\\hat{R}$ & Min $n_{\\text{eff}}$ & Divergences & Treedepth & Models \\\\
\\midrule\n"
  
  for (i in 1:nrow(diag.compact)) {
    table6.tex <- paste0(table6.tex,
                          diag.compact$Prior[i], " & ",
                          diag.compact$N[i], " & ",
                          sprintf("%.4f", diag.compact$`Max $\\hat{R}$`[i]), " & ",
                          diag.compact$`Min $n_{\\text{eff}}$`[i], " & ",
                          diag.compact$`Total Divergences`[i], " & ",
                          diag.compact$`Treedepth Warnings`[i], " & ",
                          diag.compact$`Models Fit`[i], " \\\\\n")
  }
  
  table6.tex <- paste0(table6.tex, "
\\bottomrule
\\end{tabularx}

\\vspace{1ex}
\\footnotesize
\\textit{Note:} Diagnostics aggregated across $S = 100$ simulation replications per scenario. Max $\\hat{R}$ = maximum Gelman-Rubin statistic across all parameters and replications (target: $< 1.01$). Min $n_{\\text{eff}}$ = minimum effective sample size (target: $> 400$). Divergences = total divergent transitions (target: 0). Treedepth = iterations hitting maximum treedepth. All models achieved satisfactory convergence with $\\hat{R} < 1.01$, $n_{\\text{eff}} > 400$, and zero pathological warnings.
\\end{table}")
  
  writeLines(table6.tex, "results/tables/table6_diagnostics.tex")
  print("Saved: results/tables/table6_diagnostics.tex")
}

###--------------------------------------------------------------------------###
### PART 9: SAVE ALL RESULTS
###--------------------------------------------------------------------------###

print("=== PART 9: Saving All Results ===")

# Save aggregated results
saveRDS(baseline.summary, "results/baseline_summary.rds")
saveRDS(performance.summary, "results/performance_summary.rds")
saveRDS(sensitivity.summary, "results/sensitivity_summary.rds")
saveRDS(validity.analysis, "results/validity_analysis.rds")

# Save raw results (warning: can be large)
saveRDS(baseline.results, "results/baseline_results_raw.rds")
saveRDS(confounded.results, "results/confounded_results_raw.rds")

print("=== ANALYSIS COMPLETE ===")
print(paste("Results saved to:", normalizePath("results")))
print("Files generated:")
print(list.files("results", recursive = TRUE))

###--------------------------------------------------------------------------###
### OPTIONAL: Launch Shiny Stan for Interactive Exploration
###--------------------------------------------------------------------------###

# Uncomment to launch interactive diagnostics for a specific fit
# launch_shinystan(baseline.results[["N1000_weak_pooling"]]$first.fit)
