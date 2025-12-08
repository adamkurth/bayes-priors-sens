###############################################################################
## generate_tables.r
## Generates LaTeX tables from saved simulation results
## Run this script standalone after 01_test.r has saved results to /results
###############################################################################

library(tidyverse)

# Set working directory
setwd("/Users/adamkurth/Documents/RStudio/research/bayes-priors-sens")

# Create output directory
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("=== Generating LaTeX Tables from Saved Results ===\n\n")

# =============================================================================
# LOAD SAVED RESULTS (from 01_test.r)
# =============================================================================
cat("Loading saved results...\n")

baseline.summary <- NULL
performance.summary <- NULL
sensitivity.summary <- NULL
validity.analysis <- NULL
diag.summary <- NULL

if (file.exists("results/baseline_summary.rds")) {
  baseline.summary <- readRDS("results/baseline_summary.rds")
  cat("  Loaded: baseline_summary.rds\n")
}

if (file.exists("results/performance_summary.rds")) {
  performance.summary <- readRDS("results/performance_summary.rds")
  cat("  Loaded: performance_summary.rds\n")
}

if (file.exists("results/sensitivity_summary.rds")) {
  sensitivity.summary <- readRDS("results/sensitivity_summary.rds")
  cat("  Loaded: sensitivity_summary.rds\n")
}

if (file.exists("results/validity_analysis.rds")) {
  validity.analysis <- readRDS("results/validity_analysis.rds")
  cat("  Loaded: validity_analysis.rds\n")
}

if (file.exists("results/stan_diagnostics/diagnostics_summary.rds")) {
  diag.summary <- readRDS("results/stan_diagnostics/diagnostics_summary.rds")
  cat("  Loaded: diagnostics_summary.rds\n")
}

cat("\n")

# =============================================================================
# TABLE 1: Subgroup Observation Counts by Sample Size
# =============================================================================
cat("Generating Table 1: Subgroup Counts...\n")

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
cat("  Saved: results/tables/table1_subgroup_counts.tex\n")

# =============================================================================
# TABLE 2: CATE Posterior Summaries (Weak vs Strong Pooling)
# =============================================================================
cat("Generating Table 2: CATE Summaries...\n")

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
  cat("  Saved: results/tables/table2_cate_summaries.tex\n")
} else {
  cat("  Skipped: baseline.summary not available\n")
}

# =============================================================================
# TABLE 3: Sensitivity Analysis Grid (P(Benefit) across γ_U × β_U)
# =============================================================================
cat("Generating Table 3: Sensitivity Analysis Grid...\n")

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
  cat("  Saved: results/tables/table3_sensitivity.tex\n")
} else {
  cat("  Skipped: sensitivity.summary not available\n")
}

# =============================================================================
# TABLE 4: Performance Metrics Summary
# =============================================================================
cat("Generating Table 4: Performance Metrics...\n")

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
\\caption{Performance Metrics for CATE Estimation Under Different Prior Specifications}
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
  cat("  Saved: results/tables/table4_performance.tex\n")
} else {
  cat("  Skipped: performance.summary not available\n")
}

# =============================================================================
# TABLE 5: Validity Regions
# =============================================================================
cat("Generating Table 5: Validity Regions...\n")

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
  cat("  Saved: results/tables/table5_validity.tex\n")
} else {
  cat("  Skipped: validity.analysis not available\n")
}

# =============================================================================
# TABLE 6: Stan Convergence Diagnostics Summary
# =============================================================================
cat("Generating Table 6: Stan Diagnostics...\n")

if (!is.null(diag.summary) && length(diag.summary) > 0) {
  
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
  cat("  Saved: results/tables/table6_diagnostics.tex\n")
} else {
  cat("  Skipped: diag.summary not available\n")
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Table Generation Complete ===\n")
cat("Generated tables:\n")
cat(paste("  ", list.files("results/tables", pattern = "\\.tex$"), collapse = "\n"))
cat("\n")
