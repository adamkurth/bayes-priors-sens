###--------------------------------------------------------------------------###
### Generate LaTeX Content for Write-up
###
### This script loads saved results and generates publication-ready
### tables and figures for the Results section.
###--------------------------------------------------------------------------###

rm(list = ls())
setwd("/Users/adamkurth/Documents/RStudio/research/bayes-priors-sens")

library(tidyverse)
library(xtable)
library(knitr)

###--------------------------------------------------------------------------###
### Load Results
###--------------------------------------------------------------------------###

if (!file.exists("results/baseline_summary.rds")) {
  stop("Results not found. Run 01_test.r first.")
}

baseline_summary <- readRDS("results/baseline_summary.rds")
performance_summary <- readRDS("results/performance_summary.rds")
sensitivity_summary <- readRDS("results/sensitivity_summary.rds")
validity_analysis <- readRDS("results/validity_analysis.rds")

###--------------------------------------------------------------------------###
### Table 1: Posterior Summaries - Weak vs Strong Pooling
###--------------------------------------------------------------------------###

# Format for LaTeX
table1_data <- baseline_summary %>%
  filter(n_size == 1000) %>%
  select(Subgroup = subgroup_label, 
         Prior = prior,
         Truth = true_cate,
         Mean = est_mean,
         Lower = ci_lower,
         Upper = ci_upper) %>%
  mutate(
    CrI = sprintf("[%.2f, %.2f]", Lower, Upper),
    Bias = Mean - Truth
  ) %>%
  select(Subgroup, Prior, Truth, Mean, CrI, Bias) %>%
  pivot_wider(
    names_from = Prior,
    values_from = c(Mean, CrI, Bias),
    names_glue = "{Prior}_{.value}"
  )

cat("\n=== Table 1: CATE Posterior Summaries ===\n\n")
print(kable(table1_data, format = "latex", booktabs = TRUE, digits = 2))

###--------------------------------------------------------------------------###
### Table 2: Performance Metrics
###--------------------------------------------------------------------------###

table2_data <- performance_summary %>%
  mutate(
    `Sample Size` = paste0("N = ", n_size),
    Prior = prior,
    `Mean Bias` = round(mean_bias, 3),
    RMSE = round(rmse, 3),
    `Coverage` = paste0(round(coverage_95 * 100, 1), "\\%"),
    `CI Width` = round(ci_width, 2)
  ) %>%
  select(`Sample Size`, Prior, `Mean Bias`, RMSE, Coverage, `CI Width`)

cat("\n=== Table 2: Performance Comparison ===\n\n")
print(xtable(table2_data,
             caption = "Performance metrics for CATE estimation. Coverage refers to the proportion of 95\\% credible intervals that contain the true CATE.",
             label = "tab:performance"),
      sanitize.text.function = identity,
      include.rownames = FALSE)

###--------------------------------------------------------------------------###
### Table 3: Sensitivity Analysis Grid
###--------------------------------------------------------------------------###

table3_data <- sensitivity_summary %>%
  filter(prior == "Weak Pooling") %>%
  arrange(n_size, gamma, beta) %>%
  mutate(
    N = n_size,
    `$\\gamma_U$` = gamma,
    `$\\beta_U$` = beta,
    `Naive Bias` = round(naive_bias, 2),
    `Model Bias` = round(mean_bias, 2),
    `P(Benefit|F)` = round(prob_benefit_female, 2),
    `P(Benefit|M)` = round(prob_benefit_male, 2)
  ) %>%
  select(N, `$\\gamma_U$`, `$\\beta_U$`, `Naive Bias`, `Model Bias`, 
         `P(Benefit|F)`, `P(Benefit|M)`)

cat("\n=== Table 3: Sensitivity Analysis ===\n\n")
print(xtable(table3_data,
             caption = "Sensitivity analysis across confounding scenarios. $\\gamma_U$ represents selection bias (U$\\to$A) and $\\beta_U$ represents outcome bias (U$\\to$Y). P(Benefit) is the posterior probability that the CATE is negative (beneficial).",
             label = "tab:sensitivity"),
      sanitize.text.function = identity,
      include.rownames = FALSE)

###--------------------------------------------------------------------------###
### Table 4: Validity Regions
###--------------------------------------------------------------------------###

# Identify validity regions where conclusions are robust
validity_table <- sensitivity_summary %>%
  filter(prior == "Weak Pooling", n_size == 1000) %>%
  mutate(
    Female_Robust = ifelse(prob_benefit_female > 0.95, "\\checkmark", ""),
    Male_Robust = ifelse(prob_benefit_male > 0.95, "\\checkmark", ""),
    Conclusion_Robust = ifelse(prob_benefit_female > 0.95 & prob_benefit_male > 0.95, 
                               "Both", 
                               ifelse(prob_benefit_female > 0.95, "Female only",
                                      ifelse(prob_benefit_male > 0.95, "Male only", "Neither")))
  ) %>%
  select(`$\\gamma_U$` = gamma, `$\\beta_U$` = beta, 
         `Female Robust` = Female_Robust, 
         `Male Robust` = Male_Robust,
         Conclusion = Conclusion_Robust)

cat("\n=== Table 4: Validity Regions ===\n\n")
print(xtable(validity_table,
             caption = "Validity regions where posterior probability of benefit exceeds 0.95. Checkmark indicates the subgroup conclusion is robust to the specified level of unmeasured confounding.",
             label = "tab:validity"),
      sanitize.text.function = identity,
      include.rownames = FALSE)

###--------------------------------------------------------------------------###
### Generate Summary Statistics for Text
###--------------------------------------------------------------------------###

cat("\n=== Key Statistics for Write-up ===\n\n")

# Average CATE estimates
if (!is.null(baseline_summary) && nrow(baseline_summary) > 0) {
  
  female_cates <- baseline_summary %>%
    filter(n_size == 1000, prior == "Weak Pooling", sex == 1)
  
  male_cates <- baseline_summary %>%
    filter(n_size == 1000, prior == "Weak Pooling", sex == 0)
  
  cat(sprintf("Female subgroups - Mean CATE: %.2f (True: %.2f)\n",
              mean(female_cates$est_mean), mean(female_cates$true_cate)))
  
  cat(sprintf("Male subgroups - Mean CATE: %.2f (True: %.2f)\n",
              mean(male_cates$est_mean), mean(male_cates$true_cate)))
}

# Sensitivity findings
if (!is.null(sensitivity_summary) && nrow(sensitivity_summary) > 0) {
  
  # Where does benefit become uncertain?
  uncertain_regions <- sensitivity_summary %>%
    filter(prior == "Weak Pooling", n_size == 1000) %>%
    filter(prob_benefit_female < 0.95 | prob_benefit_male < 0.95)
  
  if (nrow(uncertain_regions) > 0) {
    cat("\nConfounding scenarios where benefit becomes uncertain:\n")
    print(uncertain_regions %>% select(gamma, beta, prob_benefit_female, prob_benefit_male))
  }
  
  # Maximum bias from confounding
  max_naive_bias <- max(abs(sensitivity_summary$naive_bias), na.rm = TRUE)
  cat(sprintf("\nMaximum naive bias from confounding: %.2f MELD points\n", max_naive_bias))
}

###--------------------------------------------------------------------------###
### LaTeX Code Snippets for Results Section
###--------------------------------------------------------------------------###

cat("\n\n=== LaTeX Snippets for Results Section ===\n\n")

# Posterior summaries paragraph
cat("% Posterior Summaries Paragraph\n")
cat("\\paragraph{Posterior Summaries}\n")
cat("Table~\\ref{tab:cate_summaries} presents the posterior means and 95\\% credible intervals\n")
cat("for subgroup-specific CATEs under both weak and strong pooling priors.\n")
cat("Under weak pooling ($\\tau \\sim \\text{Half-Normal}(0, 5)$), the model\n")
cat("allows substantial heterogeneity across subgroups, while strong pooling\n")
cat("($\\tau \\sim \\text{Half-Normal}(0, 0.5)$) shrinks estimates toward the global mean.\n\n")

# Sensitivity analysis paragraph  
cat("% Sensitivity Analysis Paragraph\n")
cat("\\paragraph{Sensitivity to Unmeasured Confounding}\n")
cat("Table~\\ref{tab:sensitivity} displays how posterior conclusions vary across\n")
cat("a grid of confounding scenarios defined by selection bias $\\gamma_U$\n")
cat("and outcome bias $\\beta_U$. The validity region (Table~\\ref{tab:validity})\n")
cat("identifies combinations where the probability of benefit remains above 0.95.\n\n")

print("=== LaTeX Generation Complete ===")
