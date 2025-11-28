###--------------------------------------------------------------------------###
### Project: Quantifying Uncertainty in Causal Effects
### Step 1: Model Validation (The "Gold Standard" Check)
### 
### Objective: 
### 1. Generate a "Clean" dataset (Large N, Randomized Trial).
### 2. Fit the Partial Pooling Stan model.
### 3. Verify that the model recovers the TRUE CATEs defined in the simulation.
###--------------------------------------------
rm(list=ls())
library(tidyverse)
library(latex2exp)
library(rstan)
library(shinystan)
library(bayesplot)
setwd("/Users/adamkurth/Documents/RStudio/research/bayes-priors-sens")

# load data generating function 
source("00_data.r")

# optimize for parallel processing
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

###--------------------------------------------------------------------------###
### 1. Generate "Clean" Data
###--------------------------------------------------------------------------###
# We use N=2000 and scenario="clean" (Randomized) to ensure the model 
# has the best possible chance of succeeding. If it fails here, it's broken.

sim <- simulate.liver.data(N = 2000, scenario = "clean", seed = 123)
d.clean <- sim$d
stan.data <- sim$stan.data


sim <- simulate.liver.data(N = 2000, scenario = "confounded", seed = 123)
d.clean <- sim$d
stan.data <- sim$stan.data

# data check 
print("--- Data Check--- ")
print(paste("Sample Size", stan.data$N))
print("Subgroup Counts:")
print(stan.data$n_v)

###--------------------------------------------------------------------------###
### 2. Configure Priors for "Weak Pooling" (Baseline)
###--------------------------------------------------------------------------###
# For the validation run we use "weak pooling"
# We set WIDE prior on tau (heterogenity)
# Tells model: subgroups might be very different, dont shrink too much 

stan.data$prior_tau_sigma <- 5      # Wide half-normal prior on heterogentity
stan.data$prior_mean_mu <- 0        # Weakly informative prior on ATE 

###--------------------------------------------------------------------------###
### 3. Run Stan Model
###--------------------------------------------------------------------------###
print("--- Compiling and Running Stan Model ---")

# Ensure the .stan file is in the working directory
fit <- stan(
    file = "partial_pool_continuous.stan",
    data=stan.data, 
    iter = 2000,
    warmup = 1000, 
    chains = 2, # use 2 for speed in testing 
    seed = 2025
)

# rstan::check_all_diagnostics(fit)
rstan::traceplot(fit, pars = c("mu", "tau"))
# launch_shinystan(fit)


###--------------------------------------------------------------------------###
### DIAGNOSTIC 1: Convergence (Traceplots)
###--------------------------------------------------------------------------###
# We look for "Fuzzy Caterpillars" (good mixing)
color_scheme_set("blue")
p.trace <- mcmc_trace(fit, pars = c("mu", "tau", "sigma")) +
  ggtitle("Diagnostic 1: Chain Convergence (Traceplots)")
print(p.trace)

###--------------------------------------------------------------------------###
### DIAGNOSTIC 2: Model Fit (Posterior Predictive Check)
###--------------------------------------------------------------------------###
# Does the model generate data that looks like our actual MELD scores?
# y = Real Data, y_rep = Model Simulated Data
y.rep <- as.matrix(fit, pars = "y.rep")
y.obs <- stan.data$Y

# Compare density of observed vs 50 simulated datasets
p.ppc <- ppc_dens_overlay(y.obs, y.rep[1:50, ]) +
  ggtitle("Diagnostic 2: Posterior Predictive Check",
          subtitle = "Black = Real MELD Data, Blue Lines = Model Simulations")
print(p.ppc)

###--------------------------------------------------------------------------###
### VISUAL 3: The "Sex Disparity" (Research Question Check)
###--------------------------------------------------------------------------###
# We want to clearly show that Females get MORE benefit (more negative effect)
# than Males. We extract the CATEs and plot their full densities.

posterior_cates <- as.matrix(fit, pars = "CATE")
# Map column names (CATE[1], CATE[2]...) to Subgroup Labels
colnames(posterior_cates) <- levels(as.factor(d.clean$subgroup_label))

# Convert to Long Format for Plotting
d_post_long <- as.data.frame(posterior_cates) %>%
  pivot_longer(everything(), names_to = "Subgroup", values_to = "CATE_Estimate") %>%
  mutate(Sex = ifelse(grepl("^F", Subgroup), "Female", "Male"))

# Plot Densities
p_disparity <- ggplot(d_post_long, aes(x = CATE_Estimate, fill = Sex)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = -3, linetype="dashed", color="black") +
  annotate("text", x=-2.8, y=0.5, label="Avg Effect (-3)", angle=90) +
  labs(
    title = "Visual 3: The Sex Disparity in Treatment Effect",
    subtitle = "Posterior distributions show Females (Red) have a stronger response than Males (Blue).\nThis validates that the model recovered the heterogeneity.",
    x = "Estimated Change in MELD Score (Treatment Effect)",
    y = "Posterior Density"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8"))

print(p_disparity)

###--------------------------------------------------------------------------###
### VISUAL 4: Robust Forest Plot (Subgroup vs Truth)
###--------------------------------------------------------------------------###
# Extract Summary
res <- summary(fit, pars = "CATE")$summary %>% as.data.frame()
subgroup.map <- d.clean %>% 
  select(subgroup_id, subgroup_label, true.cate, sex, ascites) %>%
  distinct() %>%
  arrange(subgroup_id) %>%
  mutate(Sex = ifelse(sex==1, "Female", "Male"))

results <- cbind(subgroup.map, res)

# Forest Plot
p_forest <- ggplot(results, aes(x = subgroup_label, y = mean, color = Sex)) +
  # 95% CI (Thin line)
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1, linewidth = 0.5) +
  # 50% CI (Thick line - where the mass is)
  geom_errorbar(aes(ymin = `25%`, ymax = `75%`), width = 0, linewidth = 1.5) +
  # The Estimate
  geom_point(size = 3) +
  # The Truth (Red X)
  geom_point(aes(y = true.cate), color = "black", shape = 4, size = 4, stroke=1.5) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype="solid", color="grey") +
  labs(
    title = "Visual 4: CATE Forest Plot (Estimates vs Truth)",
    subtitle = "Black X = True Simulation Parameter. Dots = Model Estimate.\nThick bar = 50% CI, Thin bar = 95% CI.",
    y = "Treatment Effect (MELD Reduction)",
    x = ""
  ) +
  theme_minimal() + 
  scale_color_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8"))

print(p_forest)


###--------------------------------------------------------------------------###
### 4. Compare Estimates vs Truth
###--------------------------------------------------------------------------###

# extract CATE estimates (theta)
# summary includes Mean, SD, Quantiles (2.5%, 97.5%)
res <- summary(fit, pars = "CATE")$summary %>% as.data.frame()

# add subgroup labels to results 
# map to numeric indices (1,...8) back to "F_Ascites_Old" etc.
subgroup.map <- d.clean %>% 
    select(subgroup_id, subgroup_label, true.cate) %>%
    distinct() %>%
    arrange(subgroup_id)

results <- cbind(subgroup.map, res)

print("--- Validation Results ---")
print(results %>% select(subgroup_label, true.cate, mean, `2.5%`, `97.5%`))

# VISUAL CHECK
# Plot True CATE (Red X) vs Estimated 95% CI (Black bars)
p <- ggplot(results, aes(x = subgroup_label, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  geom_point(aes(y = true.cate), color = "red", shape = 4, size = 3) +
  geom_hline(yintercept = -3, linetype = "dashed", color = "blue", alpha=0.5) + # Global ATE roughly
  coord_flip() +
  labs(
    title = "Validation: Recovering True CATEs (MELD 3.0 Simulation)",
    subtitle = "Red X = Truth, Black Dot = Estimate. (Scenario: Clean RCT, Weak Pooling)\nFemales (F_...) should be around -5. Males (M_...) around -3.",
    y = "Change in MELD Score (Treatment Effect)",
    x = "Subgroup"
  ) +
  theme_minimal()

print(p)
