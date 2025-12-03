###--------------------------------------------------------------------------###
### Project: Quantifying Uncertainty in Causal Effects
### Step 0: Data Simulation
### 
### CONTEXT: MELD 3.0 vs. MELDNa (Sex Disparity in Liver Transplant)
###
### MOTIVATION: 
###   The MELDNa score historically underestimated disease severity in women 
###   because it relies on Creatinine. Women have lower muscle mass (Sarcopenia),
###   leading to lower Creatinine, lower MELD scores, and higher waitlist mortality 
###   (Kim et al., 2021). MELD 3.0 adds a +1.3 point boost for females.
###
###   This simulation tests a hypothetical "Intensive Nutrition Protocol" (Treatment)
###   designed to counteract Sarcopenia.
###
### OBJECTIVE:
###   1. Part 1 (Priors): Can informative priors on Sex-Subgroups (borrowing 
###      strength from Kim et al.) better estimate the benefit of Nutrition?
###   2. Part 2 (Sensitivity): If we cannot measure Sarcopenia (Confounder U), 
###      how wrong are our conclusions?
###--------------------------------------------------------------------------###
library(tidyverse)
library(latex2exp)

# Set seed for reproducibility
set.seed(2025)
# setwd("/Users/adamkurth/Documents/RStudio/research/bayes-priors-sens")


# Define sample size (Registry data is usually large, but we use N=1000 for demo)

###--------------------------------------------------------------------------###
#' Simulate Liver Data with Tunable Confounding
#'
#' @param N Sample Size
#' @param scenario "confounded" or "clean"
#' @param gamma_u Effect of U on Treatment (Log-Odds). Higher = Stronger Selection Bias.
#' @param beta_u Effect of U on Outcome (MELD Points). Higher = Stronger Outcome Bias.

simulate.liver.data <- function(N = 1000, 
                                scenario = "confounded", 
                                gamma.u = 1.5, # Default: Strong Selection
                                beta.u = 4,    # Default: Strong Outcome Bias
                                seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  ###--------------------------------------------------------------------------###
  ### 1. Define Sample and Covariates ("The Patient Profile")
  ###--------------------------------------------------------------------------###

  # Create 3 covariates: SEX, ASCITES, AGE
  # These define out "subgroups" for part 1.

  # A. SEX (The central variable of interest)
  # 0 = Male, 1 = Female
  # We assume a roughly balanced waitlist.
  sex <- rbinom(N, size = 1, prob = 0.45) 

  # B. ASCITES (Covariate for Subgrouping)
  # 0 = Absent/Controlled, 1 = Moderate/Severe Refractory
  # Ascites is a strong driver of decompensation.
  # We assume females might have slightly higher rates in this autoimmune-heavy sample.
  prob.ascites <- 0.3 + (0.1 * sex) 
  ascites <- rbinom(N, size = 1, prob = prob.ascites)

  # C. AGE GROUP (Covariate for Subgrouping)
  # 0 = <60 years, 1 = >=60 years
  # Older patients are often frailer.
  age.high <- rbinom(N, size = 1, prob = 0.4)

  # dataframe for observed covariates 
  d <- data.frame(
      id = 1:N, 
      sex = sex, 
      ascites = ascites, 
      age.high = age.high
  )

  ###--------------------------------------------------------------------------###
  ### 2. Define Subgroups (Strata) for Partial Pooling (Part 1)
  ###--------------------------------------------------------------------------###

  # We stratify patients into unique groups based on Sex, Ascites, and Age.
  # This creates 8 subgroups (2x2x2).
  # Variable 'V' in Oganisian & Roy logic.

  d <- d %>% 
    mutate(subgroup_label = paste0(
      ifelse(sex==1, "F", "M"), "_",
      ifelse(ascites==1, "Ascites", "NoAsc"), "_",
      ifelse(age.high==1, "Old", "Young")
    )) %>%
      mutate(subgroup_id = as.numeric(as.factor(subgroup_label)))

  # CHECK SPARSITY:
  # In Part 1, we argue that estimates for sparse subgroups (e.g., Female+Ascites+Old)
  # are unstable using standard frequentist methods. 
  # We will use Partial Pooling to "borrow strength" from other groups.
  # Skeptical Prior = Borrow strength from Men (assume F=M).
  # Informative Prior = Borrow strength from Women (assume F diff from M).

  # print("Sample Size per Subgroup:")
  # print(table(d$subgroup_label))

  ###--------------------------------------------------------------------------###
  ### 3. Simulate Unmeasured Confounding: SARCOPENIA (Part 2)
  ###--------------------------------------------------------------------------###

  # Variable 'U': Sarcopenia / Frailty.
  # This is the crux of the Kim et al. argument: Women have less muscle mass.
  # In registry data (UNOS), we rarely have psoas muscle area scans. It is UNMEASURED.

  # We generate U such that it is higher in Females (reflecting higher risk of muscle-based bias).
  # U ~ Normal(0, 1) shifted by sex.
  d$U <- rnorm(N, mean = 0, sd = 1) + (0.5 * d$sex)

  ###--------------------------------------------------------------------------###
  ### 4. Simulate Treatment Assignment: NUTRITION PROTOCOL (A)
  ###--------------------------------------------------------------------------###

  # Treatment A = 1 (Intensive Pre-Transplant Nutrition) vs A = 0 (Standard Care).
  # Hypothesis: Nutrition helps stabilize muscle mass and lower MELD scores.

  # PROPENSITY SCORE (Confounding Mechanism):
  # Who gets prescribed intensive nutrition?
  # 1. Patients with Ascites (visible malnutrition).
  # 2. Patients who look Frail/Sarcopenic (High U).
  # This creates CONFOUNDING BY INDICATION: The sickest people get the treatment.

  if(scenario == "confounded") {
      # Confounding: Frail (High U) and Ascites get treatment
      z.score <- -0.5 + (1.2 * d$ascites) + (1.5 * d$U) - (0.5 * d$age.high)
      prob.A <- plogis(z.score)
      d$A <- rbinom(N, size = 1, prob = prob.A)
  } else {
      # Randomized (Clean)
      z.score <- -0.5 + (1.2 * d$ascites) - (0.5 * d$age.high)
      prob.A <- plogis(z.score)
      d$A <- rbinom(N, size = 1, prob = prob.A)
  }

  # Check overlap (Positivity)
  # print("Probability of Treatment by Assignment:")
  # print(tapply(prob.A, d$A, summary))

  ###--------------------------------------------------------------------------###
  ### 5. Simulate Outcome: MELD 3.0 SCORE at 90 Days (Y)
  ###--------------------------------------------------------------------------###

  # Outcome Y: MELD 3.0 Score. 
  # NOTE: Range is 6 to 40. HIGHER SCORE = SICKER / WORSE.
  # A beneficial treatment should yield a NEGATIVE effect (Lower Score).

  # A. BASELINE MELD (The Counterfactual Y0)
  # - Base Score: 15
  # - Sex: Kim et al. show Females are naturally "sicker" than MELDNa captured. 
  #   In the "Truth" (MELD 3.0), we add 1.3 points for females.
  # - Ascites: Adds significant mortality risk (+8 points).
  # - Sarcopenia (U): Muscle wasting independently drives liver failure (+4 points).
  baseline.meld <- 15 + (1.3 * d$sex) + (8 * d$ascites) + (4 * d$U) 

  # B. HETEROGENEOUS TREATMENT EFFECTS (True CATEs)
  # This is what Part 1 tries to estimate.
  # - General Effect: Nutrition lowers MELD by 3 points (-3).
  # - INTERACTION: Females benefit MORE from nutrition because it targets their 
  #   specific deficit (sarcopenia). Effect is (-3) + (-2) = -5 points total.
  # - INTERACTION: Ascites patients are refractory, benefit LESS (+1).
  true.cate <- -3 + (-2 * d$sex) + (1 * d$ascites)

  # C. GENERATE OBSERVED Y
  # Y = Baseline + (TreatmentEffect * A) + Noise
  sigma.y <- 3 
  y.latent <- baseline.meld + (true.cate * d$A) + rnorm(N, mean=0, sd=sigma.y)  

  # Clamp to realistic MELD range (6 to 40) 
  d$Y.obs <- round(y.latent) # MELD is an integer score
  d$Y.obs <- pmax(6, pmin(40, d$Y.obs))

  # Store "True" CATE for validation later
  d$true.cate <- true.cate

  ###--------------------------------------------------------------------------###
  ### 6. Format for Stan (Bayesian Analysis)
  ###--------------------------------------------------------------------------###

  # W Matrix: Fixed Covariates we control for (Sex, Ascites, Age).
  # Corresponds to 'W' in Oganisian/Roy paper.
  W_matrix <- model.matrix(~ 1 + sex + ascites + age.high, data = d)

  # V Matrix: Subgroup Indicators.
  # Used for the Hierarchical/Partial Pooling priors.
  V_matrix <- model.matrix(~ 0 + as.factor(subgroup_id), data = d)

  # Package into list for Stan
  stan_data_list <- list(
    N = N,
    Y = d$Y.obs,                  # Outcome: MELD 3.0 Score
    A = d$A,                      # Treatment: Nutrition
    W = W_matrix,                 # Observed Confounders
    V = V_matrix,                 # Subgroup Definitions
    Pw = ncol(W_matrix),          # Number of confounders
    Pv = ncol(V_matrix),          # Number of subgroups
    n_v = as.numeric(table(d$subgroup_id)), # Count per subgroup
    # Helper for indexing in Stan (cumulative sum of counts)
    ind = c(0, cumsum(as.numeric(table(d$subgroup_id)))) 
  )

  
  # Allow passing priors dynamically later (Defaults for now)
  stan_data_list$prior_tau_sigma <- 5 
  stan_data_list$prior_mean_mu <- 0

  return(list(d=d, stan.data = stan_data_list))
} 

# END OF SCRIPT 

# 
# # Sanity Check Usage
# sim.output <- simulate.liver.data(N = 1000, scenario = "confounded", seed = 2025)
# d <- sim.output$d
# stan_data_list <- sim.output$stan.data
# 
# print("--- Data Gen Complete ---")
# print(head(d))

# print("Mean True CATE by Sex (expected Females to be more negative): ")
# print(aggregate(true_cate ~ sex, data=sim.output$d, mean))

###--------------------------------------------------------------------------###
### 7. Save Data for Analysis Steps
###--------------------------------------------------------------------------###

# 1. Save the FULL dataset with Truth (U). 
# We use this to check if our sensitivity analysis covered the true bias.
# saveRDS(d, "simulated_liver_data_FULL.rds")

# 2. Save the Stan Input. 
# This represents the "Real World" data where Sarcopenia (U) is missing/unmeasured.
# saveRDS(stan_data_list, "stan_liver_data_input.rds")

###--------------------------------------------------------------------------###
### 8. Narrative Sanity Checks (Print to Console)
###--------------------------------------------------------------------------###

# print("--- SIMULATION REPORT ---")

# 1. Naive Effect (The Confounded Estimate)
# Because Sarcopenic people (U) got treated more often, and Sarcopenia causes high MELD,
# the nutrition protocol might look like it DOESN'T work (or works poorly) in a naive view.

# print("Naive Average Treatment Effect (Observed Difference):")
# print(mean(d$Y.obs[d$A==1]) - mean(d$Y.obs[d$A==0]))

# 2. The Truth (The Causal Effect)
# We know nutrition helps (-3 base). 

# print("True Average Treatment Effect (Counterfactual):")
# print(mean(d$true.cate))

# 3. The Sex Disparity (The Motivation for Part 1)
# Females should show a stronger benefit (-5ish) compared to Males (-3ish).
# The Part 1 Models (Priors) will struggle to see this if the subgroup N is small.

# print("True Benefit by Sex (Females should benefit more):")
# print(aggregate(true.cate ~ sex, data=d, mean))

# hist(d$Y.obs, main="Simulated MELD 3.0 Scores", xlab="MELD Score", col="salmon", border="white")