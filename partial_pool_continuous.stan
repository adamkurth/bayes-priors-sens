data {
  int<lower=0> N;          // Number of observations
  int<lower=0> Pv;         // Number of subgroups
  int<lower=0> Pw;         // Number of fixed covariates
  
  matrix[N, Pv] V;         // Matrix of subgroup indicators
  matrix[N, Pw] W;         // Matrix of fixed covariates (intercept, age, etc)
  vector[N] A;             // Treatment assignment (0 or 1)
  vector[N] Y;             // Outcome (Continuous MELD Score) allows to use normal() likelihood
  
  // -- HYPERPARAMETERS FOR PRIORS (Passed from R) --
  // This allows switching between Strong/Weak/Informative priors in R
  real prior_tau_sigma;    // Width of the prior on heterogeneity (tau)
  real prior_mean_mu;      // Prior mean for the Average Treatment Effect
}

parameters {
  vector[Pw] beta_w;       // Fixed covariate effects
  vector[Pv] beta_v;       // Subgroup main effects (intercept shifts)
  
  // CATE PARAMETERS
  real mu;                 // Average Treatment Effect (Global)
  real<lower=0> tau;       // Standard Deviation of CATEs (Heterogeneity)
  vector[Pv] theta_raw;    // Non-centered parameterization for better sampling
  real<lower=0> sigma;     // Residual Error (Sigma_Y) - NEEDED FOR CONTINUOUS MODEL
}

transformed parameters {
  // Construct the actual CATEs (theta) from the global mean + deviation
  // theta[v] = mu + tau * raw_deviation
  vector[Pv] theta;
  theta = mu + tau * theta_raw;
}

model {
  // --- 1. PRIORS ---
  
  // Nuisance parameters (Weakly informative)
  beta_w ~ normal(0, 10);
  beta_v ~ normal(0, 10);
  sigma ~ cauchy(0, 5);
  
  // PARTIAL POOLING PRIORS (The core of your research)
  
  // The Average Effect
  mu ~ normal(prior_mean_mu, 10); 
  
  // The Heterogeneity (Tau) - This controls the pooling strength
  // Small tau -> Strong Pooling (Everything shrinks to mean)
  // Large tau -> Weak Pooling (Subgroups can vary wildly)
  tau ~ normal(0, prior_tau_sigma); 
  
  // The deviations for each subgroup
  theta_raw ~ std_normal(); 
  
  // --- 2. LIKELIHOOD ---
  // Linear Regression: Y = W*beta_w + V*beta_v + (V*theta)*A
  // Note: (V * theta) picks the correct CATE for the subject's subgroup
  
  vector[N] treatment_effect_per_person;
  treatment_effect_per_person = (V * theta) .* A;     // .* is element-wise multiplication
  
  Y ~ normal( W*beta_w + V*beta_v + treatment_effect_per_person, sigma );
}

generated quantities {
  // recover the CATEs for easy extraction in R
  vector[Pv] CATE;
  CATE = theta;
  
  // Posterior predictive checks  
  vector[N] y_rep;
  for (n in 1:N) {
     real y_hat = (W[n]*beta_w + V[n]*beta_v + (V[n]*theta)*A[n])[1];
     y_rep[n] = normal_rng(y_hat, sigma);
  }

}




