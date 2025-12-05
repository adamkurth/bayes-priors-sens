###--------------------------------------------------------------------------###
### Project: Quantifying Uncertainty in Causal Effects
### Utility Functions
###--------------------------------------------------------------------------###

#' Extract comprehensive Stan diagnostics
#' 
#' @param fit A stanfit object
#' @return A list with diagnostic information
extract_stan_diagnostics <- function(fit) {
  
  # Get sampler parameters
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  
  # Divergences
  n_divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  
  # Max treedepth
  n_max_treedepth <- sum(sapply(sampler_params, function(x) 
    sum(x[, "treedepth__"] >= 10)))
  
  # Energy Bayesian Fraction of Missing Information (E-BFMI)
  ebfmi <- sapply(sampler_params, function(x) {
    energy <- x[, "energy__"]
    if (length(energy) > 1) {
      numer <- sum(diff(energy)^2) / length(energy)
      denom <- var(energy)
      return(numer / denom)
    }
    return(NA)
  })
  
  # Summary statistics
  fit_summary <- summary(fit)$summary
  
  # R-hat
  rhat_vals <- fit_summary[, "Rhat"]
  rhat_max <- max(rhat_vals, na.rm = TRUE)
  rhat_warn <- sum(rhat_vals > 1.01, na.rm = TRUE)
  
  # Effective sample size
  neff_vals <- fit_summary[, "n_eff"]
  neff_min <- min(neff_vals, na.rm = TRUE)
  neff_ratio <- neff_vals / (nrow(as.matrix(fit)))
  neff_warn <- sum(neff_ratio < 0.1, na.rm = TRUE)
  
  return(list(
    n_divergent = n_divergent,
    n_max_treedepth = n_max_treedepth,
    ebfmi = ebfmi,
    ebfmi_low = sum(ebfmi < 0.2, na.rm = TRUE),
    rhat_max = rhat_max,
    rhat_warn = rhat_warn,
    neff_min = neff_min,
    neff_warn = neff_warn,
    all_ok = (n_divergent == 0) & (n_max_treedepth == 0) & 
             (rhat_max < 1.01) & (sum(ebfmi < 0.2, na.rm = TRUE) == 0)
  ))
}

#' Create a summary table for LaTeX
#' 
#' @param df Data frame to format
#' @param digits Number of decimal places
#' @return Formatted data frame
format_for_latex <- function(df, digits = 2) {
  df %>%
    mutate(across(where(is.numeric), ~round(., digits)))
}

#' Calculate posterior probability of direction
#' 
#' @param samples Vector of posterior samples
#' @param direction Either "positive" or "negative"
#' @return Probability that samples are in the specified direction
prob_direction <- function(samples, direction = "negative") {
  if (direction == "negative") {
    mean(samples < 0)
  } else {
    mean(samples > 0)
  }
}

#' Calculate the Highest Density Interval (HDI)
#' 
#' @param samples Vector of posterior samples
#' @param prob Probability mass (default 0.95)
#' @return Named vector with lower and upper bounds
hdi <- function(samples, prob = 0.95) {
  sorted_samples <- sort(samples)
  n <- length(sorted_samples)
  n_included <- ceiling(prob * n)
  n_ci <- n - n_included + 1
  
  ci_widths <- sorted_samples[(n_included):n] - sorted_samples[1:n_ci]
  min_idx <- which.min(ci_widths)
  
  c(lower = sorted_samples[min_idx], 
    upper = sorted_samples[min_idx + n_included - 1])
}

#' Create forest plot data from stanfit
#' 
#' @param fit A stanfit object
#' @param truth_data Data frame with true CATEs
#' @param par_name Name of the parameter (default "CATE")
#' @return Data frame ready for ggplot
prepare_forest_plot_data <- function(fit, truth_data, par_name = "CATE") {
  
  # Extract summary
  cate_summary <- summary(fit, pars = par_name)$summary %>%
    as.data.frame() %>%
    tibble::rownames_to_column("param")
  
  # Get subgroup mapping
  subgroup_map <- truth_data %>%
    dplyr::select(subgroup_id, subgroup_label, true.cate, sex, ascites, age.high) %>%
    dplyr::distinct() %>%
    dplyr::arrange(subgroup_id) %>%
    dplyr::mutate(
      Sex = ifelse(sex == 1, "Female", "Male"),
      Ascites = ifelse(ascites == 1, "Present", "Absent"),
      Age = ifelse(age.high == 1, "≥60", "<60")
    )
  
  # Combine
  result <- cbind(subgroup_map, cate_summary)
  
  return(result)
}

#' Compute sensitivity bounds for unmeasured confounding
#' 
#' @param naive_estimate The naive (potentially biased) estimate
#' @param gamma_range Vector of selection bias values
#' @param beta_range Vector of outcome bias values
#' @return Data frame with adjusted estimates
compute_sensitivity_bounds <- function(naive_estimate, gamma_range, beta_range) {
  
  expand.grid(gamma = gamma_range, beta = beta_range) %>%
    mutate(
      # Approximate bias adjustment (simplified Rosenbaum bounds)
      bias_adjustment = gamma * beta * 0.1,  # Scaling factor
      adjusted_lower = naive_estimate - bias_adjustment,
      adjusted_upper = naive_estimate + bias_adjustment,
      robustness = abs(naive_estimate) > bias_adjustment
    )
}

#' Generate a diagnostic report string
#' 
#' @param diagnostics List from extract_stan_diagnostics
#' @return Character string with formatted report
format_diagnostics_report <- function(diagnostics) {
  paste0(
    "Stan Diagnostics Summary:\n",
    "-------------------------\n",
    sprintf("Divergent transitions: %d\n", diagnostics$n_divergent),
    sprintf("Max treedepth exceeded: %d\n", diagnostics$n_max_treedepth),
    sprintf("Low E-BFMI chains: %d\n", diagnostics$ebfmi_low),
    sprintf("Max R-hat: %.4f (warning if > 1.01)\n", diagnostics$rhat_max),
    sprintf("Parameters with R-hat > 1.01: %d\n", diagnostics$rhat_warn),
    sprintf("Min n_eff: %.0f\n", diagnostics$neff_min),
    sprintf("Parameters with low n_eff ratio: %d\n", diagnostics$neff_warn),
    sprintf("\nOverall status: %s\n", 
            ifelse(diagnostics$all_ok, "✓ ALL OK", "⚠ ISSUES DETECTED"))
  )
}