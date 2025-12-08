###--------------------------------------------------------------------------###
### Project: Quantifying Uncertainty in Causal Effects
### Utility Functions
###--------------------------------------------------------------------------###

#' Extract comprehensive Stan diagnostics
#' 
#' @param fit A stanfit object
#' @return A list with diagnostic information
extract.stan.diagnostics <- function(fit) {
  
  # Get sampler parameters
  sampler.params <- get_sampler_params(fit, inc_warmup = FALSE)
  
  # Divergences
  n.divergent <- sum(sapply(sampler.params, function(x) sum(x[, "divergent__"])))
  
  # Max treedepth
  n.max.treedepth <- sum(sapply(sampler.params, function(x) 
    sum(x[, "treedepth__"] >= 10)))
  
  # Energy Bayesian Fraction of Missing Information (E-BFMI)
  ebfmi <- sapply(sampler.params, function(x) {
    energy <- x[, "energy__"]
    if (length(energy) > 1) {
      numer <- sum(diff(energy)^2) / length(energy)
      denom <- var(energy)
      return(numer / denom)
    }
    return(NA)
  })
  
  # Summary statistics
  fit.summary <- summary(fit)$summary
  
  # R-hat
  rhat.vals <- fit.summary[, "Rhat"]
  rhat.max <- max(rhat.vals, na.rm = TRUE)
  rhat.warn <- sum(rhat.vals > 1.01, na.rm = TRUE)
  
  # Effective sample size
  neff.vals <- fit.summary[, "n_eff"]
  neff.min <- min(neff.vals, na.rm = TRUE)
  neff.ratio <- neff.vals / (nrow(as.matrix(fit)))
  neff.warn <- sum(neff.ratio < 0.1, na.rm = TRUE)
  
  return(list(
    n.divergent = n.divergent,
    n.max.treedepth = n.max.treedepth,
    ebfmi = ebfmi,
    ebfmi.low = sum(ebfmi < 0.2, na.rm = TRUE),
    rhat.max = rhat.max,
    rhat.warn = rhat.warn,
    neff.min = neff.min,
    neff.warn = neff.warn,
    all.ok = (n.divergent == 0) & (n.max.treedepth == 0) & 
             (rhat.max < 1.01) & (sum(ebfmi < 0.2, na.rm = TRUE) == 0)
  ))
}

#' Create a summary table for LaTeX
#' 
#' @param df Data frame to format
#' @param digits Number of decimal places
#' @return Formatted data frame
format.for.latex <- function(df, digits = 2) {
  df %>%
    mutate(across(where(is.numeric), ~round(., digits)))
}

#' Calculate posterior probability of direction
#' 
#' @param samples Vector of posterior samples
#' @param direction Either "positive" or "negative"
#' @return Probability that samples are in the specified direction
prob.direction <- function(samples, direction = "negative") {
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
  sorted.samples <- sort(samples)
  n <- length(sorted.samples)
  n.included <- ceiling(prob * n)
  n.ci <- n - n.included + 1
  
  ci.widths <- sorted.samples[(n.included):n] - sorted.samples[1:n.ci]
  min.idx <- which.min(ci.widths)
  
  c(lower = sorted.samples[min.idx], 
    upper = sorted.samples[min.idx + n.included - 1])
}

#' Create forest plot data from stanfit
#' 
#' @param fit A stanfit object
#' @param truth.data Data frame with true CATEs
#' @param par.name Name of the parameter (default "CATE")
#' @return Data frame ready for ggplot
prepare.forest.plot.data <- function(fit, truth.data, par.name = "CATE") {
  
  # Extract summary
  cate.summary <- summary(fit, pars = par.name)$summary %>%
    as.data.frame() %>%
    tibble::rownames_to_column("param")
  
  # Get subgroup mapping
  subgroup.map <- truth.data %>%
    dplyr::select(subgroup_id, subgroup_label, true.cate, sex, ascites, age.high) %>%
    dplyr::distinct() %>%
    dplyr::arrange(subgroup_id) %>%
    dplyr::mutate(
      Sex = ifelse(sex == 1, "Female", "Male"),
      Ascites = ifelse(ascites == 1, "Present", "Absent"),
      Age = ifelse(age.high == 1, "≥60", "<60")
    )
  
  # Combine
  result <- cbind(subgroup.map, cate.summary)
  
  return(result)
}

#' Compute sensitivity bounds for unmeasured confounding
#' 
#' @param naive.estimate The naive (potentially biased) estimate
#' @param gamma.range Vector of selection bias values
#' @param beta.range Vector of outcome bias values
#' @return Data frame with adjusted estimates
compute.sensitivity.bounds <- function(naive.estimate, gamma.range, beta.range) {
  
  expand.grid(gamma = gamma.range, beta = beta.range) %>%
    mutate(
      # Approximate bias adjustment (simplified Rosenbaum bounds)
      bias.adjustment = gamma * beta * 0.1,  # Scaling factor
      adjusted.lower = naive.estimate - bias.adjustment,
      adjusted.upper = naive.estimate + bias.adjustment,
      robustness = abs(naive.estimate) > bias.adjustment
    )
}

#' Generate a diagnostic report string
#' 
#' @param diagnostics List from extract.stan.diagnostics
#' @return Character string with formatted report
format.diagnostics.report <- function(diagnostics) {
  paste0(
    "Stan Diagnostics Summary:\n",
    "-------------------------\n",
    sprintf("Divergent transitions: %d\n", diagnostics$n.divergent),
    sprintf("Max treedepth exceeded: %d\n", diagnostics$n.max.treedepth),
    sprintf("Low E-BFMI chains: %d\n", diagnostics$ebfmi.low),
    sprintf("Max R-hat: %.4f (warning if > 1.01)\n", diagnostics$rhat.max),
    sprintf("Parameters with R-hat > 1.01: %d\n", diagnostics$rhat.warn),
    sprintf("Min n_eff: %.0f\n", diagnostics$neff.min),
    sprintf("Parameters with low n_eff ratio: %d\n", diagnostics$neff.warn),
    sprintf("\nOverall status: %s\n", 
            ifelse(diagnostics$all.ok, "✓ ALL OK", "⚠ ISSUES DETECTED"))
  )
}