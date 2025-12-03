rm(list=ls())
# setwd("~/Documents/RStudio/research/bayes-priors-sens/test/") 

source("R/02_generate.r")
# generate fake data
N <- 100 
Y <- rnorm(N, 1.6, 0.2)
hist(Y)

# compile model
library(rstan)

model <- stan_model('model.stan')

# pass data to stan and run 
# data list must be named list
options(mc.cores=4)
fit <- sampling(model, list(N=N, Y=Y), iter=200, chains = 4)

print(fit)
# Inference for Stan model: anon_model.
# 4 chains, each with iter=200; warmup=100; thin=1; 
# post-warmup draws per chain=100, total post-warmup draws=400.
# 
# mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
# mu      1.60    0.00 0.02   1.56   1.59   1.60   1.62   1.64   347 1.00
# sigma   0.21    0.00 0.01   0.18   0.20   0.21   0.21   0.24   144 1.02
# lp__  106.56    0.06 0.81 104.28 106.23 106.80 107.12 107.36   177 1.00
# 
# Samples were drawn using NUTS(diag_e) at Fri Nov 21 12:46:02 2025.
# For each parameter, n_eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor on split chains (at convergence, Rhat=1).

# extract and graph parameters 
params <- extract(fit)
# [1] "mu"    "sigma" "lp__" 

hist(params$mu)
hist(params$sigma)

library(shinystan)
launch_shinystan(fit)
