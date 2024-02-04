############################################################################################
# @roman_avj                                                                    3 feb 2024
# hw2
############################################################################################
# libraries
library(tidyverse)

#### q8: estimate poisson/gamma parameters for Gilchrist data
# load data #
gilchrist <- "250231343032502313430325023134303"
gilchrist <- as.numeric(strsplit(gilchrist, "")[[1]])

# stan #
# compile
# stan_model <- file.path("hw2/stan", "gammpoisson_model.stan") |> cmdstanr::cmdstan_model()
stan_model <- file.path("hw2/stan", "gammpoisson_model_v2.stan") |> cmdstanr::cmdstan_model()
# data
stan_data <- list(N = length(gilchrist), x = gilchrist)
# fit
fit_stan <- stan_model$sample(
    data = stan_data, 
    chains = 4, 
    iter_warmup = 1000, 
    iter_sampling = 1000,
    seed = 42,
    parallel_chains = 4
    )

# diagnostics
fit_stan$cmdstan_diagnose()
fit_stan$summary()
fit_stan$draws(c("a", "b")) |> 
  posterior::as_draws_df() |> 
  ggplot(aes(.iteration, a)) +
  geom_line() +
  facet_wrap(~.chain, ncol = 1)

# draws for a and b
fit_stan$draws(c("a", "b")) |> 
    posterior::as_draws_df() |> 
    ggplot(aes(a, b)) +
    geom_point()