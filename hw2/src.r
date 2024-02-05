############################################################################################
# @roman_avj                                                                    3 feb 2024
# hw2
############################################################################################
# libraries
library(tidyverse)

#### functions ####
log_sum_exp <- function(x){
    max_x  <- max(x)
    max_x + log(sum(exp(x - max_x)))
}

#### q1: inference over poisson parameter with discrete prior ####
# Y  <- number of failures on t time, Y ~ Poisson(lambda*t)
# lambda ~ Discrete(1/2, 2/2, 3/2, 4/2, 5/2, 6/2)
# p(lambda) \in {0.1, 0.2, 0.3, 0.2, 0.15, 0.05}
table_prob_lambda  <- tibble(
    lambda = c(0.5, 1, 1.5, 2, 2.5, 3),
    prob = c(0.1, 0.2, 0.3, 0.2, 0.15, 0.05)
)

# posterior after observing y = 12 and t = 6
TIME <- 6
Y_OBS  <- 12
table_prob_lambda  <- table_prob_lambda |> 
    mutate(
        # get nuclues of posterior
        log_posterior = log(prob) + dpois(Y_OBS, lambda*TIME, log = TRUE),
        # normalize
        posterior = exp(log_posterior - log_sum_exp(log_posterior))
    )

# predictive posterior for y = 0 and t = 7
TIME_PRED <- 7
Y_PRED  <- 0
table_prob_lambda  <- table_prob_lambda |> 
    mutate(
        log_predictive = log(posterior) + dpois(Y_PRED, lambda*TIME_PRED, log = TRUE),
    )
y_predictive  <- exp(log_sum_exp(table_prob_lambda$log_predictive))
cat("Predictive posterior for Y = 0 and t = 7: ", y_predictive, "\n")
table_prob_lambda

#### q2: creditibility interval ####
# p ~ beta(23, 8), find 90% creditibility interval
# dplot of beta(23, 8)
x  <- seq(0, 1, 0.01)
y  <- dbeta(x, 23, 8)
# get mode 
mode  <- (23 - 1) / (23 + 8 - 2)

plot(x, y, type = "l", xlab = "p", ylab = "f(p)", main = "Beta(23, 8)")
# add mode
abline(v = mode, col = "red")

# trial 1: using quantiles
p_005 <- qbeta(0.05, 23, 8)
p_095 <- qbeta(0.95, 23, 8)
cat("Creditibility interval using quantiles: ", p_005, p_095, "\n")
cat("Length: ", p_095 - p_005, "\n")

# plot creditibility interval
plot(x, y, type = "l", xlab = "p", ylab = "f(p)", main = "Beta(23, 8)")
abline(v = p_005, col = "blue")
abline(v = p_095, col = "blue")
abline(v = mode, col = "red")

# trial 2: using binary search
# binary search #
find_credibility_interval <- function(alpha, beta, credibility = 0.9, tol = 1e-4, maxiter = 20) {
    # Binary search
    mode <- (alpha - 1) / (alpha + beta - 2)
    # Initial state
    left <- 0
    right <- dbeta(mode, alpha, beta)
    m <- (right - left) / 2

    # get p_005 and p_095 using binary search
    f <- function(x) dbeta(x, alpha, beta) - m
    p_005 <- uniroot(f, c(0, mode))$root
    p_095 <- uniroot(f, c(mode, 1))$root

    # Find p_005 and p_095 using binary search such that pbeta(p_095, alpha, beta) - pbeta(p_005, alpha, beta) = credibility
    i  <- 0
    # plot initial state
    plot(x, y, type = "l", xlab = "p", ylab = "f(p)", main = "Beta(23, 8)")
    abline(v = p_005, col = "blue")
    abline(v = p_095, col = "blue")
    abline(h = m, col = "blue")
    text(0.1, m, i)
    abline(v = mode, col = "red")

    # flags
    ci  <- pbeta(p_095, alpha, beta) - pbeta(p_005, alpha, beta)
    while ((abs(ci - credibility) > tol) & (maxiter > i)) {
        # Update line
        cat("----------------\n")
        cat("Iteration: ", i + 1, "\n")
        cat("Previous credibility interval: " , ci, "\n")
        if (ci > credibility) {
            cat("above \n")
            left <- m
            m <- m + (right - left) / 2
        } else {
            cat("below \n")
            right <- m
            m <- m - (right - left) / 2
        }

        # Update flags
        f <- function(x) dbeta(x, alpha, beta) - m
        p_005 <- uniroot(f, c(0, mode))$root
        p_095 <- uniroot(f, c(mode, 1))$root
        ci  <- pbeta(p_095, alpha, beta) - pbeta(p_005, alpha, beta)
        i  <- i + 1

        # plot 
        x <- seq(0, 1, length.out = 1000)
        y <- dbeta(x, alpha, beta)
        abline(v = p_005, col = "blue")
        abline(v = p_095, col = "blue")
        abline(h = m, col = "blue")
        text(0.1, m, i)
    }
    cat(" ================ \n")
    cat("Num iter: ", i, "\n")
    cat("New Probability interval: ", pbeta(p_095, alpha, beta) - pbeta(p_005, alpha, beta), "\n")
    cat("Credibility interval using quantiles: ", p_005, p_095, "\n")
    cat("Length: ", p_095 - p_005, "\n") 

    return(list(lower = p_005, upper = p_095))
}

# Usage
credibility_interval <- find_credibility_interval(23, 8, maxiter = 20)

## compare with quantiles
p_005 <- qbeta(0.05, 23, 8)
p_095 <- qbeta(0.95, 23, 8)
cat("Creditibility interval using quantiles: ", p_005, p_095, "\n")
cat("Length: ", p_095 - p_005, "\n")
cat("Creditibility interval using binary search: ", credibility_interval$lower, credibility_interval$upper, "\n")
cat("Length: ", credibility_interval$upper - credibility_interval$lower, "\n")

# plot creditibility interval
plot(x, y, type = "l", xlab = "p", ylab = "f(p)", main = "Beta(23, 8)") 
abline(v = p_005, col = "blue")
abline(v = p_095, col = "blue")
abline(v = credibility_interval$lower, col = "black")
abline(v = credibility_interval$upper, col = "black")
abline(v = mode, col = "red")

## proba of p exceeding 0.6
p_exceed_06  <- 1 - pbeta(0.6, 23, 8)
cat("Probability of p exceeding 0.6: ", p_exceed_06, "\n")

#### q3: normal distribution inference ####
# y|mu ~ N(mu, 100)
y_obs  <- c(38.6, 42.4, 57.5, 40.5, 51.7, 67.1, 33.4, 60.9, 64.1, 40.1, 40.7, 6.4)

table_prob_mu  <- tibble(
    mu = seq(20, 70, 10),
    prob = c(0.1, 0.15, 0.25, 0.25, 0.15, 0.1)
    ) |> group_by(
        # group by mu
        mu
    ) |> mutate(
        # posterior
        log_posterior = log(prob) + sum(dnorm(y_obs, mean=mu, sd=sqrt(100), log = TRUE)),
    ) |> ungroup(
        
    ) |> 
        mutate(
            # normalize
            posterior = exp(log_posterior - log_sum_exp(log_posterior))
        )

# q3.1: get mean
mean_mu  <- sum(exp(log(table_prob_mu$mu) + log(table_prob_mu$posterior)))
cat("Mean: ", mean_mu, "\n")

# q3.2: look probabilities
table_prob_mu

# q3.3 IC at 80%
# plot posterior as lolliplot
plot(table_prob_mu$mu, table_prob_mu$posterior, type = "h", xlab = "mu", ylab = "posterior", main = "Posterior distribution")
# add points
points(table_prob_mu$mu, table_prob_mu$posterior, pch = 19)

#### q4: cauchy distribution inference ####
# y|theta ~ Cauchy(theta, 1)
# theta ~ Unif(-2, 12)

# q4.1: grid for theta
u1  <- -2
u2  <- 12
theta  <- seq(u1, u2, 0.1)
proba_theta  <- 1 / (u2 - u1)
# plot
plot(theta, rep(proba_theta, length(theta)), type = "h", xlab = "theta", ylab = "p(theta)", main = "Uniform distribution")

# 4.2 posterior over grid
y_obs  <- c(0,10,9,8,11,3,3,8,8,11)
table_prob_theta_cauchy  <- tibble(
        theta = theta,
        prob = rep(proba_theta, length(theta))
    ) |> group_by(
        theta
    ) |> mutate(
        # posterior
        log_posterior = log(prob) + sum(dcauchy(y_obs, location = theta, scale = 1, log = TRUE))
    ) |> ungroup(

    ) |> mutate(
        posterior = exp(log_posterior - log_sum_exp(log_posterior))
    )
table_prob_theta_cauchy

# plot density for theta posterior
plot(table_prob_theta_cauchy$theta, table_prob_theta_cauchy$posterior, type = "h", xlab = "theta", ylab = "posterior", main = "Posterior distribution")
# add points
points(table_prob_theta_cauchy$theta, table_prob_theta_cauchy$posterior, pch = 19)

# 4.4: get posterior mean
mean_theta  <- sum(exp(log(table_prob_theta_cauchy$theta) + log(table_prob_theta_cauchy$posterior)), na.rm = TRUE)
cat("Mean: ", mean_theta, "\n")

# 4.5: get posterior std
std_theta  <- sqrt(sum(exp(2 * log(table_prob_theta_cauchy$theta) + log(table_prob_theta_cauchy$posterior)), na.rm = TRUE) - mean_theta^2)
cat("Std: ", std_theta, "\n")

#### q5: bayesian robustness ####
# prior1 ~ beta(100, 100)
# prior2 ~ 0.9 * beta(500, 500) + 0.1 * beta(1, 1)

# q5.1: simulate priors
set.seed(42)
N  <- 1000
p1  <- rbeta(N, 100, 100)

beta1  <- rbeta(N, 500, 500)
beta2  <- rbeta(N, 1, 1)
p2  <- 0.9 * beta1 + 0.1 * beta2
# plot both priors
tibble(p1, p2) |> 
    pivot_longer(cols = c(p1, p2)) |> 
    ggplot(aes(value, fill = name)) +
    geom_density(alpha = 0.5) +
    ggtitle("Priors") +
    theme_minimal()

# look mean for both priors
tibble(p1, p2) |> colMeans()

# prop of pi between 0.44 and 0.56
p1_44_56 <- mean((p1 > 0.44) & (p1 < 0.56))
p2_44_56 <- mean((p2 > 0.44) & (p2 < 0.56))
cat("Proportion of p1 between 0.44 and 0.56: ", p1_44_56, "\n")
cat("Proportion of p2 between 0.44 and 0.56: ", p2_44_56, "\n")

# q5.2: update priors with data
sum_yobs  <- 45
n_trials  <- 100

# update priors
p1_post  <- rbeta(N, 100 + sum_yobs, 100 + n_trials - sum_yobs)
p2_post  <- 0.9 * rbeta(N, 500 + sum_yobs, 500 + n_trials - sum_yobs) + 0.1 * rbeta(N, 1 + sum_yobs, 1 + n_trials - sum_yobs)
# plot both posteriors
tibble(p1_post, p2_post) |> 
    pivot_longer(cols = c(p1_post, p2_post)) |> 
    ggplot(aes(value, fill = name)) +
    geom_density(alpha = 0.5) +
    ggtitle("Posteriors") +
    theme_minimal()

# get posterior probability intervals at 90%
tibble(
    interval = c("lower", "upper"),
    p1_post = quantile(p1_post, c(0.05, 0.95)),
    p2_post = quantile(p2_post, c(0.05, 0.95))
)

# q5.3: update priors with data, sum_yobs = 30
sum_yobs  <- 30
n_trials  <- 100

# update priors
p1_post  <- rbeta(N, 100 + sum_yobs, 100 + n_trials - sum_yobs)
p2_post  <- 0.9 * rbeta(N, 500 + sum_yobs, 500 + n_trials - sum_yobs) + 0.1 * rbeta(N, 1 + sum_yobs, 1 + n_trials - sum_yobs)
# plot both posteriors
tibble(p1_post, p2_post) |> 
    pivot_longer(cols = c(p1_post, p2_post)) |> 
    ggplot(aes(value, fill = name)) +
    geom_density(alpha = 0.5) +
    ggtitle("Posteriors") +
    theme_minimal()

# get posterior probability intervals at 90%
tibble(
    interval = c("lower", "upper"),
    p1_post = quantile(p1_post, c(0.05, 0.95)),
    p2_post = quantile(p2_post, c(0.05, 0.95))
)

#### q6: grouped data ####
# X|mu ~ N(mu, 100)
# R|p ~ Binomial(n, p) with p = Pr(X > 70)

## trial 1: prior over p, not considering mu
# suppose p ~ beta(1, 1)
set.seed(42)
N  <- 10000
p1_sample <- rbeta(N, 1, 1)
hist(p1_sample, main = "Prior for p")

# get mu for p1
mu_sample  <- 70 - 10 * qnorm(p1_sample, mean = 0, sd = 1)
hist(mu_sample, main = "Prior for mu")

# observed data
num_rebased  <- 1
num_cars  <- 18
p1_posterior  <- rbeta(N, num_rebased + 1, num_cars - num_rebased + 1)
mu_posterior  <- 70 - 10 * qnorm(p1_posterior, mean = 0, sd = 1)
hist(p1_posterior, main = "Posterior for p")
hist(mu_posterior, main = "Posterior for mu")

# get mean of mu_posterior
mean_mu  <- mean(mu_posterior)
cat("Mean of mu_posterior: ", mean_mu, "\n")

# get predictive posterior for X > 80
p_x80  <- mean(mu_posterior > 80)
cat("Predictive posterior for X > 80: ", p_x80, "\n")

## trial 2: prioer over mu
# plain mu, mu ~ Unif(m), m > 0
# observed data
num_rebased  <- 1
num_cars  <- 18
mu_prior  <- 60:120
p1_given_mu <- pnorm(70, mean = mu_prior, sd = 10)
loglike_mu <- num_rebased * log(p1_given_mu) + (num_cars - num_rebased) * log(1 - p1_given_mu)
posterior_mu  <- exp(loglike_mu - log_sum_exp(loglike_mu))

# plot posterior
plot(mu_prior, posterior_mu, type = "h", xlab = "mu", ylab = "posterior", main = "Posterior distribution")
# add points
points(mu_prior, posterior_mu, pch = 19)
# plot max of posterior and add text
max_posterior  <- mu_prior[which.max(posterior_mu)]
abline(v = max_posterior, col = "red")
text(max_posterior + 5, max(posterior_mu), str_glue("Mode: {max_posterior}"))

# get mean of posterior
mean_posterior  <- sum(mu_prior * posterior_mu)
cat("Mean of posterior: ", mean_posterior, "\n")

# get predictive posterior for X > 80
p_x80  <- sum(posterior_mu[mu_prior > 80])
cat("Predictive posterior for X > 80: ", p_x80, "\n")

#### q7: behrens-fisher problem ####
# X ~ N(mu1, sigma1^2) and Y ~ N(mu2, sigma2^2) and X, Y independent
# prior of mu1, mu2, sigma1, sigma2 prop to 1/(sigma1^2 * sigma2^2)
x_obs  <- c(120, 107, 110, 116, 114, 111, 113, 117, 114, 112)
y_obs  <- c(110, 111, 107, 108, 110, 105, 107, 106, 111, 111)
mu1_prior  <- seq(100, 120, 1)
mu2_prior  <- seq(100, 120, 1)
sigma1_prior  <- seq(1, 6, 0.1)
sigma2_prior  <- seq(1, 6, 0.1)

# get posterior for mu1, mu2, sigma1, sigma2
n1 <- length(x_obs)
n2 <- length(y_obs)

# posterior for theta = (theta1, theta2) = ((mu1, sigma1), (mu2, sigma2))
loglike_theta_i  <- function(mu, sigma, z) {
    sum(dnorm(z, mean=mu, sd=sigma, log = TRUE)) - 2 * log(sigma)
}

# posterior for mu1, sigma1:
loglike_theta1 <- sapply(mu1_prior, function(mu1) {
    sapply(sigma1_prior, function(sigma1) {
        loglike_theta_i(mu1, sigma1, x_obs)
    })
}) |> apply(c(1, 2), sum)
posterior_theta1  <- exp(loglike_theta1 - log_sum_exp(loglike_theta1))

# posterior for mu2, sigma2:
loglike_theta2 <- sapply(mu2_prior, function(mu2) {
    sapply(sigma2_prior, function(sigma2) {
        loglike_theta_i(mu2, sigma2, y_obs)
    })
}) |> apply(c(1, 2), sum)
posterior_theta2  <- exp(loglike_theta2 - log_sum_exp(loglike_theta2))

# plot curves for mu1, sigma1 using loglike_theta1 and ggplot2
table_level_curves_theta1  <- tibble(
    mu1 = rep(mu1_prior, each = length(sigma1_prior)),
    sigma1 = rep(sigma1_prior, length(mu1_prior)),
    posterior = posterior_theta1 |> as.vector()
)

table_level_curves_theta1 |> ggplot(aes(mu1, sigma1)) +
    geom_raster(aes(fill = posterior)) +
    geom_contour(aes(z = posterior), color = "white") +
    theme_minimal() +
    ggtitle("Posterior for mu1, sigma1")

# plot curves for mu2, sigma2 using loglike_theta2 and ggplot2
table_level_curves_theta2  <- tibble(
    mu2 = rep(mu2_prior, each = length(sigma2_prior)),
    sigma2 = rep(sigma2_prior, length(mu2_prior)),
    posterior = posterior_theta2 |> as.vector()
)

table_level_curves_theta2 |> ggplot(aes(mu2, sigma2)) +
    geom_raster(aes(fill = posterior)) +
    geom_contour(aes(z = posterior), color = "white") +
    theme_minimal() +
    ggtitle("Posterior for mu2, sigma2")

# plot both curves in the same ggplot without raster
table_level_curves_theta1 |> 
    mutate(type = "male") |> 
    rename(mu = mu1, sigma = sigma1) |>
    bind_rows(
        table_level_curves_theta2 |> 
            mutate(type = "female") |> 
            rename(mu = mu2, sigma = sigma2)
    ) |> 
    ggplot(aes(mu, sigma)) +
    geom_contour(aes(z = posterior, color = type)) +
    theme_minimal() +
    ggtitle("Posterior of mu & sigma for male & female")

#### q8: estimate poisson/gamma parameters for Gilchrist data ####
## trial 1: using mcmc ##
# data
gilchrist <- "250231343032502313430325023134303"
gilchrist <- as.numeric(strsplit(gilchrist, "")[[1]])
# stan
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

## trial 2: using custom simulation (the goodone) ##
# poissongamma function
dpoisson_gamma <- function(x, a, b, log = FALSE) {
    if (log) {
        return((lgamma(x + a) - lgamma(a) - lgamma(x + 1)) + (a * log(b) - (x + a) * log(b + 1)))
    } else {
        return((gamma(x + a) / (gamma(a) * gamma(x + 1))) * (b^a / (b + 1)^(x + 1)))
    }
}

# data
gilchrist <- "250231343032502313430325023134303"
gilchrist <- as.numeric(strsplit(gilchrist, "")[[1]])

# prior of a, b  prop to 1/(a*b)^2
# transform theta1 = log(a), theta2 = log(b)
# posterior of theta1, theta2 prop to prod_i p(x_i|theta1, theta2) * 1/(a*b)^2, with a = exp(theta1), b = exp(theta2)
# prior grid
prior_theta1 <- seq(-2.5, 2.5, 0.01)
prior_theta2 <- seq(-2.5, 2.5, 0.01)

# get posterior for theta1, theta2
loglike_theta1_theta2 <- sapply(prior_theta1, function(theta1) {
    sapply(prior_theta2, function(theta2) {
        sum(dpoisson_gamma(gilchrist, exp(theta1), exp(theta2), log = TRUE)) - 2 * (theta1 + theta2)
    })
}) |> apply(c(1, 2), sum)

# get posterior
posterior_theta1_theta2  <- exp(loglike_theta1_theta2 - log_sum_exp(loglike_theta1_theta2))

# plot curves for theta1, theta2 using posterior
table_level_curves_theta1_theta2  <- tibble(
    theta1 = rep(prior_theta1, each = length(prior_theta2)),
    theta2 = rep(prior_theta2, length(prior_theta1)),
    posterior = posterior_theta1_theta2 |> as.vector()
)

table_level_curves_theta1_theta2 |> ggplot(aes(theta1, theta2)) +
    geom_raster(aes(fill = posterior)) +
    geom_contour(aes(z = posterior), color = "white") +
    theme_minimal() +
    ggtitle("Posterior for theta1, theta2")

# get credible interval at 90% using bonferroni correction
alpha  <- 0.1
n  <- length(posterior_theta1_theta2)
n_intervals  <- 2
bonferroni_correction  <- alpha / n_intervals
# get credible interval
posterior_theta1  <- posterior_theta1_theta2 |> apply(2, sum)
posterior_theta2  <- posterior_theta1_theta2 |> apply(1, sum)
# plot probability density for theta1, theta2
plot(prior_theta1, posterior_theta1, type = "h", xlab = "theta1", ylab = "posterior", main = "Posterior distribution")
plot(prior_theta2, posterior_theta2, type = "h", xlab = "theta2", ylab = "posterior", main = "Posterior distribution")

# get credible interval for theta1
theta1_005 <- prior_theta1[which.max(cumsum(posterior_theta1) > bonferroni_correction)]
theta1_095 <- prior_theta1[which.max(cumsum(posterior_theta1) > (1 - bonferroni_correction))]
cat("Credible interval for theta1: ", theta1_005, theta1_095, "\n")
cat("Creditibility interval for a: ", exp(theta1_005), exp(theta1_095), "\n")
# get credible interval for theta2
theta2_005 <- prior_theta2[which.max(cumsum(posterior_theta2) > bonferroni_correction)]
theta2_095 <- prior_theta2[which.max(cumsum(posterior_theta2) > (1 - bonferroni_correction))]
cat("Credible interval for theta2: ", theta2_005, theta2_095, "\n")
cat("Creditibility interval for b: ", exp(theta2_005), exp(theta2_095), "\n")

# add bonferroni intervals to plot
table_level_curves_theta1_theta2 |> ggplot(aes(theta1, theta2)) +
    geom_raster(aes(fill = posterior)) +
    geom_contour(aes(z = posterior), color = "white") +
    geom_vline(xintercept = c(theta1_005, theta1_095), color = "red") +
    geom_hline(yintercept = c(theta2_005, theta2_095), color = "red") +
    theme_minimal() +
    ggtitle("Posterior for theta1, theta2")





































