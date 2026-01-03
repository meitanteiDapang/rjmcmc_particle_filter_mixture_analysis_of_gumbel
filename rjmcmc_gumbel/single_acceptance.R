source(file.path(dir_path, "rjmcmc_common/utilities.r"))
library(truncnorm)
library(evd)


mu_acceptance <- function(mu, new_mu, X, n, tau, xi, kappa) {
  # Log-likelihood ratio for Gumbel
  log_likelihood_ratio <-
    -tau * (sum(X) - n * new_mu) - sum(exp(-tau * (X - new_mu))) +
    tau * (sum(X) - n * mu) + sum(exp(-tau * (X - mu)))
  # Prior ratio (normal)
  log_prior_ratio <-
    dnorm(new_mu, xi, 1 / sqrt(kappa), log = TRUE) -
    dnorm(mu, xi, 1 / sqrt(kappa), log = TRUE)
  log_r <- log_likelihood_ratio + log_prior_ratio
  r <- min(1, exp(log_r))
  return(r)
}

tau_acceptance <- function(tau, new_tau, X, n, mu, a, b, MHsd_tau) {
  # Proposal density correction (truncated normal)
  log_J_ratio <- log(1 - pnorm(-tau / MHsd_tau)) -
    log(1 - pnorm(-new_tau / MHsd_tau))
  # Log-likelihood ratio for Gumbel
  log_likelihood_ratio <-
    n * log(new_tau) - new_tau * (sum(X) - n * mu) - sum(exp(-new_tau * (X - mu))) -
    n * log(tau) + tau * (sum(X) - n * mu) + sum(exp(-tau * (X - mu)))
  # Prior ratio (gamma)
  log_prior_ratio <-
    dgamma(new_tau, a, b, log = TRUE) -
    dgamma(tau, a, b, log = TRUE)
  log_r <- log_likelihood_ratio + log_prior_ratio - log_J_ratio
  r <- min(1, exp(log_r))
  return(r)
}
