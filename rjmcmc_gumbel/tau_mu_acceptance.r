source(file.path(dir_path, "rjmcmc_common/utilities.r"))




mu_acceptance <- function(mu, new_mu, X, n, tau, xi, kappa) {
  # Likelihood ratio
  log_likelihood_ratio <- -tau * (sum(X) - n * new_mu) + tau * (sum(X) - n * mu)
  for (xt in X) {
    log_likelihood_ratio <- log_likelihood_ratio -
      min(max(exp(-tau * (xt - new_mu)), threshold), 1 / threshold) +
      min(max(exp(-tau * (xt - mu)), threshold), 1 / threshold)
  }
  # Prior ratio
  log_prior_ratio <- dnorm(new_mu, xi, 1 / sqrt(kappa), log = TRUE) -
                     dnorm(mu, xi, 1 / sqrt(kappa), log = TRUE)
  # Accept ratio
  log_r <- log_likelihood_ratio + log_prior_ratio
  r <- min(1, exp(log_r))
  return(r)
}

tau_acceptance <- function(tau, new_tau, X, n, mu, a, b, MHsd_tau) {
  # Proposal ratio (Jacobian)
  log_J_ratio <- log(1 - pnorm(-tau / MHsd_tau)) - log(1 - pnorm(-new_tau / MHsd_tau))
  # Likelihood ratio
  log_likelihood_ratio <- n * log(new_tau) - n * log(tau) -
    min(max(new_tau * (sum(X) - n * mu), threshold), 1 / threshold) +
    min(max(tau * (sum(X) - n * mu), threshold), 1 / threshold)
  for (xt in X) {
    log_likelihood_ratio <- log_likelihood_ratio -
      exp(-new_tau * (xt - mu)) +
      exp(-tau * (xt - mu))
  }
  # Prior ratio
  log_prior_ratio <- dgamma(new_tau, a, b, log = TRUE) -
                     dgamma(tau, a, b, log = TRUE)
  # Accept ratio
  log_r <- log_likelihood_ratio + log_prior_ratio - log_J_ratio
  r <- min(1, exp(log_r))
  return(r)
}
