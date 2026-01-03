source(file.path(dir_path, "rjmcmc_common/utilities.r"))

split_combine_acceptance <- function(is_split, X, k,
                                     Z, mu, tau, Xv, nv, pi,
                                     new_Z, new_mu, new_tau, new_Xv, new_nv, new_pi,
                                     xi, kappa, a, b, lambda, delta,
                                     u1, u2, u3,
                                     j1, j2, j.star,
                                     birth_pr, death_pr) {
  # Compute log acceptance ratio for split/combine move
  log_a <- 0
  n <- length(Z)

  # Likelihood ratio for affected components
  if (is_split) {
    for (t in c(1:n)) {
      if (new_Z[t] == j1 || new_Z[t] == j2) {
        log_a <- log_a + dnorm(X[t],
          mean = new_mu[new_Z[t]],
          sd = 1 / sqrt(new_tau[new_Z[t]]), log = TRUE
        )
      }
      if (Z[t] == j.star) {
        log_a <- log_a - dnorm(X[t],
          mean = mu[Z[t]],
          sd = 1 / sqrt(tau[Z[t]]), log = TRUE
        )
      }
    }
  } else {
    for (t in c(1:n)) {
      if (Z[t] == j1 || Z[t] == j2) {
        log_a <- log_a + dnorm(X[t],
          mean = mu[Z[t]],
          sd = 1 / sqrt(tau[Z[t]]), log = TRUE
        )
      }
      if (new_Z[t] == j.star) {
        log_a <- log_a - dnorm(X[t],
          mean = new_mu[new_Z[t]],
          sd = 1 / sqrt(new_tau[new_Z[t]]), log = TRUE
        )
      }
    }
  }


  # Extract relevant parameters for the move
  if (is_split) {
    pi_j1 <- new_pi[j1]
    pi_j2 <- new_pi[j2]
    pi_jstar <- pi[j.star]
    n_j1 <- new_nv[j1]
    n_j2 <- new_nv[j2]
    n_jstar <- nv[j.star]
    mu_j1 <- new_mu[j1]
    mu_j2 <- new_mu[j2]
    mu_jstar <- mu[j.star]
    tau_j1 <- new_tau[j1]
    tau_j2 <- new_tau[j2]
    tau_jstar <- tau[j.star]
  } else {
    pi_j1 <- pi[j1]
    pi_j2 <- pi[j2]
    pi_jstar <- new_pi[j.star]
    n_j1 <- nv[j1]
    n_j2 <- nv[j2]
    n_jstar <- new_nv[j.star]
    mu_j1 <- mu[j1]
    mu_j2 <- mu[j2]
    mu_jstar <- new_mu[j.star]
    tau_j1 <- tau[j1]
    tau_j2 <- tau[j2]
    tau_jstar <- new_tau[j.star]
  }

  # Prior ratio: mixture weights
  line1 <- log(lambda) - log(k + 1) - lbeta(delta, k * delta)
  line1 <- line1 + (delta + n_j1 - 1) * log(pi_j1)
  line1 <- line1 + (delta + n_j2 - 1) * log(pi_j2)
  line1 <- line1 - (delta + n_jstar - 1) * log(pi_jstar)
  log_a <- log_a + line1

  # Prior ratio: means
  line2 <- log(k + 1) + (1 / 2) * log(kappa / (2 * 3.1415))
  line2 <- line2 - (1 / 2) * kappa * ((mu_j1 - xi)^2 +
    (mu_j2 - xi)^2 - (mu_jstar - xi)^2)
  log_a <- log_a + line2

  # Prior ratio: precisions
  line3 <- a * log(b) - lgamma(a)
  line3 <- line3 + (a - 1) * log(tau_j1 * tau_j2 / tau_jstar)
  line3 <- line3 - b * (tau_j1 + tau_j2 - tau_jstar)
  log_a <- log_a + line3

  # Proposal ratio and Jacobian term
  line4 <- log(death_pr[k + 1] / birth_pr[k])
  line4 <- line4 - dbeta(u1, 2, 2, log = TRUE)
  line4 <- line4 - dbeta(u2, 2, 2, log = TRUE)
  line4 <- line4 - dbeta(u3, 1, 1, log = TRUE)
  log_a <- log_a + line4

  this_pi <- if (is_split) new_pi else pi
  this_mu <- if (is_split) new_mu else mu
  this_tau <- if (is_split) new_tau else tau
  this_Z <- if (is_split) new_Z else Z
  lineP <- 0
  for (t in seq_len(n)) {
    if (this_Z[t] == j1) {
      part1 <- this_pi[j1] * max(dnorm(X[t],
        mean = this_mu[j1],
        sd = 1 / sqrt(this_tau[j1])
      ), threshold)
      part2 <- this_pi[j2] * max(dnorm(X[t],
        mean = this_mu[j2],
        sd = 1 / sqrt(this_tau[j2])
      ), threshold)
      lineP <- lineP + log(part1 + part2) - log(part1)
    }
    if (this_Z[t] == j2) {
      part1 <- this_pi[j1] * max(dnorm(X[t],
        mean = this_mu[j1],
        sd = 1 / sqrt(this_tau[j1])
      ), threshold)
      part2 <- this_pi[j2] * max(dnorm(X[t],
        mean = this_mu[j2],
        sd = 1 / sqrt(this_tau[j2])
      ), threshold)
      lineP <- lineP + log(part1 + part2) - log(part2)
    }
  }
  log_a <- log_a + lineP



  # Jacobian, reparametrization and symmetry corrections
  line5 <- log(pi_jstar)
  line5 <- line5 - log(u2 * u3 * (1 - u2^2) * (1 - u3))
  line5 <- line5 + log(abs(mu_j2 - mu_j1)) + log(tau_j1 * tau_j2 / tau_jstar)
  log_a <- log_a + line5

  aaa <- exp(log_a)
  final_a <- ifelse(is_split, aaa, 1 / aaa)
  acceptance_r <- min(1, final_a)
  return(acceptance_r)
}

birth_death_acceptance <- function(
    is_birth, k, k0,
    delta, lambda,
    pi_jstar, n,
    birth_pr, death_pr) {
  # Compute log acceptance ratio for birth/death move
  log_a <- 0

  # Prior ratio for mixture weights
  part1 <- log(lambda) - lbeta(delta, k * delta)
  log_a <- log_a + part1

  # Prior and likelihood terms for new/deleted component
  part2 <- (delta - 1) * log(pi_jstar)
  part2 <- part2 + (n + k * delta) * log(1 - pi_jstar)
  log_a <- log_a + part2

  # Proposal ratio
  part3 <- log(death_pr[k + 1] / birth_pr[k])
  part3 <- part3 - log(k0)
  part3 <- part3 - dbeta(pi_jstar, 1, k, log = TRUE)
  log_a <- log_a + part3

  aaa <- exp(log_a)
  final_a <- ifelse(is_birth, aaa, 1 / aaa)
  acceptance_r <- min(1, final_a)
  return(acceptance_r)
}
