library(evd)
library(truncnorm)
source(file.path(dir_path, "rjmcmc_common/utilities.r"))
source(file.path(dir_path, "rjmcmc_gumbel/single_acceptance.r"))



mcmc_single_gumbel <- function(iter, X,
                               xi, kappa, a, b,
                               mu, tau,
                               MHsd_mu = 0.5, MHsd_tau = 0.01,
                               verbose = TRUE) {
  parameter <- matrix(0, iter, 2)
  acceptance_count_mu <- 0
  acceptance_count_tau <- 0
  n <- length(X)

  progress_v <- floor(iter * c(1:10) * 0.1)
  for (i in 1:iter) {
    if (verbose) {
      if (i %in% progress_v) {
        tcat("----------------\n")
        cat("Iter:", i, "\n")
      }
    }
    # update mu
    new_mu <- rnorm(1, mean = mu, sd = MHsd_mu)
    if(is.na(new_mu)){
      cat("i: ", i, " mean:", mu, " sd:", MHsd_mu, "\n")
    } 
    acceptance_rate <- mu_acceptance(mu, new_mu, X, n, tau, xi, kappa)
    if (!is.na(acceptance_rate) && runif(1) < acceptance_rate) {
      mu <- new_mu
      acceptance_count_mu <- acceptance_count_mu + 1
    }

    # update tau
    new_tau <- rtruncnorm(
      n = 1, a = left_boundary_TN, b = Inf,
      mean = tau, sd = MHsd_tau
    )
    acceptance_rate <- tau_acceptance(tau, new_tau, X, n, mu, a, b, MHsd_tau)
    if (!is.na(acceptance_rate) && runif(1) < acceptance_rate) {
      tau <- new_tau
      acceptance_count_tau <- acceptance_count_tau + 1
    }

    parameter[i, ] <- c(mu, tau)
  }

  return(list(
    parameter = parameter,
    mu_acc_rate = acceptance_count_mu / iter,
    tau_acc_rate = acceptance_count_tau / iter
  ))
}
