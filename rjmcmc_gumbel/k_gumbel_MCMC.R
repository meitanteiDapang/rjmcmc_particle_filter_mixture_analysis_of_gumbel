source(file.path(dir_path, "rjmcmc_gumbel/tau_mu_acceptance.r"))

mcmc_mixture_gumbel_k <- function(iter, X, input_k,
                                  xi, kappa, a, b, delta,
                                  mu, tau, pi, Z,
                                  MHsd_mu, MHsd_tau,
                                  verbose = TRUE) {
  n <- length(X)
  k <- input_k
  parameter <- matrix(0, iter, 3 * k)
  acceptance_count <- rep(0, 2 * k)
  Xv <- vector("list", k)
  nv <- rep(0, k)
  progress_v <- floor(iter * seq_len(10) / 10)

  for (i in seq_len(iter)) {
    if (verbose && i %in% progress_v) {
      cat("----------------\n")
      cat("Iter:", i, "\n")
      if (i > 1) {
        cat("Acceptance rate: ", acceptance_count / (i - 1), "\n")
      }
    }

    # Split X by label
    for (j in seq_len(k)) {
      Xv[[j]] <- X[Z == j]
      nv[j] <- length(Xv[[j]])
    }

    # Update pi (mixture weights)
    updating_dirichlet_shape <- nv + delta
    pi <- as.numeric(rdirichlet(1, updating_dirichlet_shape))

    # Update mu (location parameters, order constrained)
    for (j in seq_len(k)) {
      new_mu <- rnorm(1, mean = mu[j], sd = MHsd_mu[j])
      acceptance_rate <- mu_acceptance(
        mu[j], new_mu, Xv[[j]], nv[j], tau[j], xi, kappa
      )
      if (runif(1) < acceptance_rate) {
        mu[j] <- new_mu
        acceptance_count[j] <- acceptance_count[j] + 1
      }
    }

    sorted_indices <- order(mu)
    Z <- sorted_indices[Z]
    mu <- mu[sorted_indices]
    tau <- tau[sorted_indices]
    Xv <- Xv[sorted_indices]
    nv <- nv[sorted_indices]

    # Update tau (scale parameters)
    for (j in seq_len(k)) {
      new_tau <- rtruncnorm(1, a = left_boundary_TN, b = Inf, mean = tau[j], sd = MHsd_tau[j])
      acceptance_rate <- tau_acceptance(
        tau[j], new_tau, Xv[[j]], nv[j], mu[j], a, b, MHsd_tau[j]
      )
      if (runif(1) < acceptance_rate) {
        tau[j] <- new_tau
        acceptance_count[k + j] <- acceptance_count[k + j] + 1
      }
    }

    # Update Z (labels)
    for (t in seq_len(n)) {
      probs <- sapply(seq_len(k), function(j) {
        pi[j] * max(threshold, dgumbel(X[t], mu[j], 1 / tau[j]))
      })
      Z[t] <- sample(seq_len(k), size = 1, prob = probs)
    }

    # Split X by label
    for (j in seq_len(k)) {
      Xv[[j]] <- X[Z == j]
      nv[j] <- length(Xv[[j]])
    }

    # Store draws
    parameter[i, ] <- c(mu, tau, pi)
  }

  list(
    parameter = parameter,
    acc_rate = acceptance_count / iter
  )
}
