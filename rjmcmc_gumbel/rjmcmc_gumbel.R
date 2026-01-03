source(file.path(dir_path, "rjmcmc_common", "utilities.r"))
source(file.path(dir_path, "rjmcmc_gumbel", "tau_mu_acceptance.r"))
source(file.path(dir_path, "rjmcmc_gumbel", "rjmcmc_acceptance.r"))
source(file.path(dir_path, "rjmcmc_gumbel", "rjmcmc_split_combine.r"))
source(file.path(dir_path, "rjmcmc_gumbel", "rjmcmc_birth_death.r"))


rjmcmc_gumbel <- function(iter, X, k, # iteration count, data
                          max.K, min.K, birth_pr, death_pr,
                          xi, kappa, a, b, delta, lambda, # prior parameters
                          mu, tau, pi, Z, # initial value
                          MHsd_mu_list, MHsd_tau_list, left_boundary_TN,
                          verbose = TRUE) {
  n <- length(X)
  K <- max.K # Maximum allowed components

  # Store: k, mu, tau,pi (with zero padding for unused components)
  parameter <- matrix(0, iter, 1 + 3 * K)
  all_Z <- matrix(1, iter, n)

  acceptance_count_mu_tau <- rep(0, 2)
  whole_candidate_count <- 0
  Xv <- vector("list", k)
  nv <- rep(0, k)

  progress_v <- floor(iter * seq_len(100) / 100)

  for (i in seq_len(iter)) {
    if (verbose && i %in% progress_v) {
      cat("Iter:", i, "k:", k, "\n")
    }

    whole_candidate_count <- whole_candidate_count + k
    MHsd_mu <- MHsd_mu_list[[k]]
    MHsd_tau <- MHsd_tau_list[[k]]


    # Split X by label
    for (j in seq_len(k)) {
      Xv[[j]] <- X[Z == j]
      nv[j] <- length(Xv[[j]])
    }


    # Update mu (location parameters, order constrained)
    for (j in seq_len(k)) {
      new_mu <- rnorm(1, mean = mu[j], sd = MHsd_mu[j])
      acceptance_rate <- mu_acceptance(
        mu[j], new_mu, Xv[[j]], nv[j], tau[j], xi, kappa
      )
      if (runif(1) < acceptance_rate) {
        mu[j] <- new_mu
        acceptance_count_mu_tau[1] <- acceptance_count_mu_tau[1] + 1
      }
    }

    sorted_indices <- order(mu)
    Z <- sorted_indices[Z]
    mu <- mu[sorted_indices]
    tau <- tau[sorted_indices]
    Xv <- Xv[sorted_indices]
    nv <- nv[sorted_indices]

    # Update tau (scale)
    for (j in seq_len(k)) {
      new_tau <- rtruncnorm(1,
        a = left_boundary_TN, b = Inf, mean = tau[j], sd = MHsd_tau[j]
      )
      acceptance_rate <- tau_acceptance(
        tau[j],
        new_tau, Xv[[j]], nv[j], mu[j], a, b, MHsd_tau
      )
      if (runif(1) < acceptance_rate) {
        tau[j] <- new_tau
        acceptance_count_mu_tau[2] <- acceptance_count_mu_tau[2] + 1
      }
    }

    # Updatepi (mixture weights)
    pi <- as.numeric(rdirichlet(1, nv + delta))

    # Update Z (labels)
    for (t in seq_len(n)) {
      probs <- sapply(seq_len(k), function(j) {
        pi[j] * max(threshold, dgumbel(X[t], mu[j], 1 / tau[j]))
      })
      Z[t] <- sample(seq_len(k), 1, prob = probs)
    }

    # Split X by label
    for (j in seq_len(k)) {
      Xv[[j]] <- X[Z == j]
      nv[j] <- length(Xv[[j]])
    }


    # Split & combine step
    lst <- generate_lst(
      X, Z, K, k, pi, mu, tau, n, Xv, nv,
      xi, kappa, a, b, lambda, delta,
      birth_pr, death_pr
    )
    lst <- split_combine(lst)
    Z <- lst$Z
    k <- lst$k
    pi <- lst$pi
    mu <- lst$mu
    tau <- lst$tau
    Xv <- lst$Xv
    nv <- lst$nv

    # Sort by mu if needed
    if (any(diff(mu) <= 0)) {
      sorted_indices <- order(mu)
      Z <- sorted_indices[Z]
      mu <- mu[sorted_indices]
      tau <- tau[sorted_indices]
      Xv <- Xv[sorted_indices]
      nv <- nv[sorted_indices]
    }

    # Birth & death step
    lst <- generate_lst(
      X, Z, K, k, pi, mu, tau, n, Xv, nv,
      xi, kappa, a, b, lambda, delta,
      birth_pr, death_pr
    )
    lst <- birth_death(lst)
    Z <- lst$Z
    k <- lst$k
    pi <- lst$pi
    mu <- lst$mu
    tau <- lst$tau
    Xv <- lst$Xv
    nv <- lst$nv

    if (any(diff(mu) <= 0)) {
      sorted_indices <- order(mu)
      Z <- sorted_indices[Z]
      mu <- mu[sorted_indices]
      tau <- tau[sorted_indices]
      Xv <- Xv[sorted_indices]
      nv <- nv[sorted_indices]
    }

    # Store draws, padpiith 0 for unused parameters
    parameter[i, ] <- c(k, mu, tau, pi, rep(0, 3 * (K - k)))
    all_Z[i, ] <- Z
  }

  list(
    parameter = parameter,
    acc_rate_mu_tau = acceptance_count_mu_tau / whole_candidate_count,
    all_Z = all_Z
  )
}




rjmcmc_gumbel_random_b <- function(iter, X, k, # iteration count, data
                          max.K, min.K, birth_pr, death_pr,
                          xi, kappa, a, g, h, delta, lambda, # prior parameters
                          mu, tau, pi, Z, b, # initial value
                          MHsd_mu_list, MHsd_tau_list, left_boundary_TN,
                          verbose = TRUE) {
  n <- length(X)
  K <- max.K # Maximum allowed components

  # Store: k, mu, tau,pi (with zero padding for unused components)
  parameter <- matrix(0, iter, 1 + 3 * K)
  all_Z <- matrix(1, iter, n)

  acceptance_count_mu_tau <- rep(0, 2)
  whole_candidate_count <- 0
  Xv <- vector("list", k)
  nv <- rep(0, k)

  progress_v <- floor(iter * seq_len(100) / 100)

  for (i in seq_len(iter)) {
    if (verbose && i %in% progress_v) {
      cat("Iter:", i, "k:", k, "\n")
    }

    whole_candidate_count <- whole_candidate_count + k
    MHsd_mu <- MHsd_mu_list[[k]]
    MHsd_tau <- MHsd_tau_list[[k]]


    # Split X by label
    for (j in seq_len(k)) {
      Xv[[j]] <- X[Z == j]
      nv[j] <- length(Xv[[j]])
    }


    # Update mu (location parameters, order constrained)
    for (j in seq_len(k)) {
      new_mu <- rnorm(1, mean = mu[j], sd = MHsd_mu[j])
      acceptance_rate <- mu_acceptance(
        mu[j], new_mu, Xv[[j]], nv[j], tau[j], xi, kappa
      )
      if (runif(1) < acceptance_rate) {
        mu[j] <- new_mu
        acceptance_count_mu_tau[1] <- acceptance_count_mu_tau[1] + 1
      }
    }

    sorted_indices <- order(mu)
    Z <- sorted_indices[Z]
    mu <- mu[sorted_indices]
    tau <- tau[sorted_indices]
    Xv <- Xv[sorted_indices]
    nv <- nv[sorted_indices]

    # Update tau (scale)
    for (j in seq_len(k)) {
      new_tau <- rtruncnorm(1,
        a = left_boundary_TN, b = Inf, mean = tau[j], sd = MHsd_tau[j]
      )
      acceptance_rate <- tau_acceptance(
        tau[j],
        new_tau, Xv[[j]], nv[j], mu[j], a, b, MHsd_tau
      )
      if (runif(1) < acceptance_rate) {
        tau[j] <- new_tau
        acceptance_count_mu_tau[2] <- acceptance_count_mu_tau[2] + 1
      }
    }

    # Updatepi (mixture weights)
    pi <- as.numeric(rdirichlet(1, nv + delta))

    # Update Z (labels)
    for (t in seq_len(n)) {
      probs <- sapply(seq_len(k), function(j) {
        pi[j] * max(threshold, dgumbel(X[t], mu[j], 1 / tau[j]))
      })
      Z[t] <- sample(seq_len(k), 1, prob = probs)
    }

    # Split X by label
    for (j in seq_len(k)) {
      Xv[[j]] <- X[Z == j]
      nv[j] <- length(Xv[[j]])
    }

    # update b
    b <- rgamma(1, shape = g + k * a, rate =  (h + sum(tau)))

    # Split & combine step
    lst <- generate_lst(
      X, Z, K, k, pi, mu, tau, n, Xv, nv,
      xi, kappa, a, b, lambda, delta,
      birth_pr, death_pr
    )
    lst <- split_combine(lst)
    Z <- lst$Z
    k <- lst$k
    pi <- lst$pi
    mu <- lst$mu
    tau <- lst$tau
    Xv <- lst$Xv
    nv <- lst$nv

    # Sort by mu if needed
    if (any(diff(mu) <= 0)) {
      sorted_indices <- order(mu)
      Z <- sorted_indices[Z]
      mu <- mu[sorted_indices]
      tau <- tau[sorted_indices]
      Xv <- Xv[sorted_indices]
      nv <- nv[sorted_indices]
    }

    # Birth & death step
    lst <- generate_lst(
      X, Z, K, k, pi, mu, tau, n, Xv, nv,
      xi, kappa, a, b, lambda, delta,
      birth_pr, death_pr
    )
    lst <- birth_death(lst)
    Z <- lst$Z
    k <- lst$k
    pi <- lst$pi
    mu <- lst$mu
    tau <- lst$tau
    Xv <- lst$Xv
    nv <- lst$nv

    if (any(diff(mu) <= 0)) {
      sorted_indices <- order(mu)
      Z <- sorted_indices[Z]
      mu <- mu[sorted_indices]
      tau <- tau[sorted_indices]
      Xv <- Xv[sorted_indices]
      nv <- nv[sorted_indices]
    }

    # Store draws, padpiith 0 for unused parameters
    parameter[i, ] <- c(k, mu, tau, pi, rep(0, 3 * (K - k)))
    all_Z[i, ] <- Z
  }

  list(
    parameter = parameter,
    acc_rate_mu_tau = acceptance_count_mu_tau / whole_candidate_count,
    all_Z = all_Z
  )
}
