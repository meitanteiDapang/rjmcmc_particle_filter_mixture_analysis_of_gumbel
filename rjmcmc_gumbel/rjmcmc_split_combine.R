source(file.path(dir_path, "rjmcmc_gumbel", "rjmcmc_acceptance.r"))

split_combine <- function(lst) {
  k <- lst$k
  birth_pr <- lst$birth_pr
  death_pr <- lst$death_pr

  # Propose move type: split, combine, or keep (no change)
  move <-
    sample(c(1, -1, 0), 1,
      prob =
        c(birth_pr[k], death_pr[k], 1 - birth_pr[k] - death_pr[k])
    )
  if (move == 0) {
    return(lst) # No move, return input
  }

  new_k <- k + move
  is_split <- (move == 1)

  # Extract params
  pi <- lst$pi
  mu <- lst$mu
  tau <- lst$tau
  Z <- lst$Z
  X <- lst$X
  n <- lst$n
  Xv <- lst$Xv
  nv <- lst$nv

  # Pre-allocate new params
  new_mu <- rep(0, new_k)
  new_tau <- rep(0, new_k)
  new_pi <- rep(0, new_k)
  new_Z <- lst$Z
  new_Xv <- vector("list", new_k)
  new_nv <- rep(0, new_k)

  if (is_split) {
    # ----- Split Move -----
    if (exists("split_count")) {
      split_count <<- split_count + 1
    }

    # Choose a component to split
    j.star <- sample(c(1:k), size = 1)

    j1 <- j.star
    j2 <- j1 + 1

    # Copy unchanged components
    if (j1 > 1) {
      new_pi[1:(j1 - 1)] <- pi[1:(j.star - 1)]
      new_mu[1:(j1 - 1)] <- mu[1:(j.star - 1)]
      new_tau[1:(j1 - 1)] <- tau[1:(j.star - 1)]
    }
    if (j2 < new_k) {
      new_pi[(j2 + 1):new_k] <- pi[(j.star + 1):k]
      new_mu[(j2 + 1):new_k] <- mu[(j.star + 1):k]
      new_tau[(j2 + 1):new_k] <- tau[(j.star + 1):k]
    }

    # Propose new parameters for the split components
    u1 <- rbeta(1, 2, 2)
    u2 <- rbeta(1, 2, 2)
    u3 <- rbeta(1, 1, 1)

    new_pi[j1] <- pi[j.star] * u1
    new_pi[j2] <- pi[j.star] * (1 - u1)

    new_mu[j1] <- mu[j.star] -
      u2 * (1 / sqrt(tau[j.star])) * sqrt(new_pi[j2] / new_pi[j1])
    new_mu[j2] <- mu[j.star] +
      u2 * (1 / sqrt(tau[j.star])) * sqrt(new_pi[j1] / new_pi[j2])

    new_tau[j1] <- new_pi[j1] / pi[j.star] *
      tau[j.star] * (1 / (u3 * (1 - u2^2)))
    new_tau[j2] <- new_pi[j2] / pi[j.star] *
      tau[j.star] * (1 / ((1 - u3) * (1 - u2^2)))

    is_strict_inc <- function(x, tol = 0) {
      if (length(x) <= 1) {
        return(TRUE)
      }
      return(all(diff(x) > tol))
    }

    if (!is_strict_inc(new_mu)) {
      return(lst) # Reject move if new mu is not strictly increasing
    }

    # Reassign data points for the split group
    for (t in seq_len(n)) {
      if (Z[t] == j.star) {
        new_Z[t] <- sample(c(j1, j2),
          size = 1, prob = c(new_pi[j1], new_pi[j2])
        )
      } else if (Z[t] > j.star) {
        new_Z[t] <- Z[t] + 1
      } else {
        new_Z[t] <- Z[t]
      }
    }
    # Update clusters
    for (j in seq_len(new_k)) {
      new_Xv[[j]] <- X[new_Z == j]
      new_nv[j] <- length(new_Xv[[j]])
    }
  } else {
    # ----- Combine Move -----
    if (exists("combine_count")) {
      combine_count <<- combine_count + 1
    }

    # Randomly select two adjacent components to merge
    j.star <- sample(c(1:new_k), size = 1)

    # j.start: after combine, j1 and j2 will be merged into j.star
    j1 <- j.star
    j2 <- j1 + 1

    # Copy parameters not affected by combine
    if (j1 > 1) {
      new_pi[1:(j.star - 1)] <- pi[1:(j1 - 1)]
      new_mu[1:(j.star - 1)] <- mu[1:(j1 - 1)]
      new_tau[1:(j.star - 1)] <- tau[1:(j1 - 1)]
    }
    if (j2 < k) {
      new_pi[(j.star + 1):new_k] <- pi[(j2 + 1):k]
      new_mu[(j.star + 1):new_k] <- mu[(j2 + 1):k]
      new_tau[(j.star + 1):new_k] <- tau[(j2 + 1):k]
    }

    # Compute merged parameters
    new_pi[j.star] <- pi[j1] + pi[j2]
    new_mu[j.star] <- (pi[j1] * mu[j1] + pi[j2] * mu[j2]) / new_pi[j.star]

    # variance = 1/tau
    new_variance <- (pi[j1] * ((1 / tau[j1]) + mu[j1]^2) +
      pi[j2] * ((1 / tau[j2]) + mu[j2]^2)) /
      new_pi[j.star] -
      new_mu[j.star]^2

    new_tau[j.star] <- max(1 / new_variance, tau_threshold)

    # Compute reverse proposal variables
    u1 <- pi[j1] / new_pi[j.star]
    sqrt_value <- (pi[j1] * mu[j1]^2 + pi[j2] * mu[j2]^2) /
      new_pi[j.star] - new_mu[j.star]^2
    u2 <- ifelse(sqrt_value > 0, sqrt(new_tau[j.star]) *
      sqrt_value, tau_threshold)
    u3 <- min(
      (pi[j1] * new_tau[j.star]) / ((1 - u2^2) * tau[j1] * new_pi[j.star]),
      0.9999
    )

    # Relabel clusters for combine
    for (t in c(1:n)) {
      if (new_Z[t] > j.star) {
        new_Z[t] <- new_Z[t] - 1
      }
    }
    # Update clusters
    for (j in seq_len(new_k)) {
      new_Xv[[j]] <- X[new_Z == j]
      new_nv[j] <- length(new_Xv[[j]])
    }
  }

  # Acceptance probability
  acceptance_r <- split_combine_acceptance(
    is_split, X, ifelse(is_split, k, new_k), Z, mu, tau, Xv, nv, pi,
    new_Z, new_mu, new_tau, new_Xv, new_nv, new_pi,
    lst$xi, lst$kappa, lst$a, lst$b, lst$lambda, lst$delta,
    u1, u2, u3, j1, j2, j.star,
    birth_pr, death_pr
  )

  # Accept or reject
  if (runif(1) < acceptance_r) {
    lst$Z <- new_Z
    lst$k <- new_k
    lst$pi <- new_pi
    lst$mu <- new_mu
    lst$tau <- new_tau
    lst$Xv <- new_Xv
    lst$nv <- new_nv
    if (exists("split_accept_count", envir = .GlobalEnv)) {
      if (is_split) {
        split_accept_count <<- split_accept_count + 1
      } else {
        combine_accept_count <<- combine_accept_count + 1
      }
    }
  }
  return(lst)
}
