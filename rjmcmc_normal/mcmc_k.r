source(file.path(dir_path, "rjmcmc_common/utilities.r"))

mcmc_mixture_normal_k <- function(iter, X, # number of iterations, data vector
                                  xi, kappa, alpha,
                                  beta, delta, # prior parameters
                                  mu, tau, pi, Z, # initial values
                                  input_k = 3) {
  n <- length(X)
  k <- input_k

  # Matrix to store samples: mu1, tau1, mu2, tau2, mu3, tau3, pi1, pi2, pi3
  parameter <- matrix(0, iter, 3 * k)

  Xv <- vector("list", k) # List to hold data assigned to each component
  nv <- rep(0, k) # Counts for each component

  progress_v <- floor(iter * c(1:10) * 0.1)

  for (t in 1:iter) {
    if (t %in% progress_v) print(t)

    # Update Xv (grouped data) and nv (counts)
    for (j in 1:k) {
      Xv[[j]] <- X[Z == j]
      nv[j] <- length(Xv[[j]])
    }

    # ----- Update mixture weights (pi) -----
    updating_dirichlet_shape <- nv + delta
    pi <- rdirichlet(1, updating_dirichlet_shape)[1, ]

    # ----- Update component means (mu), enforce increasing order -----
    for (j in 1:k) {
      updating_tau <- nv[j] * tau[j] + kappa
      updating_mu <- (tau[j] * sum(Xv[[j]]) + xi * kappa) / updating_tau
      mu[j] <- rnorm(1, updating_mu, sqrt(1 / updating_tau))
    }

    sorted_indices <- order(mu)
    Z <- sorted_indices[Z]
    mu <- mu[sorted_indices]
    tau <- tau[sorted_indices]
    Xv <- Xv[sorted_indices]
    nv <- nv[sorted_indices]

    # ----- Update component precisions (tau) -----
    for (j in 1:k) {
      updating_alpha <- nv[j] / 2 + alpha
      updating_beta <- sum((Xv[[j]] - mu[j])^2) / 2 + beta
      tau[j] <- rgamma(1, updating_alpha, updating_beta)
    }

    # ----- Update component assignments (Z) -----
    for (i in 1:n) {
      pv <- sapply(1:k, function(j) {
        pi[j] * max(threshold, dnorm(X[i], mu[j], 1 / sqrt(tau[j])))
      })
      pv <- pv / sum(pv)
      Z[i] <- sample(1:k, size = 1, prob = pv)
    }

    # Update Xv (grouped data) and nv (counts)
    for (j in 1:k) {
      Xv[[j]] <- X[Z == j]
      nv[j] <- length(Xv[[j]])
    }

    # Store current sample
    parameter[t, ] <- c(mu, tau, pi)
  }
  return(list(parameter))
}

# mcmc_mixture_normal_k <- function(iter, X,      # number of iterations, data vector
#                                   xi, kappa, alpha, beta, delta, # prior parameters
#                                   mu, tau, pi, Z,               # initial values
#                                   input_k = 3) {

#   n <- length(X)
#   k <- input_k

#   # Matrix to store samples: mu1, tau1, mu2, tau2, mu3, tau3, pi1, pi2, pi3
#   parameter <- matrix(0, iter, 3 * k)

#   Xv <- vector("list", k)   # List to hold data assigned to each component
#   nv <- rep(0, k)           # Counts for each component

#   progress_v <- floor(iter * c(1:10) * 0.1)

#   for (t in 1:iter) {
#     if (t %in% progress_v) print(t)

#     # ----- Update component assignments (Z) -----
#     for (i in 1:n) {
#       pv <- sapply(1:k, function(j) {
#         pi[j] * max(threshold_normal, dnorm(X[i], mu[j], 1 / sqrt(tau[j])))
#       })
#       pv <- pv / sum(pv)
#       Z[i] <- sample(1:k, size = 1, prob = pv)
#     }

#     # Update Xv (grouped data) and nv (counts)
#     for (j in 1:k) {
#       Xv[[j]] <- X[Z == j]
#       nv[j] <- length(Xv[[j]])
#     }

#     # ----- Update component means (mu), enforce increasing order -----
#     while_count <- 0
#     old_mu <- mu
#     while (TRUE) {
#       while_count <- while_count + 1
#       mu <- old_mu
#       for (j in 1:k) {
#         updating_tau <- nv[j] * tau[j] + kappa
#         updating_mu <- (tau[j] * sum(Xv[[j]]) + xi * kappa) / updating_tau
#         mu[j] <- rnorm(1, updating_mu, sqrt(1 / updating_tau))
#       }
#       if (k == 1 || all(diff(mu) > 0)) break
#       if (while_count > 1000) {
#         sorted_indices <- order(mu)
#         Z <- sorted_indices[Z]
#         mu <- mu[sorted_indices]
#         tau <- tau[sorted_indices]
#         Xv <- Xv[sorted_indices]
#         nv <- nv[sorted_indices]
#         cat("while count break !!!!!!\n")
#         break
#       }
#     }

#     # ----- Update component precisions (tau) -----
#     for (j in 1:k) {
#       updating_alpha <- nv[j] / 2 + alpha
#       updating_beta <- sum((Xv[[j]] - mu[j])^2) / 2 + beta
#       tau[j] <- rgamma(1, updating_alpha, updating_beta)
#     }

#     # ----- Update mixture weights (pi) -----
#     updating_dirichlet_shape <- nv + delta
#     pi <- rdirichlet(1, updating_dirichlet_shape)[1, ]

#     # Store current sample
#     parameter[t, ] <- c(mu, tau, pi)
#   }
#   return(list(parameter))
# }
