source(file.path(dir_path, "rjmcmc_normal/split_combine_steps.r"))
source(file.path(dir_path, "rjmcmc_normal/birth_death_steps.r"))





rjmcmc_mixture_normal <- function(
    iter, X, max.K, # iteration count, data
    xi, kappa, a, b, lambda, delta, # prior parameters
    k, mu, tau, pi, Z, # initial value
    birth_pr, death_pr
    ) {
  threshold <- 1e-20
  K <- max.K
  parameter <- matrix(0, iter, 1 + 3 * K) # Store k, mu, tau, pi
  n <- length(X)
  all_Z <- matrix(1, iter, n)

  progress_v <- floor(iter * c(1:100) * 0.01)

  for (i in 1:iter) {
    if (i %in% progress_v) {
      cat("Iter:", i, " k: ", k, "\n")
    }

    # Update latent labels
    Xv <- vector("list", k)
    nv <- rep(0, k)
    for (t in c(1:n)) {
      pv <- sapply(c(1:k), function(j) {
        pi[j] * max(threshold, dnorm(X[t], mu[j], 1 / sqrt(tau[j])))
      })
      Z[t] <- sample(c(1:k), size = 1, prob = pv)
    }

    # Group observations and count per component
    for (j in c(1:k)) {
      Xv[[j]] <- X[Z == j]
      nv[j] <- length(Xv[[j]])
    }

    # Update mixture weights
    updating_dirichlet_shape <- nv + delta
    pi <- rdirichlet(1, updating_dirichlet_shape)[1, ]

    # Update means, ensure sorted order to mitigate label switching
    for (j in c(1:k)) {
      updating_tau <- (nv[j] * tau[j] + kappa)
      updating_mu <- (tau[j] * sum(Xv[[j]]) + xi * kappa) / updating_tau
      mu[j] <- rnorm(1, updating_mu, sqrt(1 / updating_tau))
    }
    sorted_indices <- order(mu)
    Z <- sorted_indices[Z]
    mu <- mu[sorted_indices]
    tau <- tau[sorted_indices]
    Xv <- Xv[sorted_indices]
    nv <- nv[sorted_indices]


    # Update precisions
    for (j in c(1:k)) {
      updating_a <- nv[j] / 2 + a
      updating_beta <- sum((Xv[[j]] - mu[j])^2) / 2 + b
      tau[j] <- max(rgamma(1, updating_a, updating_beta), tau_threshold)
    }

    # Split/combine move
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

    # Ensure means remain sorted after RJMCMC move
    sorted_indices <- order(mu)
    Z <- sorted_indices[Z]
    mu <- mu[sorted_indices]
    tau <- tau[sorted_indices]
    Xv <- Xv[sorted_indices]
    nv <- nv[sorted_indices]

    # Birth/death move
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

    # Ensure means remain sorted after RJMCMC move
    sorted_indices <- order(mu)
    Z <- sorted_indices[Z]
    mu <- mu[sorted_indices]
    tau <- tau[sorted_indices]
    Xv <- Xv[sorted_indices]
    nv <- nv[sorted_indices]



    # Store draws
    all_Z[i, ] <- Z
    parameter[i, ] <- c(k, mu, tau, pi, rep(0, 3 * (K - k)))
  }
  return(list(parameter, all_Z))
}







# rjmcmc_mixture_normal_random_b <- function(
#     iter, X, max.K, # iteration count, data
#     xi, kappa, a, g, h, lambda, delta, # prior parameters
#     k, mu, tau, pi, Z, b, # initial value
#     birth_pr, death_pr
#     ) {
#   threshold <- 1e-20
#   K <- max.K
#   parameter <- matrix(0, iter, 1 + 3 * K) # Store k, mu, tau, pi
#   n <- length(X)
#   all_Z <- matrix(1, iter, n)

#   progress_v <- floor(iter * c(1:100) * 0.01)

#   for (i in 1:iter) {
#     if (i %in% progress_v) {
#       cat("Iter:", i, " k: ", k, "\n")
#     }

#     # Update latent labels
#     Xv <- vector("list", k)
#     nv <- rep(0, k)
#     for (t in c(1:n)) {
#       pv <- sapply(c(1:k), function(j) {
#         pi[j] * max(threshold, dnorm(X[t], mu[j], 1 / sqrt(tau[j])))
#       })
#       Z[t] <- sample(c(1:k), size = 1, prob = pv)
#     }

#     # Group observations and count per component
#     for (j in c(1:k)) {
#       Xv[[j]] <- X[Z == j]
#       nv[j] <- length(Xv[[j]])
#     }

#     # Update mixture weights
#     updating_dirichlet_shape <- nv + delta
#     pi <- rdirichlet(1, updating_dirichlet_shape)[1, ]

#     # Update means, ensure sorted order to mitigate label switching
#     for (j in c(1:k)) {
#       updating_tau <- (nv[j] * tau[j] + kappa)
#       updating_mu <- (tau[j] * sum(Xv[[j]]) + xi * kappa) / updating_tau
#       mu[j] <- rnorm(1, updating_mu, sqrt(1 / updating_tau))
#     }
#     sorted_indices <- order(mu)
#     Z <- sorted_indices[Z]
#     mu <- mu[sorted_indices]
#     tau <- tau[sorted_indices]
#     Xv <- Xv[sorted_indices]
#     nv <- nv[sorted_indices]


#     # Update precisions
#     for (j in c(1:k)) {
#       updating_a <- nv[j] / 2 + a
#       updating_beta <- sum((Xv[[j]] - mu[j])^2) / 2 + b
#       tau[j] <- max(rgamma(1, updating_a, updating_beta), tau_threshold)
#     }

#     # update b
#     b <- rgamma(1, shape = g + k * a, rate =  (h + sum(tau)))


#     # Split/combine move
#     lst <- generate_lst(
#       X, Z, K, k, pi, mu, tau, n, Xv, nv,
#       xi, kappa, a, b, lambda, delta,
#       birth_pr, death_pr
#     )
#     lst <- split_combine(lst)
#     Z <- lst$Z
#     k <- lst$k
#     pi <- lst$pi
#     mu <- lst$mu
#     tau <- lst$tau
#     Xv <- lst$Xv
#     nv <- lst$nv

#     # Ensure means remain sorted after RJMCMC move
#     sorted_indices <- order(mu)
#     Z <- sorted_indices[Z]
#     mu <- mu[sorted_indices]
#     tau <- tau[sorted_indices]
#     Xv <- Xv[sorted_indices]
#     nv <- nv[sorted_indices]

#     # Birth/death move
#     lst <- generate_lst(
#       X, Z, K, k, pi, mu, tau, n, Xv, nv,
#       xi, kappa, a, b, lambda, delta,
#       birth_pr, death_pr
#     )
#     lst <- birth_death(lst)
#     Z <- lst$Z
#     k <- lst$k
#     pi <- lst$pi
#     mu <- lst$mu
#     tau <- lst$tau
#     Xv <- lst$Xv
#     nv <- lst$nv

#     # Ensure means remain sorted after RJMCMC move
#     sorted_indices <- order(mu)
#     Z <- sorted_indices[Z]
#     mu <- mu[sorted_indices]
#     tau <- tau[sorted_indices]
#     Xv <- Xv[sorted_indices]
#     nv <- nv[sorted_indices]



#     # Store draws
#     all_Z[i, ] <- Z
#     parameter[i, ] <- c(k, mu, tau, pi, rep(0, 3 * (K - k)))
#   }
#   return(list(parameter, all_Z))
# }




rjmcmc_mixture_normal_random_b <- function(
    iter, X, max.K, # iteration count, data
    xi, kappa, a, g, h, lambda, delta, # prior parameters
    k, mu, tau, pi, Z, b, # initial value
    birth_pr, death_pr
    ) {
  threshold <- 1e-20
  K <- max.K
  parameter <- matrix(0, iter, 1 + 3 * K) # Store k, mu, tau, pi
  n <- length(X)
  all_Z <- matrix(1, iter, n)

  progress_v <- floor(iter * c(1:10) * 0.1)

  for (i in 1:iter) {
    if (i %in% progress_v) {
      cat("Iter:", i, " k: ", k, "\n")
    }
    Xv <- vector("list", k)
    nv <- rep(0, k)


    # Group observations and count per component
    for (j in c(1:k)) {
      Xv[[j]] <- X[Z == j]
      nv[j] <- length(Xv[[j]])
    }

    # Update mixture weights
    updating_dirichlet_shape <- nv + delta
    pi <- rdirichlet(1, updating_dirichlet_shape)[1, ]


    # Update means, ensure sorted order to mitigate label switching
    for (j in c(1:k)) {
      updating_tau <- (nv[j] * tau[j] + kappa)
      updating_mu <- (tau[j] * sum(Xv[[j]]) + xi * kappa) / updating_tau

      mu[j] <- rnorm(1, updating_mu, sqrt(1 / updating_tau))
    }
    sorted_indices <- order(mu)
    Z <- sorted_indices[Z]
    mu <- mu[sorted_indices]
    tau <- tau[sorted_indices]
    Xv <- Xv[sorted_indices]
    nv <- nv[sorted_indices]


    # Update precisions
    for (j in c(1:k)) {
      updating_a <- nv[j] / 2 + a
      updating_beta <- sum((Xv[[j]] - mu[j])^2) / 2 + b
      tau[j] <- max(rgamma(1, updating_a, updating_beta), tau_threshold)
    }

    # Update latent labels
    for (t in c(1:n)) {
      pv <- sapply(c(1:k), function(j) {
        pi[j] * max(threshold, dnorm(X[t], mu[j], 1 / sqrt(tau[j])))
      })

      Z[t] <- sample(c(1:k), size = 1, prob = pv)
    }

    # Group observations and count per component
    for (j in c(1:k)) {
      Xv[[j]] <- X[Z == j]
      nv[j] <- length(Xv[[j]])
    }
    
    # update b
    b <- rgamma(1, shape = g + k * a, rate =  (h + sum(tau)))

    # Split/combine move
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


    # Ensure means remain sorted after RJMCMC move
    sorted_indices <- order(mu)
    Z <- sorted_indices[Z]
    mu <- mu[sorted_indices]
    tau <- tau[sorted_indices]
    Xv <- Xv[sorted_indices]
    nv <- nv[sorted_indices]


    # Birth/death move
    lst <- generate_lst(
      X, Z, K, k, pi, mu, tau, n, Xv, nv,
      xi, kappa, a, b, lambda, delta,
      birth_pr, death_pr, i
    )

    lst <- birth_death(lst)
    Z <- lst$Z
    k <- lst$k
    pi <- lst$pi
    mu <- lst$mu
    tau <- lst$tau
    Xv <- lst$Xv
    nv <- lst$nv

    # Ensure means remain sorted after RJMCMC move
    sorted_indices <- order(mu)
    Z <- sorted_indices[Z]
    mu <- mu[sorted_indices]
    tau <- tau[sorted_indices]
    Xv <- Xv[sorted_indices]
    nv <- nv[sorted_indices]


    # Store draws
    all_Z[i, ] <- Z
    parameter[i, ] <- c(k, mu, tau, pi, rep(0, 3 * (K - k)))
  }
  return(list(parameter, all_Z))
}