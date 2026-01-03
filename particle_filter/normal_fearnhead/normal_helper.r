library(evd)
library(MCMCpack)

fearnhead_gibbs_normal <- function(y, z, xi, phi, a, b, Niter = 1500, burn.in = 500) {
  # number of clusters
  k <- max(z)

  # number of observations by group
  n <- c(table(z))

  # setting the monitors
  mon_mu <- mon_tau <- array(dim = c(k, Niter))

  # initialising
  mu <- tapply(y, z, mean)
  tau <- 1 / tapply(y, z, var)
  tau[is.na(tau)] <- mean(tau, na.rm = TRUE) # for groups with membership of 1

  for (iter in 1:Niter) {
    # means
    mu_mean <- (tapply(y, z, sum) * phi + xi) / (n * phi + 1)
    mon_mu[, iter] <- mu <- rnorm(k, mu_mean,
      sd =
        sqrt(phi / (tau + tau * n * phi))
    )
    mon_tau[, iter] <- tau <-
      rgamma(
        k, a + n / 2,
        b + 1 / 2 * n * (tapply(y, z, varx) +
          (tapply(y, z, mean) - xi)^2 / (1 + n * phi))
      )
  }
  return(
    list(
      mu = mon_mu[, (burn.in + 1):Niter, drop=FALSE],
      tau = mon_tau[, (burn.in + 1):Niter, drop=FALSE]
    )
  )
}

get_estimation_in_grid <- function(x.grid, uni_particles) {
  n_samples <- 1000
  N_uni_particles <- length(uni_particles)

  # For estimation: row: particles number, column: grid
  D.est <- array(dim = c(N_uni_particles, length(x.grid)))

  for (p in 1:N_uni_particles) {
    print(p)

    # consider one particle in one time
    zz <- uni_particles[[p]]

    # the parameter vector, mm$mu, mm$tau
    mm <- fearnhead_gibbs_normal(
      y = y, z = zz,
      xi = xi, kappa = kappa, a = a, b = b,
      Niter = 1500, burn.in = 1500 - n_samples
    )

    # draw n_samples of parameters, then estimates them on
    # all grid points.
    ff <- array(0, dim = c(length(x.grid), n_samples))

    k <- max(zz)

    for (iter in 1:n_samples) {
      # normal case: multiple clusters
      for (j in 1:k) {
        ff[, iter] <- ff[, iter] +
          dnorm(
            x.grid, mm$mu[j, iter],
            1 / sqrt(mm$tau[j, iter])
          ) * mean(zz == j)
      }
    }

    # finally averaging
    D.est[p, ] <- apply(ff, 1, mean)
  }
  return(D.est)
}


get_posterior_parameters_normal <- function(particles, weights, k, samples_to_get,
xi, kappa, a, b,y) {
  n_samples <- 1000

  mm_list <- list()
  pi_list <- list()

  for (i in seq_along(particles)) {
    zz <- particles[[i]]
    cat(i, "/", length(particles), "\n")
    # Cluster proportions
    cluster_counts <- table(factor(zz, levels = 1:k))
    pi_est <- as.numeric(cluster_counts)
    pi_list[[i]] <- rep(list(pi_est), n_samples)

    # Gibbs sampling
    mm <- fearnhead_gibbs_normal(
      y = y, z = zz,
      xi = xi, phi = 1/kappa, a = a, b = b,
      Niter = 1500, burn.in = 1500 - n_samples
    )
    mm_list[[i]] <- mm
  }

  # Output matrix: (mu_1,...,mu_k, tau_1,...,tau_k, pi_1,...,pi_k)
  output_matrix <- matrix(NA, nrow = samples_to_get, ncol = 3 * k)

  for (i in seq_len(samples_to_get)) {
    p <- sample(seq_along(particles), size = 1, prob = weights)
    iter_idx <- sample(1:n_samples, size = 1)

    mu  <- mm_list[[p]]$mu[, iter_idx]
    tau <- mm_list[[p]]$tau[, iter_idx]
    pi  <- rdirichlet(1, 1 + pi_list[[p]][[iter_idx]])

    output_matrix[i, ] <- c(mu, tau, pi)
  }

  return(output_matrix)
}
