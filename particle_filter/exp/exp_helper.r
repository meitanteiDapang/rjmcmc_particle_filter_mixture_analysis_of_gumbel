library(MCMCpack)
source(file.path(dir_path, "particle_filter/fearnhead_utilities.r"))



fearnhead_gibbs_exp <- function(y, z, a, b, Niter = 1500, burn.in = 500) {
  # Number of clusters
  k <- max(z)

  # Observations per cluster
  n <- c(table(z))

  # Storage for samples
  mon_lambda <- array(dim = c(k, Niter))

  # Initial estimate for lambda
  lambda <- rep(1 / mean(y), k)  # Initial guess for rate

  for (iter in 1:Niter) {
    # Sum of observations per cluster
    y_sum <- tapply(y, z, sum)

    # Gibbs update for lambda_j ~ Gamma(a + n_j, b + sum(y_j))
    mon_lambda[, iter] <- lambda <- rgamma(
      k,
      shape = a + n,
      rate = b + y_sum
    )
  }

  # Return posterior samples after burn-in
  return(
    list(
      lambda = mon_lambda[, (burn.in + 1):Niter, drop = FALSE]
    )
  )
}


get_estimation_in_grid_exp <- function(x.grid, uni_particles) {
  n_samples <- 1000
  N_uni_particles <- length(uni_particles)

  # For estimation: row: particles number, column: grid
  D.est <- array(dim = c(N_uni_particles, length(x.grid)))

  for (p in 1:N_uni_particles) {
    print(p)

    # consider one particle in one time
    zz <- uni_particles[[p]]

    # the parameter vector, mm$mu, mm$tau
    mm <- fearnhead_gibbs_exp(
      y = y, z = zz, a = a, b = b,
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
          dexp(
            x.grid, mm$lambda[j, iter]
          ) * mean(zz == j)
      }
    }

    # finally averaging
    D.est[p, ] <- apply(ff, 1, mean)
  }
  return(D.est)
}



get_posterior_parameters_exp <- function(particles, weights, k, samples_to_get) {
  n_samples <- 1000

  mm_list <- list()
  pi_list <- list()

  for (i in seq_along(particles)) {
    cat(i, "/", length(particles), "\n")
    zz <- particles[[i]]

    # Cluster proportions
    cluster_counts <- table(factor(zz, levels = 1:k))
    pi_est <- as.numeric(cluster_counts) 
    pi_list[[i]] <- rep(list(pi_est), n_samples)

    # Gibbs sampling
    mm <- fearnhead_gibbs_exp(
      y = y, z = zz, a = a, b = b,
      Niter = 1500, burn.in = 1500 - n_samples
    )
    mm_list[[i]] <- mm
  }

  # Output matrix: (lambda_1,...,lambda_k, pi_1,...,pi_k)
  output_matrix <- matrix(NA, nrow = samples_to_get, ncol = 2 * k)

  for (i in seq_len(samples_to_get)) {
    p <- sample(seq_along(particles), size = 1, prob = weights)
    iter_idx <- sample(1:n_samples, size = 1)

    lambda  <- mm_list[[p]]$lambda[, iter_idx]
    pi  <- rdirichlet(1, 1 + pi_list[[p]][[iter_idx]])

    output_matrix[i, ] <- c(lambda, pi)
  }

  return(output_matrix)
}