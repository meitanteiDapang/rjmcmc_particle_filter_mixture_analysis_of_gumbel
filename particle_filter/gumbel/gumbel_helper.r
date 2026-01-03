library(MCMCpack)
source(file.path(dir_path, "particle_filter/fearnhead_utilities.r"))
source(file.path(dir_path, "particle_filter/fc_resample.r"))
source(file.path(dir_path, "particle_filter/pf_utilities.r"))
source(file.path(dir_path, "rjmcmc_gumbel/single_gumbel.r"))







get_posterior_parameters_gumbel <- function(particles, weights, k, samples_to_get, y,
                                            eta, phi, a, b, MHsd_mu_l = rep(0.5, k), MHsd_tau_l = rep(0.01, k),
                                            total_iter = 1500, burn_in = 500) {
  n_samples <- total_iter - burn_in

  mm_list <- list()
  pi_list <- list()

  for (i in seq_along(particles)) {
    zz <- particles[[i]]
    cat(i, "/", length(particles), "\n")

    # Cluster proportions
    cluster_counts <- table(factor(zz, levels = 1:k))
    pi_list[[i]] <- rep(list(as.numeric(cluster_counts)), n_samples)

    # Prepare matrices to store posterior samples for each cluster
    mu_matrix <- matrix(NA, nrow = k, ncol = n_samples)
    tau_matrix <- matrix(NA, nrow = k, ncol = n_samples)

    for (j in 1:k) {
      y_j <- y[zz == j]
      # cat("y_j: ", y_j, "\n")
      if (length(y_j) == 0) {
        next
      }

      # Initial values for mu and tau
      gamma_Euler <- 0.5772
      var_y_j <- var(y_j)
      if (is.na(var_y_j) || var_y_j == 0) {
        var_y_j <- 1
      }
      beta_hat <- sqrt(6 / real_pi^2 * var_y_j)
      mu_init <- mean(y_j) - beta_hat * gamma_Euler
      tau_init <- 1 / beta_hat
      # cat("mu_init: ", mu_init, " tau_init: ", tau_init, "\n")
      draws <- mcmc_single_gumbel(
        iter = total_iter,
        X = y_j,
        xi = eta, kappa = 1 / phi, # kappa being the precivsion, and phi being the standard deviation
        a = a, b = b,
        mu = mu_init, tau = tau_init,
        MHsd_mu = MHsd_mu_l[j],
        MHsd_tau = MHsd_tau_l[j],
        verbose = FALSE
      )

      mu_matrix[j, ] <- draws$parameter[(burn_in + 1):total_iter, 1]
      tau_matrix[j, ] <- draws$parameter[(burn_in + 1):total_iter, 2]
    }

    mm_list[[i]] <- list(mu = mu_matrix, tau = tau_matrix)
  }

  # Final matrix: [mu_1,...,mu_k, tau_1,...,tau_k, pi_1,...,pi_k]
  output_matrix <- matrix(NA, nrow = samples_to_get, ncol = 3 * k)

  for (i in seq_len(samples_to_get)) {
    p <- sample(seq_along(particles), size = 1, prob = weights)
    iter_idx <- sample(1:n_samples, size = 1)

    mu_sample <- mm_list[[p]]$mu[, iter_idx]
    tau_sample <- mm_list[[p]]$tau[, iter_idx]
    pi_sample <- rdirichlet(1, 1 + pi_list[[p]][[iter_idx]])

    output_matrix[i, ] <- c(mu_sample, tau_sample, pi_sample)
  }

  return(output_matrix)
}
