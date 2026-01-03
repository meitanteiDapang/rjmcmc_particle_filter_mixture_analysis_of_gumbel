library(evd)
source(file.path(dir_path, "particle_filter/pf_utilities.r"))
source(file.path(dir_path, "particle_filter/fearnhead_utilities.r"))
source(file.path(dir_path, "particle_filter/fc_resample.r"))
library(dplyr)


compute_log_L <- function(
    y_1t, z, a, b, xi, phi) {
  # Get basic cluster stats
  cluster_ids <- sort(unique(z))
  k <- length(cluster_ids)
  n <- length(y_1t)


  log_likelihood <- 0
  for (j in cluster_ids) {
    y_j <- y_1t[z == j]
    n_j <- length(y_j)
    y_bar_j <- mean(y_j)
    s2_j <- mean((y_j - y_bar_j)^2)


    # Posterior contribution from cluster j (log-scale)
    term1 <- a * log(b)
    term2 <- lgamma(a + n_j / 2) - lgamma(a)
    term3 <- -0.5 * log(n_j * phi + 1)

    quad <- s2_j + (y_bar_j - xi)^2 / (1 + n_j * phi)
    term4 <- -(a + n_j / 2) * log(b + 0.5 * n_j * quad)

    log_likelihood <- log_likelihood + term1 + term2 + term3 + term4
  }

  # Return unnormalized posterior (in log-scale and normal scale)
  return(list(
    log_likelihood = log_likelihood,
    likelihood = exp(log_likelihood)
  ))
}

pf_normal <- function(Y, N,
                      a, b, xi, phi, alpha,
                      initial_size = N,
                      is_add_term = FALSE,
                      verbose = TRUE) {
  # observation length
  n <- length(Y)

  # Initialisation
  Z <- rep(list(c(1)), initial_size)
  w <- rep(1 / initial_size, initial_size)
  k.curr <- rep(1, initial_size)
  last_L_v <- rep(0, initial_size)

  cat("First data: ", Y[1], "\n")
  progress_v <- floor(n * c(1:100) * 0.01)

  # Loop
  for (t in c(2:n)) {
    if (verbose) {
      if (t %in% progress_v) {
        cat("\nIter", t, ": ")
      }
    }

    y_1t <- Y[1:t]

    # Intialise empty new Z
    new_Z <- list()
    new_w <- c()
    new_Z_count <- 0
    new_L_v <- c()
    # generate all new Z and w
    # for each particle
    for (j in c(1:length(Z))) {
      # the i-th particle after t-1 th step
      this_Zj <- Z[[j]]

      # get the possible next allocation at t-th step
      possible_next_allocation <- 1 + k.curr[j]
      # for each possible allocation
      for (i in c(1:possible_next_allocation)) {
        new_Z_count <- new_Z_count + 1

        # generate the new particle for y_1t
        this_new_Z <- c(this_Zj, i)

        # store the new particle
        new_Z[[new_Z_count]] <- this_new_Z


        if (i == possible_next_allocation) {
          log_prior <- log((alpha) / (t - 1 + alpha))
        } else {
          # get nj for t-1 step particles
          ni <- as.numeric(table(this_Zj)[i])
          log_prior <- log((ni) / (t - 1 + alpha))
        }

        # compute the new log L for t-th step
        log_new_L <- compute_log_L(
          y_1t, this_new_Z, a, b, xi, phi
        )$log_likelihood

        new_L_v <- c(new_L_v, log_new_L)

        this_new_w <- log(w[j]) + log_prior +
          log_new_L -
          last_L_v[j]

        new_w <- c(new_w, this_new_w)
      }
    }

    # softmax
    new_w <- normalize_log_weights(new_w)

    # resampling
    if (length(new_w) > N) {
      which.in <- fc_resample(new_w, N)
      which.in <- which.in[[1]]
      Z <- new_Z[which.in]
      w <- new_w[which.in]
      last_L_v <- new_L_v[which.in]
    } else {
      Z <- new_Z
      w <- new_w
      last_L_v <- new_L_v
    }

    w <- w / sum(w)
    k.curr <- sapply(Z, max)

    cat("data: ", Y[t])
    cat("mean: ", mean(k.curr))
    cat("| weighted mean: ", sum(w * k.curr), "\n")
  }





  return(list(Z, w))
}
