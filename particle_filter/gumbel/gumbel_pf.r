
source(file.path(dir_path, "particle_filter/gumbel/gumbel_MC.r"))
source(file.path(dir_path, "particle_filter/fearnhead_utilities.r"))
source(file.path(dir_path, "particle_filter/fc_resample.r"))
source(file.path(dir_path, "particle_filter/pf_utilities.r"))

pf_gumbel <- function(Y, N,
                      a, b, eta, phi, alpha,
                      initial_size = N,
                      S = 10000, S_U = 50000,
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
        cat("\nIter", t, ": \n")
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
      flush.console()
      # the i-th particle after t-1 th step
      this_Zj <- Z[[j]]

      # get the possible next allocation at t-th step
      possible_next_allocation <- 1 + k.curr[j]

      tcat("\nthis_Zj: ", this_Zj, " log: ", last_L_v[j], " last_log_w:", log(w[j]), "-----------------------\n")
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
        log_new_L <- estimate_log_likelihood_given_Z_by_MC(
          y_1t, this_new_Z,
          eta, phi, a, b,
          S = S,S_U = S_U,
          plot = FALSE
        )$log_likelihood



        new_L_v <- c(new_L_v, log_new_L)

        this_new_w <- log(w[j]) + log_prior +
          log_new_L -
          last_L_v[j]

        tcat("The assign: ", i, " log_prior: ", log_prior, " log_new_L: ", log_new_L, "new log w: ", this_new_w, "\n")
        new_w <- c(new_w, this_new_w)
      }

      cat(sprintf("\rProcessing: %d / %d, ", j, length(Z)))
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

    tcat("---final w: ", w, "\n")

    cluster_counts <- sapply(Z, function(z) length(unique(z)))
    cat("data: ", Y[t])
    cat(", mean: ", mean(k.curr))
    cat("| weighted mean: ", sum(w * k.curr), "\n")
  }





  return(list(Z, w))
}


pf_gumbel_laplace <- function(Y, N,
                      a, b, eta, phi, alpha,
                      initial_size = N,
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
        cat("\nIter", t, ": \n")
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
      flush.console()
      # the i-th particle after t-1 th step
      this_Zj <- Z[[j]]

      # get the possible next allocation at t-th step
      possible_next_allocation <- 1 + k.curr[j]

      tcat("\nthis_Zj: ", this_Zj, " log: ", last_L_v[j], " last_log_w:", log(w[j]), "-----------------------\n")
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
        log_new_L <- estimate_log_likelihood_given_Z_by_laplace(
          y_1t, this_new_Z,
          eta, phi, a, b,
          plot = FALSE
        )$log_likelihood



        new_L_v <- c(new_L_v, log_new_L)

        this_new_w <- log(w[j]) + log_prior +
          log_new_L -
          last_L_v[j]

        tcat("The assign: ", i, " log_prior: ", log_prior, " log_new_L: ", log_new_L, "new log w: ", this_new_w, "\n")
        new_w <- c(new_w, this_new_w)
      }

      cat(sprintf("\rProcessing: %d / %d, ", j, length(Z)))
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

    tcat("---final w: ", w, "\n")

    cluster_counts <- sapply(Z, function(z) length(unique(z)))
    cat("data: ", Y[t])
    cat(", mean: ", mean(k.curr))
    cat("| weighted mean: ", sum(w * k.curr), "\n")
  }





  return(list(Z, w))
}
