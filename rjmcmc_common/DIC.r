compute_DIC <- function(k, est, X, burn.in, iter) {
  est <- est[(burn.in + 1):iter, , drop = FALSE]
  N <- nrow(est) # N: number of MCMC samples after burn-in
  n <- length(X) # n: sample size

  # --- D.bar ---
  sum_D <- 0
  for (i in 1:N) { # i = MCMC draw index
    sum_Di <- 0
    mu <- as.numeric(est[i, 1:k, drop = FALSE])
    tau <- as.numeric(est[i, (k + 1):(2 * k), drop = FALSE])
    pi <- as.numeric(est[i, (2 * k + 1):(3 * k), drop = FALSE])
    for (t in 1:n) { # t = data index
      dens <- pi * dnorm(X[t], mean = mu, sd = 1 / sqrt(tau))
      sum_Di <- sum_Di - 2 * log(sum(dens))
    }
    sum_D <- sum_D + sum_Di
  }
  D.bar <- sum_D / N

  # --- D(theta_hat) ---
  mu.bar <- as.numeric(colMeans(est[, 1:k, drop = FALSE]))
  tau.bar <- as.numeric(colMeans(est[, (k + 1):(2 * k), drop = FALSE]))
  pi.bar <- as.numeric(colMeans(est[, (2 * k + 1):(3 * k), drop = FALSE]))
  D.theta.hat <- 0
  for (t in 1:n) { # t = data index
    dens <- pi.bar * dnorm(X[t], mean = mu.bar, sd = 1 / sqrt(tau.bar))
    D.theta.hat <- D.theta.hat - 2 * log(sum(dens))
  }

  PD <- D.bar - D.theta.hat
  cat("PD =", PD, "\n")
  DIC <- D.bar + PD
  return(DIC)
}



compute_k_gumbel_DIC <- function(k, est, Y, burn.in, iter) {
  # Take post–burn-in draws; robust to thinning
  est <- est[(burn.in + 1):iter, , drop = FALSE]
  N   <- nrow(est)
  n   <- length(Y)

  # Structural check: columns must be [mu(1:k), tau(1:k), w(1:k)]
  if (ncol(est) < 3 * k) {
    stop("est must have at least 3*k columns: mu(1:k), tau(1:k), w(1:k).")
  }

  # --- D.bar ---
  sum_D <- 0
  for (i in 1:N) {
    mu  <- as.numeric(est[i, 1:k, drop = FALSE])
    tau <- as.numeric(est[i, (k + 1):(2 * k), drop = FALSE])
    w   <- as.numeric(est[i, (2 * k + 1):(3 * k), drop = FALSE])

    Di <- 0
    sc <- 1 / tau
    for (t in 1:n) {
      dens_sum <- sum(w * dgumbel(Y[t], loc = mu, scale = sc))
      Di <- Di - 2 * log(dens_sum)
    }
    sum_D <- sum_D + Di
  }
  D.bar <- sum_D / N

  # --- D(theta_hat) ---
  mu.hat  <- as.numeric(colMeans(est[, 1:k, drop = FALSE]))
  tau.hat <- as.numeric(colMeans(est[, (k + 1):(2 * k), drop = FALSE]))
  w.hat   <- as.numeric(colMeans(est[, (2 * k + 1):(3 * k), drop = FALSE]))

  sc.hat <- 1 / tau.hat

  D.theta.hat <- 0
  for (t in 1:n) {
    dens_sum <- sum(w.hat * dgumbel(Y[t], loc = mu.hat, scale = sc.hat))
    D.theta.hat <- D.theta.hat - 2 * log(dens_sum)
  }

  pd  <- D.bar - D.theta.hat
  cat("PD =", pd, "\n")
  DIC <- D.bar + pd
  return(DIC)
}
