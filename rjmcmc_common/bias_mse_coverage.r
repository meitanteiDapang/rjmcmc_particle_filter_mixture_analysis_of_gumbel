source(file.path(dir_path, "rjmcmc_common/utilities.r"))


get_bias_single <- function(est_list, N, n, index,
                            true_value) {
  sum <- 0
  for (i in c(1:N)) {
    for (j in c(1:n)) {
      sum <- sum + (est_list[[i]][j, index] - true_value) / (N * n)
    }
  }
  return(sum)
}

get_mse_single <- function(est_list, N, n, index,
                           true_value) {
  sum <- 0
  for (i in c(1:N)) {
    est_para <- 0
    for (j in c(1:n)) {
      est_para <- est_para + est_list[[i]][j, index] / n
    }
    sum <- sum + (est_para - true_value)^2 / N
  }
  return(sum)
}

get_coverage_table <- function(est_list, N, n, index, CL,
                               true_value) {
  df <- data.frame(
    x = 1:N,
    lower = rep(0, N),
    upper = rep(0, N),
    m <- rep(0, N)
  )

  for (i in c(1:N)) {
    this_data <- est_list[[i]][, index]
    qs <- quantile(this_data, probs = c((1 - CL) / 2, 1 - (1 - CL) / 2))
    df$lower[df$x == i] <- qs[1]
    df$upper[df$x == i] <- qs[2]
    df$m[df$x == i] <- mean(this_data)
  }

  return(df)
}


analyse_gumbel <- function(sim_mu, sim_tau, sim_pi,
                           est_list, N, n, k, burn.in) {
  order_index <- order(sim_mu)
  ordered_mu <- sim_mu[order_index]
  ordered_tau <- sim_tau[order_index]
  ordered_pi <- sim_pi[order_index]

  cat("burn.in:", burn.in, "\n")
  cat("n:", n, "\n")

  ## ----- Bias -----
  mu_bias <- sapply(1:k, function(j) {
    get_bias_single(est_list, N, (iter - burn.in), j, ordered_mu[j])
  })
  tau_bias <- sapply(1:k, function(j) {
    get_bias_single(est_list, N, (iter - burn.in), j + k, ordered_tau[j])
  })
  pi_bias <- sapply(1:k, function(j) {
    get_bias_single(est_list, N, (iter - burn.in), j + 2 * k, ordered_pi[j])
  })

  cat("mu bias: ", mu_bias, "\n")
  cat("tau bias: ", tau_bias, "\n")
  cat("pi bias: ", pi_bias, "\n")

  ## ----- MSE -----
  mu_mse <- sapply(1:k, function(j) {
    get_mse_single(est_list, N, (iter - burn.in), j, ordered_mu[j])
  })
  tau_mse <- sapply(1:k, function(j) {
    get_mse_single(est_list, N, (iter - burn.in), j + k, ordered_tau[j])
  })
  pi_mse <- sapply(1:k, function(j) {
    get_mse_single(est_list, N, (iter - burn.in), j + 2 * k, ordered_pi[j])
  })

  cat("mu mse: ", mu_mse, "\n")
  cat("tau mse: ", tau_mse, "\n")
  cat("pi mse: ", pi_mse, "\n")

  ## ----- Coverage -----
  mu_coverage_table_list <- lapply(1:k, function(j) {
    get_coverage_table(est_list, N, (iter - burn.in), j, 0.95, ordered_mu[j])
  })
  mu_coverage <- sapply(1:k, function(j) {
    sum(ordered_mu[j] >= mu_coverage_table_list[[j]]$lower &
          ordered_mu[j] <= mu_coverage_table_list[[j]]$upper) / N
  })

  tau_coverage_table_list <- lapply(1:k, function(j) {
    get_coverage_table(est_list, N, (iter - burn.in), j + k, 0.95, ordered_tau[j])
  })
  tau_coverage <- sapply(1:k, function(j) {
    sum(ordered_tau[j] >= tau_coverage_table_list[[j]]$lower &
          ordered_tau[j] <= tau_coverage_table_list[[j]]$upper) / N
  })

  pi_coverage_table_list <- lapply(1:k, function(j) {
    get_coverage_table(est_list, N, (iter - burn.in), j + 2 * k, 0.95, ordered_pi[j])
  })
  pi_coverage <- sapply(1:k, function(j) {
    sum(ordered_pi[j] >= pi_coverage_table_list[[j]]$lower &
          ordered_pi[j] <= pi_coverage_table_list[[j]]$upper) / N
  })

  cat("coverage of mu: ", mu_coverage, "\n")
  cat("coverage of tau: ", tau_coverage, "\n")
  cat("coverage of pi: ", pi_coverage, "\n")

  return(list(
    mu_bias = mu_bias,
    tau_bias = tau_bias,
    pi_bias = pi_bias,
    mu_mse = mu_mse,
    tau_mse = tau_mse,
    pi_mse = pi_mse,
    mu_coverage = mu_coverage,
    tau_coverage = tau_coverage,
    pi_coverage = pi_coverage,
    mu_coverage_table_list = mu_coverage_table_list,
    tau_coverage_table_list = tau_coverage_table_list,
    pi_coverage_table_list = pi_coverage_table_list
  ))
}
