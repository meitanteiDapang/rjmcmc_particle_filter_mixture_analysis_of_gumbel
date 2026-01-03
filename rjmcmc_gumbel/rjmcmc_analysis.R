library(evd)
library(truncnorm)


get_bias_k <- function(est_list, N, n, index,
                       true_value) {
  sum <- 0
  for (i in c(1:N)) {
    for (j in c(1:n)) {
      sum <- sum + (est_list[[i]][j, index] - true_value) / (N * n)
    }
  }
  return(sum)
}

get_mse_k <- function(est_list, N, n, index,
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

get_coverage_table_k <- function(est_list, N, n, index, CL,
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

analyse_k_gumbel <- function(sim_mu, sim_tau, my_list1, N, iter, burn.in, k) {
  # Order sim_mu and reorder sim_mu and sim_tau accordingly
  order_index <- order(sim_mu)
  ordered_mu <- sim_mu[order_index]
  ordered_tau <- sim_tau[order_index]

  # Calculate bias for mu
  mu_bias <- sapply(1:k, function(j) {
    get_bias_k(my_list1, N, (iter - burn.in), j, ordered_mu[j])
  })
  cat("mu bias: ", mu_bias, "\n")

  # Calculate bias for tau
  tau_bias <- sapply(1:k, function(j) {
    get_bias_k(my_list1, N, (iter - burn.in), j + k, ordered_tau[j])
  })
  cat("tau bias: ", tau_bias, "\n")

  # Calculate Mean Squared Error (MSE) for mu
  mu_mse <- sapply(1:k, function(j) {
    get_mse_k(my_list1, N, (iter - burn.in), j, ordered_mu[j])
  })
  cat("mu mse: ", mu_mse, "\n")

  # Calculate Mean Squared Error (MSE) for tau
  tau_mse <- sapply(1:k, function(j) {
    get_mse_k(my_list1, N, (iter - burn.in), j + k, ordered_tau[j])
  })
  cat("tau mse: ", tau_mse, "\n")

  # Source the file that contains the definition of get_coverage_table
  source("single_gumbel_analysis.R")

  # Calculate coverage tables and coverage probabilities for mu
  mu_coverage_table_list <- lapply(1:k, function(j) {
    get_coverage_table(my_list1, N, (iter - burn.in), j, 0.95, ordered_mu[j])
  })
  mu_coverage <- sapply(1:k, function(j) {
    sum(ordered_mu[j] >= mu_coverage_table_list[[j]]$lower &
      ordered_mu[j] <= mu_coverage_table_list[[j]]$upper) / N
  })

  # Calculate coverage tables and coverage probabilities for tau
  tau_coverage_table_list <- lapply(1:k, function(j) {
    get_coverage_table(my_list1, N, (iter - burn.in), j + k, 0.95, ordered_tau[j])
  })
  tau_coverage <- sapply(1:k, function(j) {
    sum(ordered_tau[j] >= tau_coverage_table_list[[j]]$lower &
      ordered_tau[j] <= tau_coverage_table_list[[j]]$upper) / N
  })

  cat("coverage of mu: ", mu_coverage, "\ncoverage of tau: ", tau_coverage, "\n")

  # Return all computed results including coverage tables
  return(list(
    mu_bias = mu_bias,
    tau_bias = tau_bias,
    mu_mse = mu_mse,
    tau_mse = tau_mse,
    mu_coverage = mu_coverage,
    tau_coverage = tau_coverage,
    mu_coverage_table_list = mu_coverage_table_list,
    tau_coverage_table_list = tau_coverage_table_list
  ))
}
