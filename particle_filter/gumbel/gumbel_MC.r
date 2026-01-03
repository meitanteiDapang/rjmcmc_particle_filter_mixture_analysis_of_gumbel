library(patchwork)
library(ggplot2)
library(cowplot)
my_integrand <- function(mu, tau, y_i, eta, phi, a, b) {
  if (tau <= 0) {
    return(0)
  }

  n_i <- length(y_i)
  y_bar <- mean(y_i)

  # Likelihood parts
  log_likelihood <- n_i * log(tau) -
    n_i * tau * (y_bar - mu) -
    sum(exp(-tau * (y_i - mu)))

  # Prior on mu: Normal(eta, phi)
  log_prior_mu <- -0.5 * log(2 * pi * phi) - (mu - eta)^2 / (2 * phi)

  # Prior on tau: Gamma(a, b)
  log_prior_tau <- a * log(b) - lgamma(a) + (a - 1) * log(tau) - b * tau

  # Total log density
  log_total <- log_likelihood + log_prior_mu + log_prior_tau

  # Return unnormalized density
  return(log_total)
}


wrapped_integrand <- function(x, ...) {
  mu <- x[1]
  tau <- x[2]
  exp(my_integrand(mu, tau, ...))
}









laplace_approximation <- function(y_i, eta, phi, a, b, verbose = TRUE) {
  library(numDeriv)

  log_f <- function(theta) {
    mu <- theta[1]
    tau <- theta[2]
    if (tau <= 0 || is.nan(tau) || is.infinite(tau)) return(-Inf)
    if (is.nan(mu) || is.infinite(mu)) return(-Inf)

    n_i <- length(y_i)
    y_bar <- mean(y_i)

    log_likelihood <- n_i * log(tau) -
      n_i * tau * (y_bar - mu) -
      sum(exp(-tau * (y_i - mu)))

    log_prior_mu <- -0.5 * log(2 * pi * phi) - (mu - eta)^2 / (2 * phi)
    log_prior_tau <- a * log(b) - lgamma(a) + (a - 1) * log(tau) - b * tau

    return(log_likelihood + log_prior_mu + log_prior_tau)
  }

  
  v <- var(y_i)
  tau_candidates <- c(1, a / b)
  mu_candidates <- c(mean(y_i))

  init_grid <- expand.grid(mu = mu_candidates, tau = tau_candidates)
  best_result <- NULL
  best_logf <- -Inf

  for (i in 1:nrow(init_grid)) {
    init <- as.numeric(init_grid[i, ])
    if (verbose) cat("Trying init:", init, "\n")

    result <- tryCatch({
      optim(
        par = init,
        fn = function(x) -log_f(x),
        method = "L-BFGS-B",
        lower = c(-Inf, 1e-6),
        control = list(fnscale = 1)
      )
    }, error = function(e) NULL)

    if (!is.null(result) && result$convergence == 0) {
      logf_val <- -result$value
      if (logf_val > best_logf) {
        best_logf <- logf_val
        best_result <- result
      }
    }
  }

  if (is.null(best_result)) {
    warning("All optimization attempts failed.")
    return(NA)
  }

  hat_mu <- best_result$par[1]
  hat_tau <- best_result$par[2]
  log_f_hat <- -best_result$value

  
  H <- hessian(log_f, x = c(hat_mu, hat_tau))
  det_H <- det(H)

  if (det_H <= 0 || is.nan(det_H)) {
    warning("Hessian not positive definite. Result may be invalid.")
    return(NA)
  }

  
  laplace_val <- exp(log_f_hat) * (2 * pi) / sqrt(det_H)
  log_laplace_val <- log_f_hat + log(2 * pi) - 0.5 * log(det_H)

  if (verbose) {
    cat("Mode (mu, tau):", hat_mu, hat_tau, "\n")
    cat("Log f at mode:", log_f_hat, "\n")
    cat("Determinant of Hessian:", det_H, "\n")
    cat("Laplace Approximation:", laplace_val, "\n")
    cat("Log-Laplace Approximation:", log_laplace_val, "\n")
  }

  return(list(
    mu_hat = hat_mu,
    tau_hat = hat_tau,
    log_f_hat = log_f_hat,
    det_H = det_H,
    laplace = laplace_val,
    log_laplace = log_laplace_val
  ))
}


estimate_log_likelihood_given_Z_by_laplace <- function(
    Y, Z,
    eta, phi, a, b,
    plot = FALSE) {
  cluster_ids <- sort(unique(Z))
  total_log_likelihood <- 0
  cluster_results <- list()

  for (cluster_id in cluster_ids) {
    y_i <- Y[Z == cluster_id]
    res <- laplace_approximation(y_i, eta, phi, a, b, verbose = FALSE)
    total_log_likelihood <- total_log_likelihood + res$log_laplace
    cluster_results[[as.character(cluster_id)]] <- res
  }


  return(list(
    log_likelihood = total_log_likelihood,
    cluster_details = cluster_results
  ))
}













estimate_cluster_likelihood_gumbel_MC_uniform <- function(
    y_i, eta, phi, a, b, S = 10000,
    mu_min = -10, mu_max = 10,
    tau_min = 0.01, tau_max = 10) {
  n_i <- length(y_i)
  y_bar_i <- mean(y_i)
  log_weighted_terms <- numeric(S)

  log_q_mu <- -log(mu_max - mu_min)
  log_q_tau <- -log(tau_max - tau_min)

  for (s in 1:S) {
    # Uniform proposal sampling
    mu_s <- runif(1, mu_min, mu_max)
    tau_s <- runif(1, tau_min, tau_max)

    # Likelihood terms
    log_tau_n <- n_i * log(tau_s)
    log_part1 <- -(n_i * tau_s * y_bar_i - n_i * tau_s * mu_s)
    log_part2 <- -sum(exp(-tau_s * (y_i - mu_s)))

    # Prior densities
    log_prior_mu <- dnorm(mu_s, mean = eta, sd = sqrt(phi), log = TRUE)
    log_prior_tau <- dgamma(tau_s, shape = a, rate = b, log = TRUE)

    # Proposal densities (log-uniform)
    log_proposal_mu <- log_q_mu
    log_proposal_tau <- log_q_tau
  if (length(log_tau_n) != 1)   print(paste("log_tau_n length:", length(log_tau_n)))
  if (length(log_part1) != 1)   print(paste("log_part1 length:", length(log_part1)))
  if (length(log_part2) != 1)   print(paste("log_part2 length:", length(log_part2)))
  if (length(log_prior_mu) != 1)  print(paste("log_prior_mu length:", length(log_prior_mu)))
  if (length(log_prior_tau) != 1) print(paste("log_prior_tau length:", length(log_prior_tau)))
  if (length(log_proposal_mu) != 1) print(paste("log_proposal_mu length:", length(log_proposal_mu)))
  if (length(log_proposal_tau) != 1) print(paste("log_proposal_tau length:", length(log_proposal_tau)))
    # Log weight (unnormalised integrand / proposal)
    log_weighted_terms[s] <- log_tau_n + log_part1 + log_part2 +
      log_prior_mu + log_prior_tau -
      log_proposal_mu - log_proposal_tau
  }

  # Only use finite terms to avoid -Inf or NaN
  finite_logs <- log_weighted_terms[is.finite(log_weighted_terms)]


  # hist(finite_logs, breaks = 100, main = "Log weighted terms")


  if (length(finite_logs) == 0) {
    return(list(
      likelihood = 0,
      log_likelihood = -Inf,
      standard_error = NA,
      relative_error = Inf
    ))
  }

  # log-sum-exp trick with stable relative error
  max_log <- max(finite_logs)
  shifted_exp <- exp(finite_logs - max_log)
  mean_shifted <- mean(shifted_exp)
  sd_shifted <- sd(shifted_exp)

  log_L_hat <- max_log + log(mean_shifted)
  L_hat <- exp(log_L_hat)

  se_hat <- sd_shifted / sqrt(length(finite_logs)) * exp(max_log)
  rel_error <- sd_shifted / (sqrt(length(finite_logs)) * mean_shifted)

  return(list(
    likelihood = L_hat,
    log_likelihood = log_L_hat,
    standard_error = se_hat,
    relative_error = rel_error
  ))
}



estimate_likelihood_gumbel_MC_uniform_total <- function(
    Y, Z, eta, phi, a, b, S = 10000,
    mu_min = -10, mu_max = 10,
    tau_min = 0.01, tau_max = 10) {
  cluster_ids <- sort(unique(Z))
  k <- length(cluster_ids)

  log_likelihood_total <- 0
  likelihood_total <- 1
  relative_errors <- numeric(k)

  for (i in seq_along(cluster_ids)) {
    cluster <- cluster_ids[i]
    y_i <- Y[Z == cluster]

    est <- estimate_cluster_likelihood_gumbel_MC_uniform(
      y_i = y_i,
      eta = eta,
      phi = phi,
      a = a,
      b = b,
      S = S,
      mu_min = mu_min,
      mu_max = mu_max,
      tau_min = tau_min,
      tau_max = tau_max
    )

    # Accumulate results
    log_likelihood_total <- log_likelihood_total + est$log_likelihood
    likelihood_total <- likelihood_total * est$likelihood
    relative_errors[i] <- est$relative_error
  }
  cat("relative errors:", relative_errors, "\n")
  return(list(
    likelihood = likelihood_total,
    log_likelihood = log_likelihood_total,
    mean_relative_error = mean(relative_errors)
  ))
}




















estimate_cluster_likelihood_gumbel_MC_normal_gamma <- function(
    y_i, eta, phi, a, b, S = 10000,
    mu_p = 0, sigma_p = 1,
    a_p = 2, b_p = 2,
    plot = FALSE) {
  n_i <- length(y_i)
  y_bar_i <- mean(y_i)
  log_weighted_terms <- numeric(S)

  for (s in 1:S) {
    # Sample from proposal: Normal-Gamma
    mu_s <- rnorm(1, mean = mu_p, sd = sigma_p)
    tau_s <- rgamma(1, shape = a_p, rate = b_p)

    # Likelihood components
    log_tau_n <- n_i * log(tau_s)
    log_part1 <- -n_i * tau_s * (y_bar_i - mu_s)
    log_part2 <- -sum(exp(-tau_s * (y_i - mu_s)))

    # Prior densities
    log_prior_mu <- dnorm(mu_s, mean = eta, sd = sqrt(phi), log = TRUE)
    log_prior_tau <- dgamma(tau_s, shape = a, rate = b, log = TRUE)

    # Proposal densities
    log_q_mu <- dnorm(mu_s, mean = mu_p, sd = sigma_p, log = TRUE)
    log_q_tau <- dgamma(tau_s, shape = a_p, rate = b_p, log = TRUE)

  if (length(log_tau_n) != 1)   print(paste("log_tau_n length:", length(log_tau_n)))
  if (length(log_part1) != 1)   print(paste("log_part1 length:", length(log_part1)))
  if (length(log_part2) != 1)   print(paste("log_part2 length:", length(log_part2)))
  if (length(log_prior_mu) != 1)  print(paste("log_prior_mu length:", length(log_prior_mu)))
  if (length(log_prior_tau) != 1) print(paste("log_prior_tau length:", length(log_prior_tau)))
  if (length(log_q_mu) != 1){
    print(paste("log_q_mu length:", length(log_q_mu)))
    cat("mu_s:", mu_s, " mu_p:", mu_p, " sigma_p:", sigma_p, "\n")
    cat("a_p:", a_p, " b_p:", b_p, " tau_s:", tau_s, "\n")
  } 
  if (length(log_q_tau) != 1) print(paste("log_q_tau length:", length(log_q_tau)))
    # Importance weight (log)
    log_weighted_terms[s] <- log_tau_n + log_part1 + log_part2 +
      log_prior_mu + log_prior_tau - log_q_mu - log_q_tau
  }

  # Filter finite values
  finite_logs <- log_weighted_terms[is.finite(log_weighted_terms)]
  if (plot) hist(finite_logs, breaks = 100, main = "Log weighted terms")

  if (length(finite_logs) == 0) {
    return(list(
      likelihood = 0,
      log_likelihood = -Inf,
      standard_error = NA,
      relative_error = Inf
    ))
  }

  # log-sum-exp
  max_log <- max(finite_logs)
  shifted_exp <- exp(finite_logs - max_log)
  mean_shifted <- mean(shifted_exp)
  sd_shifted <- sd(shifted_exp)

  log_L_hat <- max_log + log(mean_shifted)
  L_hat <- exp(log_L_hat)

  se_hat <- sd_shifted / sqrt(length(finite_logs)) * exp(max_log)
  rel_error <- sd_shifted / (sqrt(length(finite_logs)) * mean_shifted)


  # cat("Range of log_weighted_terms: ", range(finite_logs), "\n")
  # cat("Mean of log_weighted_terms: ", mean(finite_logs), "\n")
  # cat("SD of log_weighted_terms: ", sd(finite_logs), "\n")
  # cat("Summary of log_weighted_terms:\n")
  # print(summary(finite_logs))


  return(list(
    likelihood = L_hat,
    log_likelihood = log_L_hat,
    standard_error = se_hat,
    relative_error = rel_error
  ))
}



estimate_gamma_proposal <- function(Y, my_rate, a, b, base_fraction = 0.5) {
  NN <- length(Y)
  if (NN < 2) {
    tau_hat <- a / b
  } else {
    sigma_hat_2 <- 1 / NN * sum((Y - mean(Y))^2)
    beta_hat <- sqrt(sigma_hat_2 * 6 / real_pi^2)
    tau_hat <- 1 / beta_hat
  }

  tau_sd <- base_fraction * tau_hat / NN^(my_rate)
  a_p <- (tau_hat / tau_sd)^2
  b_p <- a_p / tau_hat
  return(list(a_p = a_p, b_p = b_p, tau_mean = tau_hat, tau_sd = tau_sd))
}



estimate_log_likelihood_given_Z_by_MC <- function(
    Y, Z,
    eta, phi, a, b,
    S = 10000, S_U = 50000,
    plot = FALSE) {
  cluster_ids <- sort(unique(Z))
  total_log_likelihood <- 0
  cluster_results <- list()

  for (cluster_id in cluster_ids) {
    y_i <- Y[Z == cluster_id]
    # cat("y_i:", y_i, "\n")
    NN <- length(y_i)
    my_rate <- 0.25

    # Step 2: proposal parameters

    # If NN < 4, use the uniform proposal
    if (NN < 5) {
      mu_hat <- mean(y_i)
      mu_range <- 2 * abs(mu_hat - eta)
      mu_min <- mu_hat - mu_range
      mu_max <- mu_hat + mu_range


      tau_min <- 0.01
      tau_max <- 2 * a/b
      # cat(mu_min, mu_max, tau_min, tau_max, "\n")
      res <- estimate_cluster_likelihood_gumbel_MC_uniform(
        y_i = y_i,
        eta = eta, phi = phi,
        a = a, b = b,
        mu_min = mu_min, mu_max = mu_max,
        tau_min = tau_min, tau_max = tau_max,
        S = S_U
      )
    } else {
      # If NN >= 4, use the normal-gamma proposal
      sigma_hat_2 <- 1 / NN * sum((y_i - mean(y_i))^2)
      beta_hat <- sqrt(sigma_hat_2 * 6 / real_pi^2)
      mu_hat <- mean(y_i) - beta_hat * 0.5772
      mu_p <- mu_hat
      sigma_p <- 1 / NN^(my_rate)
      gamma_params <- estimate_gamma_proposal(y_i, my_rate, a, b)
      a_p <- gamma_params$a_p
      b_p <- gamma_params$b_p
      # cat("a_p:", a_p, " b_p:", b_p, "\n")
      # cat("sigma_hat_2:", sigma_hat_2, " beta_hat:", beta_hat, "\n")
      res <- estimate_cluster_likelihood_gumbel_MC_normal_gamma(
        y_i = y_i,
        eta = eta, phi = phi,
        a = a, b = b,
        S = S,
        mu_p = mu_p, sigma_p = sigma_p,
        a_p = a_p, b_p = b_p,
        plot = plot
      )
    }




    total_log_likelihood <- total_log_likelihood + res$log_likelihood
    cluster_results[[as.character(cluster_id)]] <- res
    # print(res)
  }

  return(list(
    log_likelihood = total_log_likelihood,
    cluster_details = cluster_results
  ))
}


run_mc_experiment_prior <- function(NN_values, tau_values, S_values, n_rep, eta, phi, a, b) {
  all_results <- list()
  plot_data_list <- list()
  for (NN in NN_values) {
    plot_data <- data.frame()
    for (tau_true in tau_values) {
      Y <- rgumbel(NN, loc = 1, scale = 1 / tau_true)
      results <- data.frame(
        NN = NN,
        tau_true = tau_true,
        S = S_values,
        log_likelihood = NA,
        relative_error = NA
      )
      for (i in seq_along(S_values)) {
        S <- S_values[i]
        logs <- numeric(n_rep)
        rels <- numeric(n_rep)
        for (j in 1:n_rep) {
          res <- estimate_cluster_likelihood_gumbel_MC_prior(
            y_i = Y,
            eta = eta, phi = phi,
            a = a, b = b,
            S = S
          )
          logs[j] <- res$log_likelihood
          rels[j] <- res$relative_error
        }
        results$log_likelihood[i] <- mean(logs, na.rm = TRUE)
        results$relative_error[i] <- mean(rels, na.rm = TRUE)
      }
      all_results[[paste0("NN_", NN, "_tau_", tau_true)]] <- results
      plot_data <- rbind(plot_data, results)
    }
    plot_data_list[[as.character(NN)]] <- plot_data
  }
  return(list(all_results = all_results, plot_data_list = plot_data_list))
}

run_mc_experiment_normal_gamma <- function(NN_values, tau_values, S_values, n_rep, eta, phi, a, b, my_rate) {
  all_results <- list()
  plot_data_list <- list()
  for (NN in NN_values) {
    plot_data <- data.frame()
    for (tau_true in tau_values) {
      Y <- rgumbel(NN, loc = 1, scale = 1 / tau_true)
      beta_hat <- sqrt(var(Y) * 6 / real_pi^2)
      mu_hat <- mean(Y) - beta_hat * 0.5772
      mu_p <- mu_hat
      sigma_p <- 1 / NN^(my_rate)
      gamma_params <- estimate_gamma_proposal(Y, NN)
      results <- data.frame(
        NN = NN,
        tau_true = tau_true,
        S = S_values,
        log_likelihood = NA_real_,
        relative_error = NA_real_
      )
      for (i in seq_along(S_values)) {
        S <- S_values[i]
        logs <- numeric(n_rep)
        rels <- numeric(n_rep)
        for (j in 1:n_rep) {
          res <- estimate_cluster_likelihood_gumbel_MC_normal_gamma(
            y_i = Y,
            eta = eta, phi = phi,
            a = a, b = b,
            mu_p = mu_p, sigma_p = sigma_p,
            a_p = gamma_params$a_p, b_p = gamma_params$b_p,
            S = S
          )
          logs[j] <- res$log_likelihood
          rels[j] <- res$relative_error
        }
        results$log_likelihood[i] <- mean(logs, na.rm = TRUE)
        results$relative_error[i] <- mean(rels, na.rm = TRUE)
      }
      all_results[[paste0("NN_", NN, "_tau_", tau_true)]] <- results
      plot_data <- rbind(plot_data, results)
    }
    plot_data_list[[as.character(NN)]] <- plot_data
  }
  list(all_results = all_results, plot_data_list = plot_data_list)
}

run_mc_experiment_uniform <- function(NN_values, tau_values, S_values, n_rep, eta, phi, a, b) {
  all_results <- list()
  plot_data_list <- list()
  
  for (NN in NN_values) {
    plot_data <- data.frame()
    
    for (tau_true in tau_values) {
      Y <- rgumbel(NN, loc = 1, scale = 1 / tau_true)
      
      mu_hat <- mean(Y)
      mu_range <- 2 * abs(mu_hat - eta)
      mu_min <- mu_hat - mu_range
      mu_max <- mu_hat + mu_range
      
      tau_min <- 0.01
      tau_max <- 2 * phi
      
      results <- data.frame(
        NN = NN,
        tau_true = tau_true,
        S = S_values,
        log_likelihood = NA_real_,
        relative_error = NA_real_
      )
      
      for (i in seq_along(S_values)) {
        S <- S_values[i]
        logs <- numeric(n_rep)
        rels <- numeric(n_rep)
        
        for (j in 1:n_rep) {
          res <- estimate_cluster_likelihood_gumbel_MC_uniform(
            y_i = Y,
            eta = eta, phi = phi,
            a = a, b = b,
            mu_min = mu_min, mu_max = mu_max,
            tau_min = tau_min, tau_max = tau_max,
            S = S
          )
          logs[j] <- res$log_likelihood
          rels[j] <- res$relative_error
        }
        
        results$log_likelihood[i] <- mean(logs, na.rm = TRUE)
        results$relative_error[i] <- mean(rels, na.rm = TRUE)
      }
      
      all_results[[paste0("NN_", NN, "_tau_", tau_true)]] <- results
      plot_data <- rbind(plot_data, results)
    }
    
    plot_data_list[[as.character(NN)]] <- plot_data
  }
  
  list(all_results = all_results, plot_data_list = plot_data_list)
}



make_mc_plots <- function(plot_data_list, n_values) {
  plots <- lapply(seq_along(n_values), function(i) {
    NN <- n_values[i]
    p <- ggplot(plot_data_list[[as.character(NN)]],
                aes(x = S, y = relative_error, colour = factor(tau_true))) +
      geom_line() +
      geom_point() +
      scale_x_log10() +
      ylab("RE") +
      xlab("S") +
      labs(colour = expression(tau^{true})) +
      theme_minimal() +
      ggtitle(paste0("(", letters[i], ") n=", NN))
    return(p)
  })
  
  legend <- get_legend(
    plots[[1]] + theme(legend.position = "bottom")
  )
  
  plots <- lapply(plots, function(p) p + theme(legend.position = "none"))
  
  grid <- plot_grid(plotlist = plots, ncol = 2, labels = NULL)
  final_plot <- plot_grid(grid, legend, ncol = 1, rel_heights = c(1, 0.1))
  
  print(final_plot)
  return(final_plot)
}


run_multiple_mc_estimates <- function(
    y_i, eta, phi, a, b,
    mu_p, sigma_p,
    a_p, b_p,
    S_values = c(100, 500, 1000, 5000, 10000, 20000, 50000),
    n_rep = 10) {
  results_list <- list()

  for (S in S_values) {
    cat("\n========== S =", S, " ==========\n")
    estimates <- numeric(n_rep)

    for (i in 1:n_rep) {
      set.seed(1000 + i) 
      res <- estimate_cluster_likelihood_gumbel_MC_normal_gamma(
        y_i = y_i,
        eta = eta, phi = phi,
        a = a, b = b,
        mu_p = mu_p, sigma_p = sigma_p,
        a_p = a_p, b_p = b_p,
        S = S,
        plot = FALSE
      )
      estimates[i] <- res$likelihood
    }

    mean_est <- mean(estimates)
    sd_est <- sd(estimates)
    rel_err <- sd_est / mean_est

    results_list[[as.character(S)]] <- data.frame(
      S = S,
      mean = mean_est,
      sd = sd_est,
      relative_error = rel_err
    )
  }

  results_df <- bind_rows(results_list)
  return(results_df)
}




























estimate_cluster_likelihood_gumbel_MC_prior <- function(
    y_i, eta, phi, a, b, S = 10000) {
  n_i <- length(y_i)
  y_bar_i <- mean(y_i)
  log_terms <- numeric(S)

  for (s in 1:S) {
    # Sample from prior
    mu_s <- rnorm(1, mean = eta, sd = sqrt(phi))
    tau_s <- rgamma(1, shape = a, rate = b)

    # Compute log-likelihood-like terms
    log_term1 <- n_i * log(tau_s)
    log_term2 <- -n_i * tau_s * (y_bar_i - mu_s)
    log_term3 <- -sum(exp(-tau_s * (y_i - mu_s)))

    log_terms[s] <- log_term1 + log_term2 + log_term3
  }

  # Filter finite values
  finite_logs <- log_terms[is.finite(log_terms)]

  if (length(finite_logs) == 0) {
    return(list(
      likelihood = 0,
      log_likelihood = -Inf,
      standard_error = NA,
      relative_error = Inf
    ))
  }

  # log-sum-exp
  max_log <- max(finite_logs)
  shifted_exp <- exp(finite_logs - max_log)
  mean_shifted <- mean(shifted_exp)
  sd_shifted <- sd(shifted_exp)

  log_L_hat <- max_log + log(mean_shifted)
  L_hat <- exp(log_L_hat)

  se_hat <- sd_shifted / sqrt(length(finite_logs)) * exp(max_log)
  rel_error <- sd_shifted / (sqrt(length(finite_logs)) * mean_shifted)


  # cat("Range of log_weighted_terms: ", range(finite_logs), "\n")
  # cat("Mean of log_weighted_terms: ", mean(finite_logs), "\n")
  # cat("SD of log_weighted_terms: ", sd(finite_logs), "\n")
  # cat("Summary of log_weighted_terms:\n")
  # print(summary(finite_logs))


  return(list(
    likelihood = L_hat,
    log_likelihood = log_L_hat,
    standard_error = se_hat,
    relative_error = rel_error
  ))
}
