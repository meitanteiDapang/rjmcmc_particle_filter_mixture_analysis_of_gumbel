source(file.path(dir_path, "rjmcmc_common/utilities.r"))

is_debug <- FALSE
debug_start <- 396332
low_threshold <- 1e-5
tcat <- function(...) {
  if (is_debug) {
    cat(...)
  }
}

ess_threshold_ratio <- 0.5

asadf <- TRUE
asdf <- FALSE


get_estimate_mean <- function(mu, sigma, w, k = 1) {
  return(c(sum(mu * w), sum(sigma * w) ))
}

get_estimate_mean_k <- function(mu, sigma, pi_, w) {
  N <- nrow(mu)
  
  w <- w / sum(w)
  
  # Weighted column means
  mu_est     <- colSums(mu * w)
  sigma_est  <- colSums(sigma * w)
  pi_est     <- colSums(pi_ * w)

  return(list(mu = mu_est, sigma = sigma_est, pi_ = pi_est))
}





log_normalize <- function(log_x) {
  max_log_x <- max(log_x) 
  log_sum_exp_x <- max_log_x + log(sum(exp(log_x - max_log_x)))

  normalized_x <- exp(log_x - log_sum_exp_x)

  return(normalized_x)
}


solve_alpha <- function(n, expect) {
  # Define the function to find the root
  f <- function(alpha, n) {
    alpha * log((n + alpha) / alpha) - expect
  }

  # Use uniroot to find alpha
  solution <- uniroot(f, lower = 0.01, upper = 100, n = n)

  return(solution$root)
}


estimate_k <- function(particles, weights, N) {
  res <- 0
  for (i in c(1:N)) {
    res <- res + weights[i] * length(unique(particles[[i]]))
  }
  return(res)
}

get_all_by_k <- function(particles, weights, N, k) {
  index <- c()
  for (i in c(1:N)) {
    if (length(unique(particles[[i]])) == k) {
      index <- c(index, i)
    }
  }

  return(list(particles[index], weights[index]))
}

estimate_mean_sigma_w <- function(particles, weights, Y, k) {
  weights <- weights/sum(weights)
  res <- 0
  final_mean_v <- rep(0, k)
  final_sigma_v <- rep(0, k)
  final_weight_v <- rep(0, k)

  N <- length(particles)
  for (i in c(1:N)) {
    this_particles <- particles[[i]]
    mean_v <- rep(0, k)
    sigma_v <- rep(0, k)
    weight_v <- rep(0, k)

    component_w <- weights[i]
    for (j in c(1:k)){
      index <- which(this_particles == j)
      Y_g <- Y[index]
      mean_v[j] <- mean(Y_g)
      if (length(Y_g) > 1) {
        sigma_v[j] <- sd(Y_g)
      } else {
        sigma_v[j] <- 0.0001
      }
      weight_v[j] <- length(index)
    }
    weight_v <- weight_v/length(Y)

    # add up by weight
    final_mean_v <- final_mean_v + mean_v * component_w
    final_sigma_v <- final_sigma_v + sigma_v * component_w
    final_weight_v <- final_weight_v + weight_v * component_w
  }

  return (list(final_mean_v, final_sigma_v, final_weight_v))
}
