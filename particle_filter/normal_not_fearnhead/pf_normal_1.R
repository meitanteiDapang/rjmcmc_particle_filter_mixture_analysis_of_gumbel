source(file.path(dir_path, "particle_filter/pf_utilities.r"))



plot_resample_mu_sigma <- function(mu, sigma, w,
                                   mu_true = NULL, sigma_true = NULL,
                                   M = NULL,
                                   adjust = 1, 
                                   alpha = 0.25) {           
  # ---- Validate ----
  stopifnot(is.numeric(mu), is.numeric(sigma), is.numeric(w))
  N <- length(w)
  if (length(mu) != N || length(sigma) != N) {
    stop("mu, sigma, and w must have the same length.")
  }
  if (sum(w) <= 0) stop("Sum of weights must be > 0.")
  if (is.null(M)) M <- N
  w <- w / sum(w)

  # ---- Resample with replacement ----
  ix <- sample.int(N, size = M, replace = TRUE, prob = w)
  mu_r <- mu[ix]; sg_r <- sigma[ix]

  # ---- Build data ----
  df_mu <- data.frame(value = mu_r)
  df_sg <- data.frame(value = sg_r)

  library(ggplot2)

  # ---- Helper to plot density ----
  make_density <- function(df, x_expr, truth = NULL, panel_title) {
  g <- ggplot(df, aes(x = value)) +
    geom_density(adjust = adjust, alpha = alpha,
                 fill = "steelblue", colour = "steelblue") +
    theme_minimal(base_size = 12) +
    labs(x = x_expr, y = NULL, title = panel_title) +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  if (!is.null(truth) && is.finite(truth)) {
    g <- g + geom_vline(xintercept = truth, linetype = "dashed",
                        colour = "black", linewidth = 0.6)
  }
  g
}


  p_mu <- make_density(df_mu,
                       expression(mu),
                       mu_true,
                       expression("(a) density of " * underline(mu)))

  p_sg <- make_density(df_sg,
                       expression(sigma),
                       sigma_true,
                       expression("(b) density of " * underline(sigma)))

  combined <- cowplot::plot_grid(p_mu, p_sg, ncol = 1, align = "hv")
  print(combined)
  invisible(combined)
  combined
}




get_initialised_particles_1_normal <- function(
    N,
    mu_0, sigma_0, a, b) {
  # initial mu
  mu <- rnorm(N, mean = mu_0, sd = sigma_0)
  # initial tau
  sigma <- runif(N, a, b)
  # initial w
  w <- rep(1 / N, N)

  return(list(mu = mu, sigma = sigma, w = w))
}

get_index_systematic_resampling <- function(w, N) {
  cum_weights <- cumsum(w)
  # generate u
  u <- runif(1, min = 0, max = 1 / N)
  m <- 1
  new_samples_index <- numeric(N)

  for (j in 1:N) {
    while (u > cum_weights[m]) {
      m <- m + 1
    }
    new_samples_index[j] <- m
    u <- u + 1 / N
  }
  return(new_samples_index)
}


get_index_adaptive_resampling <- function(w, N) {
  #   # calculate ESS
  #   ess <- 1 / sum(w^2)

  #   # Adaptive resampling step based on ESS
  #   if (ess < ess_threshold_ratio * N) {
  #     cum_weights <- cumsum(w)
  #     u <- runif(1, min = 0, max = 1 / N)
  #     m <- 1
  #     new_samples_index <- numeric(N)

  #     for (j in 1:N) {
  #       while (u > cum_weights[m]) {
  #         m <- m + 1
  #       }
  #       new_samples_index[j] <- m
  #       u <- u + 1 / N
  #     }

  #     # resample particles
  #     mu <- mu[new_samples_index]
  #     sigma <- sigma[new_samples_index]
  #     w <- rep(1/N, N) # reset weights after resampling
  #   }
  # }
}


pf_normal_1 <- function(y, N,
                        mu_0, sigma_0, a, b,
                        sd_mu, sd_sigma,
                        verbose = TRUE) {
  # pre-data
  n <- length(y)
  k <- 1

  initialised <- get_initialised_particles_1_normal(N, mu_0, sigma_0, a, b)
  mu <- initialised$mu
  sigma <- initialised$sigma
  w <- initialised$w

  # start loop
  progress_v <- floor(n * c(1:100) * 0.01)
  for (t in c(1:n)) {
    # print info
    if (verbose) {
      if (t %in% progress_v) {
        cat("Iter:", t, "\n")
      }
    }

    # propagating : Nothing


    # weight updating
    likelihood <- dnorm(y[t], mean = mu, sd = sigma) + 1e-10
    w <- w * likelihood
    w <- w / sum(w)

  }
  return(list(mu = mu, sigma = sigma, w = w))
}
