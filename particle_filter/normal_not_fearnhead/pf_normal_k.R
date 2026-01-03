source(file.path(dir_path,"particle_filter/normal_not_fearnhead/pf_normal_1.r"))
source(file.path(dir_path, "particle_filter/pf_utilities.r"))






lot_resample_mix_params <- function(mu_mat, sigma_mat, pi_mat, w,
                                    mu_true = NULL, sigma_true = NULL, pi_true = NULL,
                                    M = NULL,
                                    adjust = 1,
                                    alpha = 0.35,
                                    fill_col = "steelblue",
                                    panel_label_size = 12,
                                    grid_cols = 3) {
  # ---- Validate ----
  stopifnot(is.numeric(w), length(dim(mu_mat)) == 2, length(dim(sigma_mat)) == 2, length(dim(pi_mat)) == 2)
  mu_mat    <- as.matrix(mu_mat)
  sigma_mat <- as.matrix(sigma_mat)
  pi_mat    <- as.matrix(pi_mat)

  N <- length(w)
  stopifnot(nrow(mu_mat) == N, nrow(sigma_mat) == N, nrow(pi_mat) == N)
  k <- ncol(mu_mat)
  stopifnot(ncol(sigma_mat) == k, ncol(pi_mat) == k)
  if (sum(w) <= 0) stop("Sum of weights must be > 0.")
  if (is.null(M)) M <- N

  # ---- Normalise weights & resample indices ----
  w  <- w / sum(w)
  ix <- sample.int(N, size = M, replace = TRUE, prob = w)

  # ---- Slice resampled parameters ----
  mu_r    <- mu_mat[ix, , drop = FALSE]
  sigma_r <- sigma_mat[ix, , drop = FALSE]
  pi_r    <- pi_mat[ix, , drop = FALSE]

  # ---- Helpers ----
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Need 'ggplot2'.")
  if (!requireNamespace("cowplot", quietly = TRUE)) stop("Need 'cowplot'.")

  # 标题表达式：(a) density of underline(mu[j]) / sigma / pi
  make_title <- function(sym, j, lab) {
    # sym 为 "mu"/"sigma"/"pi"
    sym_expr <- as.name(sym)
    bquote(.(lab) * " density of " * underline(.(sym_expr)[.(j)]))
  }

  # x 轴标签表达式：mu[j] / sigma[j] / pi[j]
  make_xlab <- function(sym, j) {
    sym_expr <- as.name(sym)
    bquote(.(sym_expr)[.(j)])
  }

  make_density <- function(x, x_lab, truth = NULL, panel_title, xlim01 = FALSE) {
    g <- ggplot2::ggplot(data.frame(value = x), ggplot2::aes(x = value)) +
      ggplot2::geom_density(adjust = adjust, alpha = alpha,
                            fill = fill_col, colour = fill_col, linewidth = 0.6) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(x = x_lab, y = NULL, title = panel_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))
    if (!is.null(truth) && is.finite(truth)) {
      g <- g + ggplot2::geom_vline(xintercept = truth, linetype = "dashed",
                                   colour = "black", linewidth = 0.6)
    }
    if (xlim01) g <- g + ggplot2::coord_cartesian(xlim = c(0, 1))
    g
  }

  # ---- Build panel list: mu1..k, sigma1..k, pi1..k ----
  panels  <- list()
  lab_idx <- 0L
  next_lab <- function() { lab_idx <<- lab_idx + 1L; paste0("(", letters[lab_idx], ")") }

  # mu panels
  for (j in seq_len(k)) {
    title_j <- make_title("mu", j, next_lab())
    xlab_j  <- make_xlab ("mu", j)
    truth_j <- if (!is.null(mu_true) && length(mu_true) >= j) mu_true[j] else NULL
    panels[[length(panels) + 1L]] <- make_density(mu_r[, j], xlab_j, truth = truth_j, panel_title = title_j)
  }

  # sigma panels
  for (j in seq_len(k)) {
    title_j <- make_title("sigma", j, next_lab())
    xlab_j  <- make_xlab ("sigma", j)
    truth_j <- if (!is.null(sigma_true) && length(sigma_true) >= j) sigma_true[j] else NULL
    panels[[length(panels) + 1L]] <- make_density(sigma_r[, j], xlab_j, truth = truth_j, panel_title = title_j)
  }

  # pi panels (simplex)
  for (j in seq_len(k)) {
    title_j <- make_title("pi", j, next_lab())
    xlab_j  <- make_xlab ("pi", j)
    truth_j <- if (!is.null(pi_true) && length(pi_true) >= j) pi_true[j] else NULL
    panels[[length(panels) + 1L]] <- make_density(pi_r[, j], xlab_j, truth = truth_j, panel_title = title_j, xlim01 = TRUE)
  }

  # ---- Arrange grid ----
  combined <- cowplot::plot_grid(plotlist = panels, ncol = grid_cols, align = "hv")
  print(combined)
  invisible(combined)
  combined
}






get_initialised_particles_k_normal <- function(
    N, k,
    mu_0, sigma_0, a, b, delta) {

  mu <- matrix(NA, nrow = N, ncol = k)
  mu <- t(apply(matrix(rnorm(N * k, mu_0, sigma_0), nrow = N), 1, sort))

  sigma <- matrix(NA, nrow = N, ncol = k)
  sigma <- matrix(runif(N * k, a, b), nrow = N)

  # initial pi
  pi_ <- rdirichlet(N, rep(delta, k))
  # initial w
  w <- rep(1 / N, N)

  return(list(mu = mu, sigma = sigma, pi_ = pi_, w = w))
}


pf_normal_k <- function(y, N, k,
                        a_mu, b_mu, a_sigma, b_sigma, delta,
                        is_initial = TRUE, last=list(),
                        verbose = TRUE) {
  start_time <- Sys.time()
  # pre-data
  n <- length(y)

  # initialise
  if (is_initial)
  {
    mu <- t(apply(matrix(
      rtruncnorm(n = N * k, a = a_mu, b = b_mu, mean = 0, sd = 1000),
      nrow = N
    ), 1, sort))

    sigma <- matrix(
      runif(n = N * k, a_sigma, b_sigma),
      nrow = N
    )

    pi_ <- matrix(
      runif(n = N * k, 0 + low_threshold, 1 - low_threshold),
      nrow = N
    )
    pi_ <- sweep(pi_, 1, rowSums(pi_), FUN = "/")
    w <- rep(1 / N, N)
  } else {
    mu <- last$mu
    sigma <- last$sigma
    pi_ <- last$pi_
    w <- last$w
  }


  # start loop
  progress_v <- floor(n * c(1:10) * 0.1)
  for (t in c(1:n)) {
    # print info
    if (verbose) {
      if (t %in% progress_v) {
        tcat("----------------\n")
        cat("Iter:", t, "\n")
      }
    }

    # weight updating
    likelihood <- rep(0, N)
    for (l in c(1:k)) {
      likelihood <- likelihood +
        pi_[, l] * dnorm(y[t], mean = mu[, l], sd = sigma[, l])
    }

    w <- w * likelihood
    w <- w / sum(w)
  }
  end_time <- Sys.time()
  if (verbose) {
    cat("Time used: ", end_time - start_time, "\n")
  }
  return(list(mu = mu, sigma = sigma, pi_ = pi_, w = w))
}
