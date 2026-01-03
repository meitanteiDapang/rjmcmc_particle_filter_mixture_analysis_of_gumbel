merge_same <- function(particles, weights, toPrint = FALSE) {
  # Step 1: Convert particles to strings so we can detect duplicates
  particle_strings <- sapply(particles, function(p) paste(p, collapse = ","))

  # Step 2: Find unique particles
  unique_particles_strings <- unique(particle_strings)

  # Step 3: For each unique particle, sum weights
  uni_particles <- list()
  uni_w <- numeric(length(unique_particles_strings))

  for (i in seq_along(unique_particles_strings)) {
    idx <- which(particle_strings == unique_particles_strings[i])
    uni_particles[[i]] <- particles[[idx[1]]] # representative particle
    uni_w[i] <- sum(weights[idx])
  }

  # Now sort by weight
  sorted_indices <- order(uni_w, decreasing = TRUE)

  # Apply the sorting
  uni_particles <- uni_particles[sorted_indices]
  uni_w <- uni_w[sorted_indices]

  # For print
  if (toPrint) {
    particle_strings_out <-
      sapply(uni_particles, function(p) paste(p, collapse = ","))
    for (i in c(1:length(uni_particles))) {
      cat("w: ", uni_w[i], "     Particles: ", particle_strings_out[i], "\n")
    }
  }
  return(
    list(
      uni_particles = uni_particles,
      uni_w = uni_w
    )
  )
}


normalize_log_weights <- function(log_weights) {
  max_log <- max(log_weights)
  log_sum_exp <- max_log + log(sum(exp(log_weights - max_log)))
  normalized_weights <- exp(log_weights - log_sum_exp)
  return(normalized_weights)
}
get_particle_indices_by_k <- function(uni_particles, k) {
  # Indices of uni_particles whose cluster count equals k
  which(sapply(uni_particles, function(particle) {
    length(unique(particle)) == k
  }))
}


varx <- function(x) {
  if (length(x) > 1) {
    return(var(x))
  } else {
    return(0)
  }
}


plot_threepanel_density <- function(est_k3,
                                    k,
                                    truths = list(mu = NULL, tau = NULL, pi = NULL),
                                    panel_label_size = 10,
                                    density_alpha = 0.25) {
  stopifnot(is.matrix(est_k3) || is.data.frame(est_k3))
  est_k3 <- as.data.frame(est_k3)

  mu_idx  <- 1:k
  tau_idx <- (k + 1):(2 * k)
  pi_idx  <- (2 * k + 1):(3 * k - 1)

  colnames(est_k3)[mu_idx]  <- paste0("mu",  1:k)
  colnames(est_k3)[tau_idx] <- paste0("tau", 1:k)

  if (k == 1) {
    if (is.null(est_k3[["pi1"]])) est_k3[["pi1"]] <- 1
  } else {
    if (length(pi_idx) > 0 && max(pi_idx) <= ncol(est_k3)) {
      colnames(est_k3)[pi_idx] <- paste0("pi", 1:(k - 1))
      est_k3[[paste0("pi", k)]] <- 1 - rowSums(est_k3[, paste0("pi", 1:(k - 1)), drop = FALSE])
    } else {
      for (j in 1:(k - 1)) est_k3[[paste0("pi", j)]] <- 1 / k
      est_k3[[paste0("pi", k)]] <- 1 - rowSums(est_k3[, paste0("pi", 1:(k - 1)), drop = FALSE])
    }
  }

  build_param_df <- function(prefix) {
    cols <- paste0(prefix, 1:k)
    df <- est_k3[, cols, drop = FALSE]
    df |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      tidyr::pivot_longer(cols = dplyr::all_of(cols),
                          names_to = "name", values_to = "value") |>
      dplyr::mutate(comp = factor(as.integer(sub(prefix, "", name)), levels = 1:k),
                    param = prefix) |>
      dplyr::select(value, comp, param)
  }

  df_mu  <- build_param_df("mu")
  df_tau <- build_param_df("tau")
  df_pi  <- build_param_df("pi")

  base_cols <- scales::hue_pal()(k)
  names(base_cols) <- as.character(1:k)

  add_truth_vlines <- function(p, truth_vec) {
    if (!is.null(truth_vec) && length(truth_vec) > 0) {
      for (val in truth_vec) if (!is.na(val)) {
        p <- p + ggplot2::geom_vline(xintercept = val, linetype = "dashed",
                                     colour = "black", linewidth = 0.4)
      }
    }
    p
  }

  make_dens <- function(df, x_symbol, truth_vec = NULL, panel_title) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = value, colour = comp, fill = comp)) +
      ggplot2::geom_density(alpha = density_alpha) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(x = x_symbol, y = NULL, title = panel_title) +
      ggplot2::scale_colour_manual(values = base_cols, name = "Component") +
      ggplot2::scale_fill_manual(values = base_cols, guide = "none") +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))
    if (!is.null(truth_vec)) {
        p <- add_truth_vlines(p, truth_vec)
    }
    p
  }

  # panels with subheadings
  p_mu  <- make_dens(df_mu,  expression(mu),  truths$mu,
                     expression("(a) density of " * underline(mu)))
  p_tau <- make_dens(df_tau, expression(tau), truths$tau,
                     expression("(b) density of " * underline(tau)))
  p_pi  <- make_dens(df_pi,  expression(pi),  truths$pi,
                     expression("(c) density of " * underline(pi)))

  # legend
  p_for_legend <- ggplot2::ggplot(df_mu, ggplot2::aes(x = value, colour = comp, fill = comp)) +
    ggplot2::geom_density(alpha = density_alpha) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = expression(mu), y = NULL) +
    ggplot2::scale_colour_manual(values = base_cols, name = "Component") +
    ggplot2::scale_fill_manual(values = base_cols, guide = "none") +
    ggplot2::theme(legend.position = "bottom")
  legend_g <- cowplot::get_legend(p_for_legend)

  # 3 panels stacked vertically
  grid_3 <- cowplot::plot_grid(
    p_mu, p_tau, p_pi,
    ncol = 1, align = "v"
  )

  out <- cowplot::plot_grid(
    grid_3, legend_g,
    ncol = 1, rel_heights = c(1, 0.12)
  )

  print(out)
  invisible(out)
  out
}



plot_twopanel_density <- function(est_k3,
                                  k,
                                  truths = list(lambda = NULL, pi = NULL),
                                  panel_label_size = 10,
                                  density_alpha = 0.25) {
  stopifnot(is.matrix(est_k3) || is.data.frame(est_k3))
  est_k3 <- as.data.frame(est_k3)

  lambda_idx <- 1:k
  pi_idx     <- (k + 1):(2 * k - 1)

  # ---- rename / complete columns ----
  if (!all(paste0("lambda", 1:k) %in% colnames(est_k3))) {
    if (max(lambda_idx) <= ncol(est_k3)) {
      colnames(est_k3)[lambda_idx] <- paste0("lambda", 1:k)
    }
  }

  if (k == 1) {
    if (is.null(est_k3[["pi1"]])) est_k3[["pi1"]] <- 1
  } else {
    have_pi <- any(startsWith(colnames(est_k3), "pi"))
    if (!have_pi && length(pi_idx) > 0 && max(pi_idx) <= ncol(est_k3)) {
      colnames(est_k3)[pi_idx] <- paste0("pi", 1:(k - 1))
    }
    for (j in 1:(k - 1)) if (is.null(est_k3[[paste0("pi", j)]])) est_k3[[paste0("pi", j)]] <- 1 / k
    est_k3[[paste0("pi", k)]] <- 1 - rowSums(est_k3[, paste0("pi", 1:(k - 1)), drop = FALSE])
  }

  # ---- long dfs ----
  build_param_df <- function(prefix) {
    cols <- paste0(prefix, 1:k)
    missing_cols <- setdiff(cols, colnames(est_k3))
    if (length(missing_cols)) {
      stop(sprintf("Missing columns for prefix '%s': %s", prefix, paste(missing_cols, collapse = ", ")))
    }
    est_k3 |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      tidyr::pivot_longer(cols = dplyr::all_of(cols), names_to = "name", values_to = "value") |>
      dplyr::mutate(comp = factor(as.integer(sub(prefix, "", name)), levels = 1:k),
                    param = prefix) |>
      dplyr::select(value, comp, param)
  }

  df_lambda <- build_param_df("lambda")
  df_pi     <- build_param_df("pi")

  base_cols <- scales::hue_pal()(k)
  names(base_cols) <- as.character(1:k)

  add_truth_vline <- function(p, val) {
    if (!is.null(val) && !is.na(val)) {
      p <- p + ggplot2::geom_vline(xintercept = val, linetype = "dashed",
                                   colour = "black", linewidth = 0.4)
    }
    p
  }

  # ---- letter counter for subheadings ----
  lab_idx <- 0L
  next_lab <- function() { lab_idx <<- lab_idx + 1L; paste0("(", letters[lab_idx], ")") }

  # ---- lambda panels (each component its own filled density) ----
  make_lambda_panel <- function(j, lab_text) {
    df_j <- dplyr::filter(df_lambda, comp == j)
    p <- ggplot2::ggplot(df_j, ggplot2::aes(x = value)) +
      ggplot2::geom_density(fill = base_cols[as.character(j)], alpha = density_alpha,
                            colour = base_cols[as.character(j)]) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        x = bquote(lambda[.(j)]),
        y = NULL,
        title = bquote(.(lab_text) * " density of " * underline(lambda[.(j)]))
      ) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size)
      )
    truth_val <- if (!is.null(truths$lambda) && length(truths$lambda) >= j) truths$lambda[j] else NA
    add_truth_vline(p, truth_val)
  }

  lambda_plots <- lapply(1:k, function(j) make_lambda_panel(j, next_lab()))

  # ---- pi panel (multi-component overlaid with shared palette) ----
  pi_lab <- next_lab()
  p_pi <- ggplot2::ggplot(df_pi, ggplot2::aes(x = value, colour = comp, fill = comp)) +
    ggplot2::geom_density(alpha = density_alpha) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x = expression(pi), y = NULL,
      title = bquote(.(pi_lab) * " density of " * underline(pi))
    ) +
    ggplot2::scale_colour_manual(values = base_cols, name = "Component") +
    ggplot2::scale_fill_manual(values = base_cols, guide = "none") +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size)
    )

  if (!is.null(truths$pi) && length(truths$pi) > 0) {
    for (j in seq_len(min(k, length(truths$pi)))) {
      val <- truths$pi[j]
      if (!is.na(val)) {
        p_pi <- p_pi + ggplot2::geom_vline(xintercept = val, linetype = "dashed",
                                           colour = "black", linewidth = 0.4)
      }
    }
  }

  # ---- legend (once) ----
  p_for_legend <- ggplot2::ggplot(df_pi, ggplot2::aes(x = value, colour = comp, fill = comp)) +
    ggplot2::geom_density(alpha = density_alpha) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = expression(pi), y = NULL) +
    ggplot2::scale_colour_manual(values = base_cols, name = "Component") +
    ggplot2::scale_fill_manual(values = base_cols, guide = "none") +
    ggplot2::theme(legend.position = "bottom")
  legend_g <- cowplot::get_legend(p_for_legend)

  # ---- arrange panels (no letter labels; titles already contain them) ----
  panels <- c(lambda_plots, list(p_pi))
  grid_panels <- cowplot::plot_grid(plotlist = panels, ncol = 2, align = "hv")

  out <- cowplot::plot_grid(grid_panels, legend_g, ncol = 1, rel_heights = c(1, 0.12))

  print(out)
  invisible(out)
  out
}
