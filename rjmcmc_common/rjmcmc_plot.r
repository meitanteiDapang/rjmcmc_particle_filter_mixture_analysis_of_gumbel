source(file.path(dir_path, "rjmcmc_common/utilities.r"))

plot_trace_and_density <- function(est_k3, burn.in, iter,
                                   param_type = c("mu", "tau", "pi"),
                                   k = NA,
                                   truth = NULL) {
  total_col <- ncol(est_k3)
  if (is.na(k)) {
    if (total_col == 2) {
      k <- 1
    } else {
      k <- as.integer((total_col + 1) / 3)
    }
  }

  mu_idx <- 1:k
  tau_idx <- (k + 1):(2 * k)
  pi_idx <- (2 * k + 1):(3 * k - 1)

  post <- as.data.frame(est_k3[burn.in:iter, ])
  colnames(post)[mu_idx] <- paste0("mu", 1:k)
  colnames(post)[tau_idx] <- paste0("tau", 1:k)

  if (param_type == "pi") {
    if (length(pi_idx) > 0 && max(pi_idx) <= ncol(post)) {
      colnames(post)[pi_idx] <- paste0("pi", 1:(k - 1))
      post[[paste0("pi", k)]] <- 1 - rowSums(post[, paste0("pi", 1:(k - 1)), drop = FALSE])
    } else if (k == 1) {
      post[["pi1"]] <- 1
    }
  }

  post$iter <- seq_len(nrow(post))

  if (param_type == "mu") {
    param_names <- paste0("mu", 1:k)
    greek_labels <- sapply(1:k, function(j) sprintf("mu[%d]", j))
  } else if (param_type == "tau") {
    param_names <- paste0("tau", 1:k)
    greek_labels <- sapply(1:k, function(j) sprintf("tau[%d]", j))
  } else if (param_type == "pi") {
    param_names <- paste0("pi", 1:k)
    greek_labels <- sapply(1:k, function(j) sprintf("pi[%d]", j))
  } else {
    stop("param_type must be one of 'mu', 'tau', or 'pi'.")
  }
  names(greek_labels) <- param_names

  df_long <- post %>%
    dplyr::select(iter, dplyr::all_of(param_names)) %>%
    tidyr::pivot_longer(-iter, names_to = "Parameter", values_to = "Value")

  # Trace plot
  p_trace <- ggplot2::ggplot(df_long, ggplot2::aes(x = iter, y = Value, color = Parameter)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Iteration", y = NULL) +
    ggplot2::scale_color_discrete(labels = parse(text = greek_labels)) +
    ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 12))

  if (!is.null(truth) && length(truth) > 0) {
    for (v in truth) {
      p_trace <- p_trace + ggplot2::geom_hline(
        yintercept = v, linetype = "dashed",
        color = "black", linewidth = 0.4
      )
    }
  }
  print(p_trace)

  # Density plot
  p_density <- ggplot2::ggplot(df_long, ggplot2::aes(x = Value, color = Parameter, fill = Parameter)) +
    ggplot2::geom_density(alpha = 0.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_color_discrete(labels = parse(text = greek_labels)) +
    ggplot2::scale_fill_discrete(labels = parse(text = greek_labels)) +
    ggplot2::guides(color = ggplot2::guide_legend(title = NULL), fill = "none") +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 12))

  if (!is.null(truth) && length(truth) > 0) {
    for (v in truth) {
      p_density <- p_density + ggplot2::geom_vline(
        xintercept = v, linetype = "dashed",
        color = "black", linewidth = 0.4
      )
    }
  }
  print(p_density)
}


plot_sixpanel_trace_density <- function(est_k3, burn.in, iter,
                                        k = NA,
                                        truths = list(mu = NULL, tau = NULL, pi = NULL),
                                        trace_alpha = 0.45,
                                        legend_position = "bottom",
                                        panel_label_size = 10, # smaller a-f labels
                                        mu_range = NULL,
                                        tau_range = NULL,
                                        pi_range = NULL) {
  # ---- Infer k & column indices (mu 1..k, tau 1..k, pi 1..k-1) ----
  total_col <- ncol(est_k3)
  if (is.na(k)) {
    if (total_col == 2) k <- 1 else k <- as.integer((total_col + 1) / 3)
  }
  mu_idx <- 1:k
  tau_idx <- (k + 1):(2 * k)
  pi_idx <- (2 * k + 1):(3 * k - 1) # k-1 columns, last one will be completed to sum to 1

  post <- as.data.frame(est_k3[burn.in:iter, ])
  colnames(post)[mu_idx] <- paste0("mu", 1:k)
  colnames(post)[tau_idx] <- paste0("tau", 1:k)

  # ---- Build long data for each parameter; unify colour mapping by component (1..k) ----
  build_param_df <- function(param) {
    if (param == "mu") {
      df <- post |>
        dplyr::select(dplyr::all_of(paste0("mu", 1:k))) |>
        dplyr::mutate(iter = seq_len(nrow(post))) |>
        tidyr::pivot_longer(
          cols = starts_with("mu"),
          names_to = "name", values_to = "value"
        ) |>
        dplyr::mutate(
          comp = as.integer(sub("mu", "", name)),
          param = "mu"
        )
    } else if (param == "tau") {
      df <- post |>
        dplyr::select(dplyr::all_of(paste0("tau", 1:k))) |>
        dplyr::mutate(iter = seq_len(nrow(post))) |>
        tidyr::pivot_longer(
          cols = starts_with("tau"),
          names_to = "name", values_to = "value"
        ) |>
        dplyr::mutate(
          comp = as.integer(sub("tau", "", name)),
          param = "tau"
        )
    } else { # "pi"
      # Complete pi_k so that sum_j pi_j = 1
      if (length(pi_idx) > 0 && max(pi_idx) <= ncol(post)) {
        colnames(post)[pi_idx] <- paste0("pi", 1:(k - 1))
        post[[paste0("pi", k)]] <- 1 - rowSums(post[, paste0("pi", 1:(k - 1)), drop = FALSE])
      } else if (k == 1 && is.null(post[["pi1"]])) {
        post[["pi1"]] <- 1
      }
      df <- post |>
        dplyr::select(dplyr::all_of(paste0("pi", 1:k))) |>
        dplyr::mutate(iter = seq_len(nrow(post))) |>
        tidyr::pivot_longer(
          cols = starts_with("pi"),
          names_to = "name", values_to = "value"
        ) |>
        dplyr::mutate(
          comp = as.integer(sub("pi", "", name)),
          param = "pi"
        )
    }
    df |>
      dplyr::select(iter, value, comp, param) |>
      dplyr::mutate(comp = factor(comp, levels = 1:k))
  }

  df_mu <- build_param_df("mu")
  df_tau <- build_param_df("tau")
  df_pi <- build_param_df("pi")

  # ---- Shared colours for components 1..k (single legend) ----
  base_cols <- scales::hue_pal()(k)
  names(base_cols) <- levels(df_mu$comp)

  # ---- Helpers: add all truth values regardless of length ----
  add_truth_hlines <- function(p, truth_vec) {
    if (!is.null(truth_vec) && length(truth_vec) > 0) {
      for (val in truth_vec) {
        if (!is.na(val)) {
          p <- p + geom_hline(
            yintercept = val, linetype = "dashed",
            colour = "black", linewidth = 0.4
          )
        }
      }
    }
    p
  }

  add_truth_vlines <- function(p, truth_vec) {
    if (!is.null(truth_vec) && length(truth_vec) > 0) {
      for (val in truth_vec) {
        if (!is.na(val)) {
          p <- p + geom_vline(
            xintercept = val, linetype = "dashed",
            colour = "black", linewidth = 0.4
          )
        }
      }
    }
    p
  }


  # ---- Plot builders with Greek axis labels ----
  make_trace <- function(df, y_symbol, truth_vec = NULL, y_range = NULL) {
    gg <- ggplot(df, aes(x = iter, y = value, colour = comp)) +
      geom_line(alpha = trace_alpha, linewidth = 0.5) +
      theme_minimal(base_size = 12) +
      labs(x = "Iteration", y = y_symbol) + # y in Greek
      scale_colour_manual(values = base_cols, name = "Component") +
      theme(legend.position = "none")
    if (!is.null(truth_vec)) {
      gg <- add_truth_hlines(gg, truth_vec)
    }
    if (!is.null(y_range)) {
      gg <- gg + coord_cartesian(ylim = y_range)
    }
    gg
  }
  make_dens <- function(df, x_symbol, truth_vec = NULL, x_range = NULL) {
    gg <- ggplot(df, aes(x = value, colour = comp, fill = comp)) +
      geom_density(alpha = 0.2) +
      theme_minimal(base_size = 12) +
      labs(x = x_symbol, y = NULL) + # x in Greek
      scale_colour_manual(values = base_cols, name = "Component") +
      scale_fill_manual(values = base_cols, guide = "none") +
      theme(legend.position = "none")
    if (!is.null(truth_vec)) {
      gg <- add_truth_vlines(gg, truth_vec)
    }
    if (!is.null(x_range)) {
      gg <- gg + coord_cartesian(xlim = x_range)
    }
    gg
  }

  # ---- Six panels with Greek labels and subheadings ----
  p_a <- make_trace(df_mu, expression(mu), truths$mu, mu_range) +
    labs(title = expression("(a) traceplot of " * underline(mu))) +
    theme(plot.title = element_text(hjust = 0.5, size = panel_label_size))

  p_b <- make_dens(df_mu, expression(mu), truths$mu, mu_range) +
    labs(title = expression("(b) density of " * underline(mu))) +
    theme(plot.title = element_text(hjust = 0.5, size = panel_label_size))

  p_c <- make_trace(df_tau, expression(tau), truths$tau, tau_range) +
    labs(title = expression("(c) traceplot of " * underline(tau))) +
    theme(plot.title = element_text(hjust = 0.5, size = panel_label_size))

  p_d <- make_dens(df_tau, expression(tau), truths$tau, tau_range) +
    labs(title = expression("(d) density of " * underline(tau))) +
    theme(plot.title = element_text(hjust = 0.5, size = panel_label_size))

  p_e <- make_trace(df_pi, expression(pi), truths$pi, pi_range) +
    labs(title = expression("(e) traceplot of " * underline(pi))) +
    theme(plot.title = element_text(hjust = 0.5, size = panel_label_size))

  p_f <- make_dens(df_pi, expression(pi), truths$pi, pi_range) +
    labs(title = expression("(f) density of " * underline(pi))) +
    theme(plot.title = element_text(hjust = 0.5, size = panel_label_size))


  # ---- Single legend (from a dummy plot) ----
  p_for_legend <- ggplot(df_mu, aes(x = iter, y = value, colour = comp)) +
    geom_line(alpha = trace_alpha, linewidth = 0.5) +
    theme_minimal(base_size = 12) +
    labs(x = "Iteration", y = expression(mu)) +
    scale_colour_manual(values = base_cols, name = "Component") +
    theme(legend.position = legend_position)
  legend_g <- cowplot::get_legend(p_for_legend)


  grid_6 <- cowplot::plot_grid(
    p_a, p_b, p_c, p_d, p_e, p_f,
    ncol = 2, align = "hv"
  )

  # ---- Attach legend once ----
  if (legend_position %in% c("bottom", "top")) {
    rel_heights <- if (legend_position == "bottom") c(1, 0.12) else c(0.12, 1)
    out <- cowplot::plot_grid(
      if (legend_position == "top") legend_g else grid_6,
      if (legend_position == "top") grid_6 else legend_g,
      ncol = 1, rel_heights = rel_heights
    )
  } else if (legend_position %in% c("left", "right")) {
    rel_widths <- if (legend_position == "left") c(0.22, 1) else c(1, 0.22)
    out <- cowplot::plot_grid(
      if (legend_position == "left") legend_g else grid_6,
      if (legend_position == "left") grid_6 else legend_g,
      ncol = 2, rel_widths = rel_widths
    )
  } else {
    out <- grid_6
  }

  print(out)
  invisible(out)
  return(out)
}


plot_rjmcmc_trace_and_density <- function(mat, param_type = c("mu", "tau", "pi"), k,
                                          true_values = NULL) {
  param_type <- match.arg(param_type)

  mat <- as.data.frame(mat)
  n <- nrow(mat)

  # Only take first k columns
  if (ncol(mat) < k) stop("Input matrix does not have enough columns for the chosen k.")
  mat <- mat[, 1:k, drop = FALSE]
  colnames(mat) <- paste0(param_type, 1:k)

  # Greek legend labels
  greek_labels <- sapply(1:k, function(j) sprintf("%s[%d]", param_type, j))
  names(greek_labels) <- paste0(param_type, 1:k)

  # Add iteration index
  mat$iter <- seq_len(n)

  # Reshape to long format
  df_long <- mat %>%
    pivot_longer(-iter, names_to = "Parameter", values_to = "Value")

  # Traceplot
  p_trace <- ggplot(df_long, aes(x = iter, y = Value, color = Parameter)) +
    geom_line() +
    theme_minimal() +
    labs(x = "Iteration", y = NULL) +
    scale_color_discrete(labels = parse(text = greek_labels)) +
    guides(color = guide_legend(title = NULL)) +
    theme(legend.text = element_text(size = 12))

  # Add true values as black dashed lines (if provided)
  if (!is.null(true_values)) {
    true_vals <- rep(NA, k)
    true_vals[seq_along(true_values)] <- true_values
    for (j in seq_len(k)) {
      if (!is.na(true_vals[j])) {
        p_trace <- p_trace +
          geom_hline(yintercept = true_vals[j], linetype = "dashed", color = "black", linewidth = 0.7)
      }
    }
  }
  print(p_trace)

  # Density plot (no duplicate legend)
  p_density <- ggplot(df_long, aes(x = Value, color = Parameter, fill = Parameter)) +
    geom_density(alpha = 0.2) +
    theme_minimal() +
    labs(x = NULL, y = NULL) +
    scale_color_discrete(labels = parse(text = greek_labels)) +
    scale_fill_discrete(labels = parse(text = greek_labels)) +
    guides(
      color = guide_legend(title = NULL),
      fill = "none"
    ) +
    theme(legend.text = element_text(size = 12))

  # Add true values as black vertical dashed lines
  if (!is.null(true_values)) {
    true_vals <- rep(NA, k)
    true_vals[seq_along(true_values)] <- true_values
    for (j in seq_len(k)) {
      if (!is.na(true_vals[j])) {
        p_density <- p_density +
          geom_vline(xintercept = true_vals[j], linetype = "dashed", color = "black", linewidth = 0.7)
      }
    }
  }
  print(p_density)
}


# Function: Plot coverage for two sets of results
plot_coverage <- function(df, y_label, true_value) {
  df <- df %>%
    arrange(m) %>%
    mutate(x = row_number())

  p <- ggplot(df, aes(x = x, y = m)) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
      width = 0.1, color = "red"
    ) +
    geom_point(color = "blue", size = 2) +
    geom_hline(yintercept = true_value, color = "black") +
    labs(x = "Simulations", y = y_label) +
    theme_minimal() +
    theme(legend.position = "none")

  print(p)
}


# ------------------------------------------------------------------------------
# Six-panel (a-f) plot for RJMCMC outputs: for mu, tau, pi
# Panels: (a) mu trace, (b) mu density, (c) tau trace, (d) tau density,
#         (e) pi trace, (f) pi density
# - Single shared legend (components 1..k)
# - Trace lines with alpha
# - Greek axis labels
# - All truth values are drawn (no truncation)
# ------------------------------------------------------------------------------
plot_rjmcmc_sixpanel <- function(mats,
                                 k,
                                 truths = list(mu = NULL, tau = NULL, pi = NULL),
                                 trace_alpha = 0.45,
                                 legend_position = "bottom",
                                 panel_label_size = 10,
                                 mu_range = NULL,
                                 tau_range = NULL,
                                 pi_range = NULL) {
  # ---- Validate inputs ----
  stopifnot(is.list(mats), all(c("mu", "tau", "pi") %in% names(mats)))
  Xmu <- as.data.frame(mats$mu)
  Xtau <- as.data.frame(mats$tau)
  Xpi <- as.data.frame(mats$pi)

  n_mu <- nrow(Xmu)
  n_tau <- nrow(Xtau)
  n_pi <- nrow(Xpi)
  if (length(unique(c(n_mu, n_tau, n_pi))) != 1) {
    stop("mu, tau, and pi matrices must have the same number of rows (iterations).")
  }

  # ---- Ensure first k columns exist; for pi, complete to k if only k-1 are present ----
  if (ncol(Xmu) < k || ncol(Xtau) < k) stop("mu/tau do not have at least k columns.")
  Xmu <- Xmu[, 1:k, drop = FALSE]
  Xtau <- Xtau[, 1:k, drop = FALSE]

  if (ncol(Xpi) >= k) {
    Xpi <- Xpi[, 1:k, drop = FALSE]
  } else if (ncol(Xpi) == (k - 1)) {
    pi_last <- 1 - rowSums(Xpi[, 1:(k - 1), drop = FALSE])
    Xpi <- cbind(Xpi[, 1:(k - 1), drop = FALSE], pi_last)
    colnames(Xpi) <- paste0("V", 1:k)
  } else {
    stop("pi must have k or k-1 columns.")
  }

  # ---- Name columns by parameter and component ----
  colnames(Xmu) <- paste0("mu", 1:k)
  colnames(Xtau) <- paste0("tau", 1:k)
  colnames(Xpi) <- paste0("pi", 1:k)

  # ---- Long data builders (shared 'comp' factor 1..k) ----
  build_long <- function(X, prefix) {
    X$iter <- seq_len(nrow(X))
    X %>%
      tidyr::pivot_longer(-iter, names_to = "name", values_to = "value") %>%
      dplyr::mutate(
        comp = factor(as.integer(sub(prefix, "", name)), levels = 1:k),
        param = prefix
      ) %>%
      dplyr::select(iter, value, comp, param)
  }

  df_mu <- build_long(Xmu, "mu")
  df_tau <- build_long(Xtau, "tau")
  df_pi <- build_long(Xpi, "pi")

  # ---- Shared colours for components 1..k ----
  base_cols <- scales::hue_pal()(k)
  names(base_cols) <- levels(df_mu$comp)

  # ---- Truth lines ----
  add_truth_hlines <- function(p, truth_vec) {
    if (!is.null(truth_vec) && length(truth_vec) > 0) {
      for (val in truth_vec) {
        if (!is.na(val)) {
          p <- p + ggplot2::geom_hline(
            yintercept = val, linetype = "dashed",
            colour = "black", linewidth = 0.4
          )
        }
      }
    }
    p
  }
  add_truth_vlines <- function(p, truth_vec) {
    if (!is.null(truth_vec) && length(truth_vec) > 0) {
      for (val in truth_vec) {
        if (!is.na(val)) {
          p <- p + ggplot2::geom_vline(
            xintercept = val, linetype = "dashed",
            colour = "black", linewidth = 0.4
          )
        }
      }
    }
    p
  }

  # ---- Plot constructors ----
  make_trace <- function(df, y_expr, truth_vec = NULL, y_range = NULL) {
    g <- ggplot2::ggplot(df, ggplot2::aes(x = iter, y = value, colour = comp)) +
      ggplot2::geom_line(alpha = trace_alpha, linewidth = 0.5) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(x = "Iteration", y = y_expr) +
      ggplot2::scale_colour_manual(values = base_cols, name = "Component") +
      ggplot2::theme(legend.position = "none")
    if (!is.null(y_range)) {
      g <- g + ggplot2::coord_cartesian(ylim = y_range)
    }
    add_truth_hlines(g, truth_vec)
  }
  make_density <- function(df, x_expr, truth_vec = NULL, x_range = NULL) {
    g <- ggplot2::ggplot(df, ggplot2::aes(x = value, colour = comp, fill = comp)) +
      ggplot2::geom_density(alpha = 0.2) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(x = x_expr, y = NULL) +
      ggplot2::scale_colour_manual(values = base_cols, name = "Component") +
      ggplot2::scale_fill_manual(values = base_cols, guide = "none") +
      ggplot2::theme(legend.position = "none")
    if (!is.null(x_range)) {
      g <- g + ggplot2::coord_cartesian(xlim = x_range)
    }
    add_truth_vlines(g, truth_vec)
  }

  # ---- Six panels + centered subheadings with expression(underline(.)) ----
  p_a <- make_trace(df_mu, expression(mu), truths$mu, mu_range) +
    ggplot2::labs(title = expression("(a) traceplot of " * underline(mu))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  p_b <- make_density(df_mu, expression(mu), truths$mu, mu_range) +
    ggplot2::labs(title = expression("(b) density of " * underline(mu))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  p_c <- make_trace(df_tau, expression(tau), truths$tau, tau_range) +
    ggplot2::labs(title = expression("(c) traceplot of " * underline(tau))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  p_d <- make_density(df_tau, expression(tau), truths$tau, tau_range) +
    ggplot2::labs(title = expression("(d) density of " * underline(tau))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  p_e <- make_trace(df_pi, expression(pi), truths$pi, pi_range) +
    ggplot2::labs(title = expression("(e) traceplot of " * underline(pi))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  p_f <- make_density(df_pi, expression(pi), truths$pi, pi_range) +
    ggplot2::labs(title = expression("(f) density of " * underline(pi))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  # ---- Single legend (extract once) ----
  p_for_legend <- ggplot2::ggplot(df_mu, ggplot2::aes(x = iter, y = value, colour = comp)) +
    ggplot2::geom_line(alpha = trace_alpha, linewidth = 0.5) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Iteration", y = expression(mu)) +
    ggplot2::scale_colour_manual(values = base_cols, name = "Component") +
    ggplot2::theme(legend.position = legend_position)
  legend_g <- cowplot::get_legend(p_for_legend)

  # ---- Combine panels: 2 cols * 3 rows (no a-f labels; titles already added) ----
  grid_6 <- cowplot::plot_grid(
    p_a, p_b, p_c, p_d, p_e, p_f,
    ncol = 2, align = "hv"
  )

  # ---- Attach legend once ----
  out <- switch(legend_position,
    "bottom" = cowplot::plot_grid(grid_6, legend_g, ncol = 1, rel_heights = c(1, 0.12)),
    "top"    = cowplot::plot_grid(legend_g, grid_6, ncol = 1, rel_heights = c(0.12, 1)),
    "left"   = cowplot::plot_grid(legend_g, grid_6, ncol = 2, rel_widths = c(0.22, 1)),
    "right"  = cowplot::plot_grid(grid_6, legend_g, ncol = 2, rel_widths = c(1, 0.22)),
    grid_6
  )

  print(out)
  invisible(out)
  return(out)
}


plot_rjmcmc_sixpanel_case_specific <- function(
  mats,
  k,
  truths = list(mu = NULL, tau = NULL, pi = NULL),
  trace_alpha = 0.45,
  legend_position = "bottom",
  panel_label_size = 10,
  # ---- NEW: axis ranges (use coord_cartesian to zoom without dropping data) ----
  iter_xlim = NULL, # e.g., c(1001, 5000)
  mu_trace_ylim = NULL, # e.g., c(-5, 5)
  tau_trace_ylim = NULL, # e.g., c(0, 10)
  pi_trace_ylim = c(0, 1), # sensible default for probabilities
  mu_dens_xlim = NULL, # e.g., c(-5, 5)
  tau_dens_xlim = NULL, # e.g., c(0, 10)
  pi_dens_xlim = c(0, 1) # sensible default for probabilities
) {
  # ---- Validate inputs ----
  stopifnot(is.list(mats), all(c("mu", "tau", "pi") %in% names(mats)))
  Xmu <- as.data.frame(mats$mu)
  Xtau <- as.data.frame(mats$tau)
  Xpi <- as.data.frame(mats$pi)

  n_mu <- nrow(Xmu)
  n_tau <- nrow(Xtau)
  n_pi <- nrow(Xpi)
  if (length(unique(c(n_mu, n_tau, n_pi))) != 1) {
    stop("mu, tau, and pi matrices must have the same number of rows (iterations).")
  }

  # ---- Ensure first k columns exist; for pi, complete to k if only k-1 are present ----
  if (ncol(Xmu) < k || ncol(Xtau) < k) stop("mu/tau do not have at least k columns.")
  Xmu <- Xmu[, 1:k, drop = FALSE]
  Xtau <- Xtau[, 1:k, drop = FALSE]

  if (ncol(Xpi) >= k) {
    Xpi <- Xpi[, 1:k, drop = FALSE]
  } else if (ncol(Xpi) == (k - 1)) {
    pi_last <- 1 - rowSums(Xpi[, 1:(k - 1), drop = FALSE])
    Xpi <- cbind(Xpi[, 1:(k - 1), drop = FALSE], pi_last)
    colnames(Xpi) <- paste0("V", 1:k)
  } else {
    stop("pi must have k or k-1 columns.")
  }

  # ---- Name columns by parameter and component ----
  colnames(Xmu) <- paste0("mu", 1:k)
  colnames(Xtau) <- paste0("tau", 1:k)
  colnames(Xpi) <- paste0("pi", 1:k)

  # ---- Long data builders (shared 'comp' factor 1..k) ----
  build_long <- function(X, prefix) {
    X$iter <- seq_len(nrow(X))
    X %>%
      tidyr::pivot_longer(-iter, names_to = "name", values_to = "value") %>%
      dplyr::mutate(
        comp  = factor(as.integer(sub(prefix, "", name)), levels = 1:k),
        param = prefix
      ) %>%
      dplyr::select(iter, value, comp, param)
  }

  df_mu <- build_long(Xmu, "mu")
  df_tau <- build_long(Xtau, "tau")
  df_pi <- build_long(Xpi, "pi")

  # ---- Shared colours for components 1..k ----
  base_cols <- scales::hue_pal()(k)
  names(base_cols) <- levels(df_mu$comp)

  # ---- Truth lines ----
  add_truth_hlines <- function(p, truth_vec) {
    if (!is.null(truth_vec) && length(truth_vec) > 0) {
      for (val in truth_vec) {
        if (!is.na(val)) {
          p <- p + ggplot2::geom_hline(
            yintercept = val, linetype = "dashed",
            colour = "black", linewidth = 0.4
          )
        }
      }
    }
    p
  }
  add_truth_vlines <- function(p, truth_vec) {
    if (!is.null(truth_vec) && length(truth_vec) > 0) {
      for (val in truth_vec) {
        if (!is.na(val)) {
          p <- p + ggplot2::geom_vline(
            xintercept = val, linetype = "dashed",
            colour = "black", linewidth = 0.4
          )
        }
      }
    }
    p
  }

  # ---- Plot constructors (with coord_cartesian limits) ----
  make_trace <- function(df, y_expr, truth_vec = NULL, ylim = NULL) {
    g <- ggplot2::ggplot(df, ggplot2::aes(x = iter, y = value, colour = comp)) +
      ggplot2::geom_line(alpha = trace_alpha, linewidth = 0.5) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(x = "Iteration", y = y_expr) +
      ggplot2::scale_colour_manual(values = base_cols, name = "Component") +
      ggplot2::theme(legend.position = "none") +
      ggplot2::coord_cartesian(xlim = iter_xlim, ylim = ylim)
    add_truth_hlines(g, truth_vec)
  }

  make_density <- function(df, x_expr, truth_vec = NULL, xlim = NULL) {
    if (!is.null(xlim)) {
      df <- dplyr::filter(df, value >= xlim[1], value <= xlim[2])
    }

    g <- ggplot2::ggplot(df, ggplot2::aes(x = value, colour = comp, fill = comp)) +
      ggplot2::geom_density(alpha = 0.2) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(x = x_expr, y = NULL) +
      ggplot2::scale_colour_manual(values = base_cols, name = "Component") +
      ggplot2::scale_fill_manual(values = base_cols, guide = "none") +
      ggplot2::theme(legend.position = "none")

    add_truth_vlines(g, truth_vec)
  }

  # ---- Six panels + centred subheadings with expression(underline(.)) ----
  p_a <- make_trace(df_mu, expression(mu), truths$mu, ylim = mu_trace_ylim) +
    ggplot2::labs(title = expression("(a) traceplot of " * underline(mu))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  p_b <- make_density(df_mu, expression(mu), truths$mu, xlim = mu_dens_xlim) +
    ggplot2::labs(title = expression("(b) density of " * underline(mu))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  p_c <- make_trace(df_tau, expression(tau), truths$tau, ylim = tau_trace_ylim) +
    ggplot2::labs(title = expression("(c) traceplot of " * underline(tau))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  p_d <- make_density(df_tau, expression(tau), truths$tau, xlim = tau_dens_xlim) +
    ggplot2::labs(title = expression("(d) density of " * underline(tau))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  p_e <- make_trace(df_pi, expression(pi), truths$pi, ylim = pi_trace_ylim) +
    ggplot2::labs(title = expression("(e) traceplot of " * underline(pi))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  p_f <- make_density(df_pi, expression(pi), truths$pi, xlim = pi_dens_xlim) +
    ggplot2::labs(title = expression("(f) density of " * underline(pi))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = panel_label_size))

  # ---- Single legend (extract once) ----
  p_for_legend <- ggplot2::ggplot(df_mu, ggplot2::aes(x = iter, y = value, colour = comp)) +
    ggplot2::geom_line(alpha = trace_alpha, linewidth = 0.5) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Iteration", y = expression(mu)) +
    ggplot2::scale_colour_manual(values = base_cols, name = "Component") +
    ggplot2::theme(legend.position = legend_position)
  legend_g <- cowplot::get_legend(p_for_legend)

  # ---- Combine panels: 2 cols * 3 rows (no a–f labels; titles already added) ----
  grid_6 <- cowplot::plot_grid(
    p_a, p_b, p_c, p_d, p_e, p_f,
    ncol = 2, align = "hv"
  )

  # ---- Attach legend once ----
  out <- switch(legend_position,
    "bottom" = cowplot::plot_grid(grid_6, legend_g, ncol = 1, rel_heights = c(1, 0.12)),
    "top"    = cowplot::plot_grid(legend_g, grid_6, ncol = 1, rel_heights = c(0.12, 1)),
    "left"   = cowplot::plot_grid(legend_g, grid_6, ncol = 2, rel_widths = c(0.22, 1)),
    "right"  = cowplot::plot_grid(grid_6, legend_g, ncol = 2, rel_widths = c(1, 0.22)),
    grid_6
  )

  print(out)
  invisible(out)
  return(out)
}


plot_k_trace_and_density <- function(lst_est) {
  par(mfrow = c(1, 2))
  plot(lst_est$k, main = "k (Trace)", type = "l", col = "red")
  plot(density(lst_est$k), main = "k (Density)", type = "l", col = "red")
  par(mfrow = c(1, 1))

  if (exists("split_accept_count", envir = .GlobalEnv)) {
    cat("split_accept_count: ", get("split_accept_count", envir = .GlobalEnv), "/", get("split_count", envir = .GlobalEnv), "\n")
  }
  if (exists("combine_accept_count", envir = .GlobalEnv)) {
    cat("combine_accept_count: ", get("combine_accept_count", envir = .GlobalEnv), "/", get("combine_count", envir = .GlobalEnv), "\n")
  }
  if (exists("birth_accept_count", envir = .GlobalEnv)) {
    cat("birth_accept_count: ", get("birth_accept_count", envir = .GlobalEnv), "/", get("birth_count", envir = .GlobalEnv), "\n")
  }
  if (exists("death_accept_count", envir = .GlobalEnv)) {
    cat("death_accept_count: ", get("death_accept_count", envir = .GlobalEnv), "/", get("death_count", envir = .GlobalEnv), "\n")
  }

  k_table <- table(lst_est$k)
  print(k_table)
  print(k_table / sum(k_table))
  invisible(k_table)
}
