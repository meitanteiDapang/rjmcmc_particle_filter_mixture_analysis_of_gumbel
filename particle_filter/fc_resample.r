
fc_resample <- function(candidate_weights, N) {
  M <- length(candidate_weights)

  f <- function(c) {
    sum(pmin(c * candidate_weights, 1))
  }
  # print(candidate_weights)
  # solve for c
  nonzero <- candidate_weights > 0
  if (!any(nonzero)) {
    print(candidate_weights)
    stop("all weights are zero!")
  } 
  c_max <- min(N / min(candidate_weights[nonzero]) / M, 1.797693e+300)
  # cat("c_max: ",c_max, "\n")
  c_val <- uniroot(function(c) f(c) - N, lower = 0, upper = c_max)$root
  # cat("c_val:", c_val, "\n")

  # split particles into A (keep) and B (resample)
  idx_A <- which(c_val * candidate_weights > 1)
  idx_B <- which(c_val * candidate_weights <= 1)
  N_A <- length(idx_A)
  N_B <- N - N_A
  # cat("N_A:", N_A, "\n")

  # sample from group B using stratified resampling
  selected_B <- integer(0)
  if (N_B > 0 && length(idx_B) > 0) {
    w_B <- candidate_weights[idx_B]
    p_B <- w_B / sum(w_B)
    cum_p_B <- cumsum(p_B)
    u <- (0:(N_B - 1) + runif(N_B)) / N_B
    selected_B <- idx_B[findInterval(u, cum_p_B) + 1]
  }

  # merge selected A and B
  final_idx <- c(idx_A, selected_B)

  # assign new weights
  # new_weights <- c(candidate_weights[idx_A], candidate_weights[selected_B])
  new_weights <- c(candidate_weights[idx_A], rep(1 / c_val, length(selected_B)))
  # new_weights <- rep(1/N, N)
  new_weights <- new_weights / sum(new_weights)

  return(list(indices = final_idx, weights = new_weights))
}


