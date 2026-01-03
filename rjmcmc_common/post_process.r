filter_by_component_size <- function(all_Z, lst_est1, min_count = 2) {
  n <- nrow(all_Z)
  # Check whether every component in each row appears at least min_count times
  row_good <- sapply(1:n, function(j) {
    tab <- table(factor(all_Z[j, ], levels = seq_len(lst_est1$k[j])))
    all(tab >= min_count)
  })

  # Mark rows to drop: if row_good[j] == FALSE, then drop j and j+1
  drop_rows <- rep(FALSE, n)
  for (j in 1:(n - 1)) {
    if (!row_good[j]) {
      drop_rows[j] <- TRUE
      drop_rows[j + 1] <- TRUE
    }
  }
  # If last row is not good, also drop it
  if (!row_good[n]) drop_rows[n] <- TRUE

  keep_rows <- !drop_rows

  # Filter all parameter matrices/vectors
  lst_est_after <- list()
  lst_est_after$mu  <- lst_est1$mu[keep_rows, , drop = FALSE]
  lst_est_after$tau <- lst_est1$tau[keep_rows, , drop = FALSE]
  lst_est_after$pi  <- lst_est1$pi[keep_rows, , drop = FALSE]
  lst_est_after$k   <- lst_est1$k[keep_rows]
  lst_est_after$keep_rows <- keep_rows

  return(lst_est_after)
}