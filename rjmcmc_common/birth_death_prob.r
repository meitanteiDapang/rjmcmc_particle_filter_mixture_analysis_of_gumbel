get_birth_death_pr <- function(max.K, min.K = 1, lambda, c_val = 0.5) {
  # Initialize vectors
  birth_pr <- numeric(max.K)
  death_pr <- numeric(max.K)

  for (k in min.K:max.K) {
    if (k < max.K) {
      birth_pr[k] <- c_val * min(1, lambda / (k + 1))
    } else {
      birth_pr[k] <- 0
    }

    if (k > min.K) {
      death_pr[k] <- c_val * min(1, k / lambda)
    } else {
      death_pr[k] <- 0
    }
  }

  list(
    birth_pr = birth_pr,
    death_pr = death_pr
  )
}
