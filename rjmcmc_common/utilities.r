library(evd)
library(truncnorm)
library(MCMCpack)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(tictoc)
library(kde1d)
library(gtools)
library(cowplot)

max.K <- 10
threshold <- 1e-300
tau_threshold <- 0.01
left_boundary_TN <- 0.001
birth_pr <- c(0.4, rep(0.2, (max.K - 2)), 0)
death_pr <- c(0, rep(0.2, (max.K - 2)), 0.4)

is_debug <- FALSE

tcat <- function(...) {
  if (is_debug) {
    cat(...)
  }
}

generate_lst <- function(...) {
  args <- list(...)
  names(args) <- as.character(match.call()[-1])
  return(args)
}

extract_data <- function(est, start_point, end_point) {
  lst <- list()
  lst$k <- rep(0, end_point - start_point + 1)
  lst$mu <- matrix(0, end_point - start_point + 1, max.K)
  lst$tau <- matrix(0, end_point - start_point + 1, max.K)
  lst$pi <- matrix(0, end_point - start_point + 1, max.K)
  for (t in c(1:(end_point - start_point + 1))) {
    k <- est[t + start_point - 1, 1]
    lst$k[t] <- k
    lst$mu[t, 1:k] <- est[t + start_point - 1, (2 + 0 * k):(1 * k + 1)]
    lst$tau[t, 1:k] <- est[t + start_point - 1, (2 + 1 * k):(2 * k + 1)]
    lst$pi[t, 1:k] <- est[t + start_point - 1, (2 + 2 * k):(3 * k + 1)]
  }
  return(lst)
}


remove_all_k_but <- function(lst, keep_k) {
  index_to_remove <- which(lst$k != keep_k)
  lst$k <- lst$k[-index_to_remove]
  lst$mu <- lst$mu[-index_to_remove, , drop = FALSE]
  lst$tau <- lst$tau[-index_to_remove, , drop = FALSE]
  lst$pi <- lst$pi[-index_to_remove, , drop = FALSE]
  return(lst)
}
