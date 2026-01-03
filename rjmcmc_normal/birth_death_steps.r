source(file.path(dir_path, "rjmcmc_normal/acceptance_functions.r"))

birth_death <- function(lst) {
  k <- lst$k
  birth_pr <- lst$birth_pr
  death_pr <- lst$death_pr
  iii <- lst$i
  # propose new k

  move <-
    sample(c(1, -1, 0), 1,
      prob =
        c(birth_pr[k], death_pr[k], 1 - birth_pr[k] - death_pr[k])
    )
  if (move == 0) {
    return(lst)
  }

  new_k <- lst$k + move
  is_birth <- ifelse(move == 1, TRUE, FALSE)

  # extract parameters
  pi <- lst$pi
  mu <- lst$mu
  tau <- lst$tau
  k <- lst$k
  Z <- lst$Z
  X <- lst$X
  n <- lst$n
  Xv <- lst$Xv
  nv <- lst$nv
  xi <- lst$xi
  kappa <- lst$kappa
  a <- lst$a
  b <- lst$b
  lambda <- lst$lambda
  delta <- lst$delta

  # pre-assign
  new_mu <- rep(0, new_k)
  new_tau <- rep(0, new_k)
  new_pi <- rep(0, new_k)
  new_Z <- lst$Z
  new_Xv <- vector("list", new_k)
  new_nv <- rep(0, new_k)

  # empty label
  empty_labels <- which(nv == 0)
  k0 <- length(empty_labels)

  if (is_birth) {
    # birth case
    if (exists("birth_count")) {
      birth_count <<- birth_count + 1
    }

    # draw
    pi_jstar <- rbeta(1, 1, k)
    mu_jstar <- rnorm(1, mean = xi, sd = 1 / sqrt(kappa))
    tau_jstar <- max(rgamma(1, a, b), tau_threshold)

    # find the insert pos for mu_jstar
    insert_index <- which(mu > mu_jstar)[1]
    # if mu_jstar is the biggest, put it in the last
    if (is.na(insert_index)) {
      insert_index <- length(mu) + 1
    }

    # increase k0 by 1 by definition
    k0 <- k0 + 1

    # insert
    new_mu <- append(mu, mu_jstar, after = insert_index - 1)
    # for pi, the weight of other position will mutiply 1 - pi_jstar
    new_pi <- append((1 - pi_jstar) * pi, pi_jstar, after = insert_index - 1)
    new_tau <- append(tau, tau_jstar, after = insert_index - 1)
    new_nv <- append(nv, 0, after = insert_index - 1)

    # birth change Z
    new_Z <- sapply(Z, function(z) {
      if (z >= insert_index) {
        return(z + 1)
      } else {
        return(z)
      }
    })
  } else {
    # death case
    if (exists("death_count")) {
      death_count <<- death_count + 1
    }
    # nobody to die
    if (k0 == 0) {
      return(lst)
    }
    
    # choose one to die
    die_index <- if (length(empty_labels) == 1) {
      empty_labels
    } else {
      sample(empty_labels, size = 1)
    }

    # get data
    pi_jstar <- pi[die_index]

    # remove
    new_pi <- pi[-die_index] / (1 - pi_jstar)
    new_mu <- mu[-die_index]
    new_tau <- tau[-die_index]
    new_nv <- nv[-die_index]
    Xv[die_index] <- NULL


    # dead change Z
    new_Z <- sapply(Z, function(z) {
      if (z > die_index) {
        return(z - 1)
      } else if (z == die_index) {
        cat("Cannot happen, Z should not be equal to die_index\n")
        return(NA) # This should not happen, but just in case
      } else {
        return(z)
      }
    })


  }

  acceptance_r <- birth_death_acceptance(
    is_birth, ifelse(is_birth, lst$k, new_k), k0,
    delta, lambda,
    pi_jstar, n,
    birth_pr, death_pr
  )
  if (runif(1, 0, 1) < acceptance_r) {
    # accepted
    lst$Z <- new_Z
    lst$k <- new_k
    lst$pi <- new_pi
    lst$mu <- new_mu
    lst$tau <- new_tau
    lst$Xv <- new_Xv
    lst$nv <- new_nv
    if (exists("birth_accept_count")) {
      if (is_birth) {
        # cat("birth 123123acceptanmce: ", acceptance_r, "\n")
        birth_accept_count <<- birth_accept_count + 1
      } else {
        # cat("death acceptanmce: ", acceptance_r, "\n")
        death_accept_count <<- death_accept_count + 1
      }
    }
  }
  return(lst)
}
