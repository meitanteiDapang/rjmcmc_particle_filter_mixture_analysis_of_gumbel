# Get the current script's path, even when sourced

# packages.R
packages <- c(
  "kde1d",
  "truncnorm",
  "MCMCpack",
  "vistime",
  "ggplot2",
  "here",
  "evd",
  "patchwork",
  "tictoc",
  "gtools"
)
real_pi <- pi

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# invisible(sapply(packages, install_if_missing))


library(here)



dir_path <- here::here()
pic_path <- "/Users/yangzhonghao/Desktop/master_pic/"
