# Mixture Analysis of Gumbel (RJMCMC + Particle Filter)

Code for my Master of Mathematical Science thesis at the University of Canterbury (supervised by Dr Elena Moltchanova, Dr Philipp Wacker, and Dr Christoph Bartneck) (Link will be updated later). The work focuses on mixture analysis of the Gumbel distribution, combining Reversible Jump MCMC (RJMCMC) for finite mixture models (FMM) and particle filters for Dirichlet Process Mixture Models (DPMM). Normal and exponential variants are included for comparison.

## Repository layout
- `.r` files: main function code; `.Rmd` files: runnable examples and experiment scripts (e.g., `case_study.rmd`, `time.Rmd`, `rjmcmc_gumbel/rjmcmc_gumbel.Rmd`, `particle_filter/gumbel/gumbel_pf.Rmd`).
- `rjmcmc_common/`: shared RJMCMC utilities such as proposal moves and plotting helpers.
- `rjmcmc_gumbel/`: RJMCMC for Gumbel mixtures, including acceptance logic, split/combine moves, and notebooks to reproduce experiments.
- `rjmcmc_normal/`: RJMCMC for normal mixtures plus DIC calculations.
- `particle_filter/`: sequential Monte Carlo implementations for DPMM with Gumbel, exponential, and normal variants (`gumbel/`, `exp/`, `normal_fearnhead/`, `normal_not_fearnhead/`); shared helpers live in `pf_utilities.r`, `fearnhead_utilities.r`, and `fc_resample.r`.

## Prerequisites
- R with packages used across the scripts: `evd`, `truncnorm`, `MCMCpack`, `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `tictoc`, `mixAK`, `kde1d`, `gtools`, `cowplot`, `here`, `vistime`, and related dependencies. `dir_manager.r` contains a small install helper.
- Data for the case study is not included; provide your own dataset and update paths in `case_study.rmd` (currently assumes a swimming performance CSV).


## Notes
- Code favors clarity over efficiency and has not been optimized.
- If your paths differ from the defaults in `dir_manager.r`, adjust `dir_path` and `pic_path` there before running.
