
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MetaMR

<!-- badges: start -->

<!-- badges: end -->

MetaMR is an R package for multi-ancestry Mendelian randomization. It
estimates the causal effect in a target population using exposure and
outcome GWAS summary statistics from that target population, as well as
exposure summary statistics from $K-1$ other auxiliary populations.

## Installation

MetaMR uses the `mvtnorm` package as a dependency.

    install.packages("mvtnorm")

To install `MetaMR` from GitHub:

    devtools::install_github("lijackk/MetaMR")
