#' The marginal log-likelihood for the summary statistics of a single variant under our meta-analysis model, assuming no selection, and i.i.d. normal mu and delta
#'
#' @param gamma The exposure-outcome causal effect parameter
#' @param tau_mu The variance of the normal distribution underlying mu
#' @param tau_delta The variance of the normal distributions underlying delta
#' @param r_mat a matrix of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide by default.
#' @param is_overlap whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param sumstat_beta the GWAS effect size estimates for a single variant across K+1 studies: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations.
#' @param sumstat_se the standard errors of the effect size estimates in sumstat_beta
#'
#' @returns A log-likelihood
#' @export
#'
#' @examples
#' simple_loglik_single(sumstat_beta = c(0,0,0), sumstat_se = c(1,1,1), gamma = 0, tau_mu = 1, tau_delta = 1)
#' simple_loglik_single(sumstat_beta = c(0.4, 0.6, 0.3), sumstat_se = c(0.1, 0.1, 0.05), is_overlap = TRUE, r_mat = matrix(c(1, 0.2, 0, 0.2, 1, 0, 0, 0, 1), nrow = 3, ncol = 3), gamma = 0.8, tau_mu = 1, tau_delta = 1)
simple_loglik_single <- function(sumstat_beta, sumstat_se, gamma, is_overlap = FALSE, r_mat = NA, tau_mu, tau_delta) {
  # Initial checks
  if (length(sumstat_beta) != length(sumstat_se)) {
    stop("Summary statistic effect size estimates and standard errors are not the same length!")
  }
  if (is_overlap == TRUE & identical(r_mat, NA)) {
    stop("If you assume sample overlap, please provide a residual correlation matrix r_mat!")
  }
  if (!identical(r_mat, NA)) {
    if (dim(r_mat)[1] != dim(r_mat)[2]) stop("Please provide a square residual correlation matrix r_mat!")
    if (length(sumstat_beta) != dim(r_mat)[1]) stop("r_mat and summary statistic vectors are incompatible!")
  }

  #Setting up the covariance matrix
  K <- length(sumstat_beta) - 1

  if (is_overlap == FALSE & identical(r_mat, NA)) { #the default is a standard diagonal matrix
    r_mat <- diag(K+1)
  }
  sigma_Ej <- r_mat * outer(sumstat_se, sumstat_se)

  sigma_G <- matrix(0, nrow = K + 1, ncol = K + 1)
  sigma_G[,1] <- c(gamma^2 * (tau_mu + tau_delta), gamma * (tau_mu + tau_delta), rep(gamma*tau_mu, K - 1))
  sigma_G[1,] <- c(gamma^2 * (tau_mu + tau_delta), gamma * (tau_mu + tau_delta), rep(gamma*tau_mu, K - 1))
  sigma_G[2:(K+1), 2:(K+1)] <- diag(tau_delta, K) + matrix(tau_mu, nrow = K, ncol = K)

  sigma_j <- sigma_Ej + sigma_G

  #Calculating the log-likelihood
  logLik <- -1/2 * (log(det(sigma_j)) + t(sumstat_beta) %*% solve(sigma_j) %*% t(t(sumstat_beta)) + (K+1)*log(2*pi))
  return(logLik)
}
