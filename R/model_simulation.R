#' Simulating summary statistics under our simple hierarchical model
#'
#' This model generates independent sets of summary statistics for a target population and K auxiliary populations under a simple version of our hierarchical model, where mu and delta terms are normally distributed with variances tau_mu and tau_delta and all are independent of each other. This uses the mvtnorm package.
#'
#' @param gamma The exposure-outcome causal effect
#' @param tau_mu The variance of the normal distribution underlying mu (central values)
#' @param tau_delta The variance of the normal distribution underlying the delta (deviation terms)
#' @param vars the number of variants to simulate for
#' @param pops the number of populations to simulate data from, consisting of one target population and K-1 auxiliary populations.
#' @param SE_list a list of vectors of desired standard errors for each variant
#' @param r_mat a matrix of estimated (residual) summary statistic correlations due to sample overlap between each set of GWAS summary statistics. If not provided by default, it will create a matrix assuming no sample overlap.
#' @param seed setting a seed for reproducibility
#'
#' @returns a list, consisting of the simulated summary statistics in an N x (K + 1) matrix (beta_matrix), the standard errors of each of the summary statistics (se_matrix), the simulated expectations of the summary statistics (mean_matrix), and the underlying input parameters (r_mat, gamma, tau_mu, and tau_delta) for reproducibility
#' @export
#'
#' @examples
#' SE_list <- rep(list(c(0.4, 0.4, 0.2, 0.1)), 10)
#' simplemodel_sim(gamma = 0, tau_mu = 1, tau_delta = 1, SE_list = SE_list, vars = 10, pops = 3)
#' SE_list <- rep(list(c(0.2, 0.1, 0.05, 0.05)), 50)
#' r_mat <- diag(4)
#' r_mat[1,2] <- 0.2
#' r_mat[2,1] <- 0.2
#' simplemodel_sim(gamma = 0.7, tau_mu = 0.5, tau_delta = 0.2, SE_list = SE_list, vars = 50, pops = 3, r_mat = r_mat)
simplemodel_sim <- function(gamma, tau_mu, tau_delta, vars, pops, SE_list, r_mat = NA, seed = 2025) {
  set.seed(seed)
  # Basic checks to ensure dimensions of standard errors and residual correlations are compatible with variant counts

  #If no r_mat is provided, make a diagonal matrix by default
  if (identical(r_mat, NA)) {
    r_mat <- diag(pops + 1)
  } else { #otherwise, check to make sure that r_mat has the intended input
    if (!is.matrix(r_mat)) stop("r_mat is not a matrix!")
    if (dim(r_mat)[1] != (pops + 1) & dim(r_mat)[2] != (pops + 1)) {
      stop("the dimensions of r_mat are not compatible with the number of populations!")
    }
  }

  if (length(SE_list) != vars) {
    stop("standard error list length is not compatible with number of variants!")
  }
  SE_list_length <- unlist(lapply(SE_list, length))
  SE_unique_k <- unique(SE_list_length)
  if (length(SE_unique_k) > 1) {
    stop("standard error vectors are not the same length!")
  }
  if (SE_unique_k != pops + 1) {
    stop("standard error vectors are not compatible with number of populations!")
  }

  # Simulating expectations for each variant
  mu_vector <- stats::rnorm(n = vars, mean = 0, sd = sqrt(tau_mu))
  delta_matrix <- matrix(stats::rnorm(n = vars*pops, mean = 0, sd = sqrt(tau_delta)),
                         nrow = vars, ncol = pops)
  mean_matrix_exposure <- matrix(rep(mu_vector, pops), nrow = vars, ncol = pops) + delta_matrix
  mean_vector_outcome <- matrix(gamma*(mean_matrix_exposure[,1]), nrow = vars)
  mean_matrix <- cbind(mean_vector_outcome, mean_matrix_exposure)

  # Generating covariance matrices for each variant conditional on mean values
  SE_matrices <- lapply(SE_list, function(X){r_mat * outer(X, X)})

  # Simulating values using row vectors and covariance matrices
  sumstat_matrix <- matrix(nrow = vars, ncol = (pops + 1))
  sumstat_se_matrix <- matrix(nrow = vars, ncol = (pops + 1))

  for (i in 1:vars) {
    sumstat_se_matrix[i,] <- SE_list[[i]]
    sumstat_matrix[i,] <- mvtnorm::rmvnorm(n = 1, mean = mean_matrix[i,],
                                           sigma = SE_matrices[[i]])
  }

  return(list(beta_matrix = sumstat_matrix, se_matrix = sumstat_se_matrix,
         mean_matrix = mean_matrix, r_mat = r_mat,
         gamma = gamma, tau_mu = tau_mu, tau_delta = tau_delta))
}
