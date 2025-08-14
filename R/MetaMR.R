#' MetaMR using a simple hierarchical model (NOT YET COMPLETE)
#'
#' This function estimates an exposure-outcome causal effect, \eqn{\gamma}, in the target population, using exposure and outcome summary statistics in the target population as well as exposure summary statistics from \eqn{K - 1} auxiliary populations using a simple hierarchical model. This function inputs a list of GWAS effect size estimates across all \eqn{K+1} studies for a set of \eqn{n} variants selected as instrumental variables, a list of standard errors for these variants, and a user-specified residual correlation matrix reflecting sample overlap between studies, then maximizes a likelihood that models exposure summary statistics for each variant \eqn{j} under a random effects meta-analysis framework based on central values \eqn{\mu_j} and population-specific deviations \eqn{\delta_{kj}}, where \eqn{\mu_j} are i.i.d. normal with variance \eqn{\tau_{\mu}} and \eqn{\delta_{kj}} are also i.i.d. normal with variance \eqn{\tau_{\delta}}. Note that this function does NOT yet incorporate variant selection, potential horizontal pleiotropy in the outcome, or alternative distributions for \eqn{\mu} and \eqn{\delta} values.
#' @param sumstat_beta_list a list of vectors of GWAS effect size estimates, each with length K+1: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations. Each vector in the list represents summary statistics for one of N variants.
#' @param sumstat_se_list a list of vectors of the standard errors for GWAS effect size estimates in sumstat_beta_list
#' @param is_overlap  boolean describing whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat_list a list of matrices of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap.
#'
#' @returns a list that includes point estimates for the exposure-outcome causal effect \eqn{\gamma} and variance terms \eqn{\tau_{\mu}} and \eqn{\tau_{\delta}} estimated using maximum likelihood, the covariance matrix of these point estimates (calculated using the inverse of the Hessian), and the results of a likelihood ratio test that tests against null hypothesis \eqn{\gamma = 0} (including the test statistic and p-value)
#' @export
#'
#' @examples
#' SE_list <- rep(list(c(0.2, 0.1, 0.05, 0.05)), 50)
#' r_mat <- diag(4)
#' r_mat[1,2] <- 0.2
#' r_mat[2,1] <- 0.2
#' observed_data <- simplemodel_sim(gamma = 0.7, tau_mu = 0.05, tau_delta = 0.02, SE_list = SE_list, vars = 50, pops = 3, r_mat = r_mat)
#' sumstat_beta_list <- apply(observed_data$beta_matrix, MARGIN = 1, function(x) {return(x)}, simplify = FALSE)
#' MetaMR_simplemodel(sumstat_beta_list = sumstat_beta_list, sumstat_se_list = SE_list, r_mat_list = rep(list(r_mat), 50))
#'
MetaMR_simplemodel <- function(sumstat_beta_list, sumstat_se_list, is_overlap = FALSE, r_mat_list= NA) {
  #The usual battery of checks before we proceed
  if (length(sumstat_beta_list) != length(sumstat_se_list)) {
    stop("Summary statistic effect size estimates and standard errors imply differing numbers of variants!")
  }
  k_counts_beta <- unlist(lapply(sumstat_beta_list, length))
  k_counts_se <- unlist(lapply(sumstat_se_list, length))
  unique_k_counts <- unique(c(k_counts_beta, k_counts_se))
  if (length(unique_k_counts) > 1) {
    stop("Different variants have summary statistic vectors of different lengths!")
  }

  if (is_overlap == TRUE & identical(r_mat_list, NA)) {
    stop("If you assume sample overlap, please provide residual correlation matrices in r_mat_list!")
  }

  n <- length(sumstat_beta_list) #number of variants
  k <- unique_k_counts - 1 #number of populations

  #If no residual correlation matrix is given due to a lack of overlap, generate diag(k) matrices by default
  if (is_overlap == FALSE & identical(r_mat_list, NA)) { #the default is a standard diagonal matrix for all matrices
    diag_mat <- list(diag(k + 1))
    r_mat_list <- rep(diag_mat, n)
  }

  #Optimizing the full likelihood
  MetaMR_point_est <- simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                             sumstat_se_list = sumstat_se_list,
                                             is_overlap = is_overlap,
                                             r_mat_list = r_mat_list)

  #Optimizing the constrained likelihood under the null hypothesis
  MetaMR_null_loglik <- simple_loglik_optimize_null(sumstat_beta_list = sumstat_beta_list,
                                               sumstat_se_list = sumstat_se_list,
                                               is_overlap = is_overlap,
                                               r_mat_list = r_mat_list)

  #Obtaining the variance of each of the estimated parameters using the information matrix
    #Information matrix = negative hessian of log-likelihood at MLE
  covar_matrix <- solve(MetaMR_point_est$hessian)

  #Performing a likelihood ratio test
  LRT_test_stat <- 2 * (MetaMR_null_loglik$value - MetaMR_point_est$value)
  LRT_pval <- stats::pchisq(LRT_test_stat, df = 1, lower.tail = FALSE)

  return(list(optim_object = MetaMR_point_est,
              covar = covar_matrix,
              LRT = list(test_stat = LRT_test_stat, pval = LRT_pval)))
}
