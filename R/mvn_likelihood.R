#' Log-likelihood function under simple model (single-variant)
#'
#' The marginal log-likelihood for the summary statistics of a single variant under our meta-analysis model, assuming no selection, and i.i.d. normal mu and delta
#'
#' @param gamma The exposure-outcome causal effect parameter
#' @param tau_mu The variance of the normal distribution underlying mu. It should be log-transformed unless set to be exactly zero.
#' @param tau_delta The variance of the normal distributions underlying delta. It should be log-transformed unless set to be exactly zero.
#' @param r_mat a matrix of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap.
#' @param is_overlap a boolean whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param sumstat_beta the GWAS effect size estimates for a single variant across K+1 studies: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations.
#' @param sumstat_se the standard errors of the effect size estimates in sumstat_beta
#' @param tau_mu_log Whether tau_mu represents its log-transformed value or not
#' @param tau_delta_log Whether tau_delta represents its log-transformed value or not
#'
#' @returns A log-likelihood for the summary statistics of a single variant given some parameters gamma, mu, tau
#' @export
#'
#' @examples
#' simple_loglik_single(sumstat_beta = c(0,0,0), sumstat_se = c(1,1,1),
#'                      gamma = 0, tau_mu = 1, tau_delta = 1)
#' simple_loglik_single(sumstat_beta = c(0,0,0), sumstat_se = c(1,1,1),
#'                      gamma = 0, tau_mu = 1, tau_delta = 0)
#' simple_loglik_single(sumstat_beta = c(0.4, 0.6, 0.3), sumstat_se = c(0.1, 0.1, 0.05),
#'                      is_overlap = TRUE,
#'                      r_mat = matrix(c(1, 0.2, 0, 0.2, 1, 0, 0, 0, 1), nrow = 3, ncol = 3),
#'                      gamma = 0.8, tau_mu = 1, tau_delta = 1)
#' simple_loglik_single(sumstat_beta = c(0.4, 0.6, 0.3), sumstat_se = c(0.1, 0.1, 0.05),
#'                      is_overlap = TRUE,
#'                      r_mat = matrix(c(1, 0.2, 0, 0.2, 1, 0, 0, 0, 1), nrow = 3, ncol = 3),
#'                      gamma = 0.8, tau_mu = 0, tau_delta = 0, tau_mu_log = TRUE, tau_delta_log = TRUE)
simple_loglik_single <- function(sumstat_beta, sumstat_se, gamma, is_overlap = FALSE, r_mat = NA, tau_mu, tau_delta, tau_mu_log = FALSE, tau_delta_log = FALSE) {
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

  if (tau_mu_log) { #if tau_mu is log-transformed, exponentiate
    tau_mu <- exp(tau_mu)
  }

  if (tau_delta_log) { #if tau_delta is log-transformed, exponentiate
    tau_delta <- exp(tau_delta)
  }

  sigma_G <- matrix(0, nrow = K + 1, ncol = K + 1)
  sigma_G[,1] <- c(gamma^2 * (tau_mu + tau_delta), gamma * (tau_mu + tau_delta), rep(gamma*tau_mu, K - 1))
  sigma_G[1,] <- c(gamma^2 * (tau_mu + tau_delta), gamma * (tau_mu + tau_delta), rep(gamma*tau_mu, K - 1))
  sigma_G[2:(K+1), 2:(K+1)] <- diag(tau_delta, K) + matrix(tau_mu, nrow = K, ncol = K)

  sigma_j <- sigma_Ej + sigma_G

  #Calculating the log-likelihood
  logLik <- -1/2 * (log(det(sigma_j)) + t(sumstat_beta) %*% solve(sigma_j) %*% t(t(sumstat_beta)) + (K+1)*log(2*pi))
  return(logLik)
}

#' Log-likelihood function under simple model (multi-variant)
#'
#' The marginal log-likelihood for the summary statistics of a list of variants under our meta-analysis model, assuming no selection, and i.i.d. normal mu and delta
#'
#' @param sumstat_beta_list a list of vectors of GWAS effect size estimates with length K+1: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations. Each vector in the list represents summary statistics for one of many variants.
#' @param sumstat_se_list a list of vectors of the standard errors for GWAS effect size estimates in sumstat_beta_list
#' @param gamma The exposure-outcome causal effect parameter
#' @param is_overlap a boolean describing whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat_list a list of matrices of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap
#' @param tau_mu The variance of the normal distribution underlying mu. It should be log-transformed unless set to be exactly zero.
#' @param tau_delta The variance of the normal distributions underlying delta. It should be log-transformed unless set to be exactly zero.
#' @param tau_mu_log Whether tau_mu represents its log-transformed value or not.
#' @param tau_delta_log Whether tau_delta represents its log-transformed value or not.
#'
#' @returns A log-likelihood for the summary statistics of a set of variants given some parameters gamma, mu, tau
#' @export
#'
#' @examples
#' simple_loglik_full(sumstat_beta_list = list(c(0,0,0), c(0,0,0), c(0,0,0)),
#'                      sumstat_se_list = list(c(1,1,1), c(1,1,1), c(1,1,1)),
#'                      gamma = 0, tau_mu = 1, tau_delta = 1)
#' simple_loglik_full(sumstat_beta_list = list(c(0,0,0), c(0,0,0), c(0,0,0)),
#'                      sumstat_se_list = list(c(1,1,1), c(1,1,1), c(1,1,1)),
#'                      gamma = 0, tau_mu = 1, tau_delta = 0)
#' simple_loglik_full(sumstat_beta_list = list(c(0.4, 0.6, 0.3), c(0.7,1,0.2), c(0.3,0.36,0.3)),
#'                      sumstat_se_list = list(c(0.1, 0.1, 0.05), c(0.2, 0.2, 0.1), c(0.1, 0.1, 0.05)),
#'                      is_overlap = TRUE,
#'                      r_mat_list <- rep(list(matrix(c(1, 0.2, 0, 0.2, 1, 0, 0, 0, 1), nrow = 3, ncol = 3)), 3),
#'                      gamma = 0.8, tau_mu = 1, tau_delta = 1)
#' simple_loglik_full(sumstat_beta_list = list(c(0.4, 0.6, 0.3), c(0.7,1,0.2), c(0.3,0.36,0.3)),
#'                      sumstat_se_list = list(c(0.1, 0.1, 0.05), c(0.2, 0.2, 0.1), c(0.1, 0.1, 0.05)),
#'                      is_overlap = TRUE,
#'                      r_mat_list <- rep(list(matrix(c(1, 0.2, 0, 0.2, 1, 0, 0, 0, 1), nrow = 3, ncol = 3)), 3),
#'                      gamma = 0.8, tau_mu = 0, tau_delta = 0, tau_mu_log = TRUE, tau_delta_log = TRUE)
simple_loglik_full <- function(sumstat_beta_list, sumstat_se_list, gamma, is_overlap = FALSE, r_mat_list = NA, tau_mu, tau_delta, tau_mu_log = FALSE, tau_delta_log = FALSE) {
  #Initial checks - will make more later if necessary
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

  n <- length(sumstat_beta_list)
  k <- unique_k_counts - 1

  #If no residual correlation matrix is given due to a lack of overlap, generate diag(k) matrices by default
  if (is_overlap == FALSE & identical(r_mat_list, NA)) { #the default is a standard diagonal matrix for all matrices
    diag_mat <- list(diag(k + 1))
    r_mat_list <- rep(diag_mat, n)
  }

  loglik_variants <- vector()
  for (i in 1:n) {
    loglik_variants[i] <- simple_loglik_single(sumstat_beta = sumstat_beta_list[[i]],
                                               sumstat_se = sumstat_se_list[[i]],
                                               is_overlap = TRUE,
                                               r_mat = r_mat_list[[i]],
                                               gamma = gamma,
                                               tau_mu = tau_mu,
                                               tau_delta = tau_delta,
                                               tau_mu_log = tau_mu_log,
                                               tau_delta_log = tau_delta_log)
  }
  return(sum(loglik_variants))
}

#' Optimize log-likelihood of simple model for given summary statistics
#'
#' For a list of summary statistics (effect size estimates and standard errors) for n variants across K + 1 studies (outcome from target population, exposure from target population, exposures from K-1 auxiliary populations), we use the optim() function to minimize the negative log-likelihood under the simple hierarchical model of MetaMR over the exposure-outcome causal effect (gamma), as well as the variances of mu and delta (tau_mu and tau_delta). We also allow for subsets of these parameters to be held constant.
#'
#' @param sumstat_beta_list a list of vectors of GWAS effect size estimates, each with length K+1: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations. Each vector in the list represents summary statistics for one of N variants.
#' @param sumstat_se_list a list of vectors of the standard errors for GWAS effect size estimates in sumstat_beta_list
#' @param is_overlap a boolean describing whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat_list a list of matrices of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap.
#' @param is.fixed which parameters of interested should be fixed during optimization. Follows order (gamma, tau_mu, tau_delta).
#' @param fix.params if a parameter is set as fixed in is.fixed, what it should be. Follows order (gamma, tau_mu, tau_delta). Our function will check to make sure is.fixed and fix.params are compatible with each other.
#' @param tau_mu_log whether tau_mu is log-transformed. It should be log-transformed unless set to be exactly zero.
#' @param tau_delta_log whether tau_delta is log-transformed. It should be log-transformed unless set to be exactly zero.
#'
#' @returns the output from optim(), detailing the optimized values for gamma, tau_mu and tau_delta, information about convergence, and the negative log-likelihood
#' @export
#'
#' @examples
#' SE_list <- rep(list(c(0.2, 0.1, 0.05, 0.05)), 50)
#' r_mat <- diag(4)
#' r_mat[1,2] <- 0.2
#' r_mat[2,1] <- 0.2
#' observed_data <- simplemodel_sim(gamma = 0.7, tau_mu = 0.5, tau_delta = 0.2, SE_list = SE_list, vars = 50, pops = 3, r_mat = r_mat)
#' sumstat_beta_list <- apply(observed_data$beta_matrix, MARGIN = 1, function(x) {return(x)}, simplify = FALSE)
#' simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list, sumstat_se_list = SE_list, r_mat_list = rep(list(r_mat), 50), tau_mu_log = TRUE, tau_delta_log = TRUE)
#'
#' observed_data <- simplemodel_sim(gamma = 0.7, tau_mu = 0.5, tau_delta = 0, SE_list = SE_list, vars = 50, pops = 3, r_mat = r_mat)
#' sumstat_beta_list <- apply(observed_data$beta_matrix, MARGIN = 1, function(x) {return(x)}, simplify = FALSE)
#' simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list, sumstat_se_list = SE_list, r_mat_list = rep(list(r_mat), 50), is.fixed = c(FALSE, FALSE, TRUE), fix.params = c(NA, NA, 0), tau_mu_log = TRUE, tau_delta_log = FALSE)
#'
simple_loglik_optimize <- function(sumstat_beta_list, sumstat_se_list, is_overlap = FALSE, r_mat_list= NA, is.fixed = c(FALSE, FALSE, FALSE), fix.params = c(NA, NA, NA), tau_mu_log = FALSE, tau_delta_log = FALSE) {

  #Whether is.fixed and fix.params are compatible
    #If is.fixed[i] is FALSE, is.na(fix.params)[i] should be TRUE
    #If is.fixed[i] is TRUE, is.na(fix.params)[i] should be FALSE.
  fixed_check <- xor(is.fixed, is.na(fix.params))
  if (any(fixed_check == FALSE)) {
    stop("is.fixed and fix.params are not compatible!")
  }

  #first perform some basic checks, the same as from simple_loglik_full
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

  #Setting up the function and parameters to optimize over
  init.params <- c(gamma = 0, tau_mu = 0, tau_delta = 0)[!is.fixed]
  optim_fn <- function(params) {
    gamma <- ifelse("gamma" %in% names(init.params), params["gamma"], fix.params[1])
    tau_mu <- ifelse("tau_mu" %in% names(init.params), params["tau_mu"], fix.params[2])
    tau_delta <- ifelse("tau_delta" %in% names(init.params), params["tau_delta"], fix.params[3])

    return(-simple_loglik_full(sumstat_beta_list = sumstat_beta_list,
                               sumstat_se_list = sumstat_se_list,
                               gamma = gamma,
                               is_overlap = is_overlap,
                               r_mat_list = r_mat_list,
                               tau_mu = tau_mu,
                               tau_delta = tau_delta,
                               tau_mu_log = tau_mu_log,
                               tau_delta_log = tau_delta_log))
  }

  #minimize the negative log-likelihood
  stats::optim(par = init.params,
        fn = optim_fn,
        hessian = TRUE)
}


#' Optimize log-likelihood of simple model for given summary statistics under the null hypothesis (DEPRECATED)
#'
#' Similar to simple_loglik_optimize, but now \eqn{\gamma} is constrained to be exactly zero. Used for likelihood ratio testing. NO LONGER USED DUE TO A MORE FLEXIBLE BASE OPTIMIZATION FUNCTION.
#' @param sumstat_beta_list a list of vectors of GWAS effect size estimates, each with length K+1: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations. Each vector in the list represents summary statistics for one of N variants.
#' @param sumstat_se_list a list of vectors of the standard errors for GWAS effect size estimates in sumstat_beta_list
#' @param is_overlap a boolean describing whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat_list a list of matrices of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap.
#'
#' @returns the output from optim(), detailing the optimized values for tau_mu and tau_delta, information about convergence, and the negative log-likelihood
#' @export
#'
#' @examples
#' SE_list <- rep(list(c(0.2, 0.1, 0.05, 0.05)), 50)
#' r_mat <- diag(4)
#' r_mat[1,2] <- 0.2
#' r_mat[2,1] <- 0.2
#' observed_data <- simplemodel_sim(gamma = 0.7, tau_mu = 0.5, tau_delta = 0.2, SE_list = SE_list, vars = 50, pops = 3, r_mat = r_mat)
#' sumstat_beta_list <- apply(observed_data$beta_matrix, MARGIN = 1, function(x) {return(x)}, simplify = FALSE)
#' simple_loglik_optimize_null(sumstat_beta_list = sumstat_beta_list, sumstat_se_list = SE_list, r_mat_list = rep(list(r_mat), 50))
#'
simple_loglik_optimize_null <- function(sumstat_beta_list, sumstat_se_list, is_overlap = FALSE, r_mat_list= NA) {
  #first perform some basic checks, the same as from simple_loglik_full
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

  #minimize the negative log-likelihood over gamma, tau_mu, and tau_delta
  stats::optim(par = c(tau_mu = 1, tau_delta = 1),
               fn = function(params) {
                 -simple_loglik_full(sumstat_beta_list = sumstat_beta_list,
                                     sumstat_se_list = sumstat_se_list,
                                     gamma = 0,
                                     is_overlap = is_overlap,
                                     r_mat_list = r_mat_list,
                                     tau_mu = params[1],
                                     tau_delta = params[2])
               },
               method = "L-BFGS-B", lower = c(1e-16, 1e-16), upper = c(Inf, Inf),
               hessian = TRUE)
}
