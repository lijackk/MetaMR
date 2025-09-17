#' Log-likelihood function under simple model (single-variant)
#'
#' The marginal log-likelihood for the summary statistics of a single variant under our meta-analysis model, assuming i.i.d. normal mu and delta. Now includes an option to specify whether variants were selected based on a Z-score threshold for the target exposure.
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
#' @param select_zscore if the instrument was selected because its association with the exposure in the target exceeded some Z-score threshold, what threshold this is. If no selection was performed to obtain the summary statistics, put NA.
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
#'                      gamma = 0.8, tau_mu = 0, tau_delta = 0,
#'                      tau_mu_log = TRUE, tau_delta_log = TRUE)
simple_loglik_single <- function(sumstat_beta, sumstat_se, gamma, is_overlap = FALSE, r_mat = NA, tau_mu, tau_delta, tau_mu_log = FALSE, tau_delta_log = FALSE, select_zscore = NA) {
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

  #Bringing in the possibility of selection
  #log-likelihood = joint/marginal; marginal = probability of selection, joint = same as original log-likelihood for selected variants
  if (!is.na(as.numeric(select_zscore))) {
    logptau_thresh <- log(2 * stats::pnorm(abs(select_zscore) * sumstat_se[2]/sqrt(tau_mu + tau_delta + sumstat_se[2]^2), lower.tail = FALSE))
    logLik <- logLik - logptau_thresh
  }

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
#' @param select_zscore if the instrument was selected because its association with the exposure in the target exceeded some Z-score threshold, what threshold this is. If no selection was performed to obtain the summary statistics, put NA.
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
#'                      r_mat_list <- rep(list(matrix(c(1, 0.2, 0, 0.2, 1, 0, 0, 0, 1),
#'                                        nrow = 3, ncol = 3)), 3),
#'                      gamma = 0.8, tau_mu = 1, tau_delta = 1)
#' simple_loglik_full(sumstat_beta_list = list(c(0.4, 0.6, 0.3), c(0.7,1,0.2), c(0.3,0.36,0.3)),
#'                      sumstat_se_list = list(c(0.1, 0.1, 0.05), c(0.2, 0.2, 0.1), c(0.1, 0.1, 0.05)),
#'                      is_overlap = TRUE,
#'                      r_mat_list <- rep(list(matrix(c(1, 0.2, 0, 0.2, 1, 0, 0, 0, 1),
#'                                        nrow = 3, ncol = 3)), 3),
#'                      gamma = 0.8, tau_mu = 0, tau_delta = 0,
#'                      tau_mu_log = TRUE, tau_delta_log = TRUE)
simple_loglik_full <- function(sumstat_beta_list, sumstat_se_list, gamma, is_overlap = FALSE, r_mat_list = NA, tau_mu, tau_delta, tau_mu_log = FALSE, tau_delta_log = FALSE, select_zscore = NA) {
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

  #If selection was performed, check whether inputted z-scores actually fulfill selection threshold
  if (!is.na(as.numeric(select_zscore))) {
    zscores <- unlist(lapply(sumstat_beta_list, function(X){X[2]})) / unlist(lapply(sumstat_se_list, function(X){X[2]}))
    if (any(abs(zscores) < select_zscore)) {
      stop("Not all selected variants fulfill the specified target exposure association z-score selection criteria!")
    }
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
                                               tau_delta_log = tau_delta_log,
                                               select_zscore = select_zscore)
  }
  return(sum(loglik_variants))
}

#' Setting Initial Parameters for Better Optim()
#'
#' Numerical optimization methods may be sensitive to choice of initial parameters. This function uses a method of moments-like estimator to set more plausible initial parameters for optimization of the log-likelihood function over gamma, tau_mu and delta.
#' Note that this initial parametrization does NOT factor in variant selection. I have not yet tested the reliability of this initialization in this scenario.
#'
#' @param sumstat_beta_list a list of vectors of GWAS effect size estimates with length K+1: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations. Each vector in the list represents summary statistics for one of many variants.
#' @param tau_mu_log whether tau_mu is log-transformed. It should be log-transformed unless set to be exactly zero.
#' @param tau_delta_log whether tau_delta is log-transformed. It should be log-transformed unless set to be exactly zero.
#' @param use_tau_delta Whether to only estimate tau_mu or tau_delta when no auxiliary populations are provided.
#' @param sumstat_se_list a list of vectors of the standard errors for GWAS effect size estimates in sumstat_beta_list
#'
#' @returns A named vector of initial starting parameters for optim()
#' @export
#'
#' @examples
#' nvar <- 100
#' gamma <- 0
#' h2 <- 0.20
#' persnp_h2 <- h2 / nvar
#' tau_mu <- 0 #persnp_h2 to set tau_delta = 0, 0 to set tau_delta = 1
#' tau_delta <- persnp_h2 - tau_mu
#' N_k <- c(100000, 10000, 40000, 50000)
#' sigma_Xk <- 1/sqrt(N_k)
#' SE_list <- rep(list(sigma_Xk), nvar)
#'
#' i <- 1
#' observed_data <- simplemodel_sim(gamma = gamma, tau_mu = tau_mu, tau_delta = tau_delta,
#'                                  SE_list = SE_list, vars = nvar, pops = 3, seed = i)
#' sumstat_beta_list <- apply(observed_data$beta_matrix, MARGIN = 1, function(x) {return(x)},
#'                            simplify = FALSE)
#'
#' sumstat_beta_list0 <- lapply(sumstat_beta_list, function(x){return(x[1:2])})
#' sumstat_beta_list2 <- lapply(sumstat_beta_list, function(x){return(x[1:4])})
#' SE_list0 <- lapply(SE_list, function(x){return(x[1:2])})
#' SE_list2 <- lapply(SE_list, function(x){return(x[1:4])})
#' set_initial_params(sumstat_beta_list = sumstat_beta_list0, sumstat_se_list = SE_list0,
#'                    tau_mu_log = TRUE, tau_delta_log = TRUE)
#' set_initial_params(sumstat_beta_list = sumstat_beta_list2, sumstat_se_list = SE_list2,
#'                    tau_mu_log = TRUE, tau_delta_log = TRUE)
#'
set_initial_params <- function(sumstat_beta_list, sumstat_se_list, tau_mu_log = FALSE, tau_delta_log = FALSE,
                               use_tau_delta = FALSE) {

  beta_matrix <- matrix(unlist(sumstat_beta_list), nrow = length(sumstat_beta_list), byrow = TRUE)
  se_matrix <- matrix(unlist(sumstat_se_list), nrow = length(sumstat_beta_list), byrow = TRUE)

  #gamma - simple IVW regression of outcome/exposure from the target population
  init_gamma <- unname(stats::lm(beta_matrix[,1] ~ beta_matrix[,2] - 1, weights = 1/se_matrix[,1]^2)$coefficients)

  #the covariance of exposure effect size estimates for each population
  beta_x_vars <- stats::var(beta_matrix[,-1])

  if (length(beta_x_vars) == 1) {
    #tau_mu and tau_delta are not identifiable in the presence of no auxiliary populations.
      #This means you can only estimate one.
    # ifelse(use_tau_delta, print("No auxiliary populations, only estimating tau_delta"),
    #                       print("No auxiliary populations, only estimating tau_mu"))
    init_tau_mu <- ifelse(use_tau_delta, NA, beta_x_vars)

    init_tau_delta <- ifelse(use_tau_delta, beta_x_vars, NA)

    return(c(gamma = init_gamma,
             tau_mu = ifelse(tau_mu_log, log(init_tau_mu), init_tau_mu),
             tau_delta = ifelse(tau_delta_log, log(init_tau_delta), init_tau_delta)))
  }

  #tau_mu, estimated as the average off-diagonal value of this covariance matrix
  init_tau_mu <- unname(mean(c(beta_x_vars[upper.tri(beta_x_vars, diag = FALSE)], beta_x_vars[lower.tri(beta_x_vars, diag = FALSE)])))
  if (init_tau_mu <= 0) {init_tau_mu <- 1e-16} #occasionally may fall under zero


  #when variants are standardized, tau_delta should be easy to estimate because SE is the same between variants
  se_x_squared <- se_matrix[1,-1]^2
  tau_delta_estimates <- unname(diag(beta_x_vars) - init_tau_mu - se_x_squared)

  init_tau_delta <- mean(tau_delta_estimates)
  if (init_tau_delta <= 0) {init_tau_delta <- 1e-16}

  return(c(gamma = init_gamma,
           tau_mu = ifelse(tau_mu_log, log(init_tau_mu), init_tau_mu),
           tau_delta = ifelse(tau_delta_log, log(init_tau_delta), init_tau_delta)))
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
#' @param optim_method what method to use for the optim function (Nelder-Mead (default), Brent (for 1-dimensional optimization), or L-BFGS-B (not recommended))
#' @param set.init.params Whether to manually set initial parameters from the initial data. If not, the initial parameters used are (0, 0, 0).
#' @param select_zscore if the instrument was selected because its association with the exposure in the target exceeded some Z-score threshold, what threshold this is. If no selection was performed to obtain the summary statistics, put NA.
#'
#' @returns the output from optim(), detailing the optimized values for gamma, tau_mu and tau_delta, information about convergence, and the negative log-likelihood
#' @export
#'
#' @examples
#' SE_list <- rep(list(c(0.2, 0.1, 0.05, 0.05)), 50)
#' r_mat <- diag(4)
#' r_mat[1,2] <- 0.2
#' r_mat[2,1] <- 0.2
#' observed_data <- simplemodel_sim(gamma = 0.7, tau_mu = 0.5, tau_delta = 0.2, SE_list = SE_list,
#'                                  vars = 50, pops = 3, r_mat = r_mat)
#' sumstat_beta_list <- apply(observed_data$beta_matrix, MARGIN = 1, function(x) {return(x)},
#'                            simplify = FALSE)
#' simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list, sumstat_se_list = SE_list,
#'                        r_mat_list = rep(list(r_mat), 50),
#'                        tau_mu_log = TRUE, tau_delta_log = TRUE, optim_method = "Nelder-Mead")
#' simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list, sumstat_se_list = SE_list,
#'                        r_mat_list = rep(list(r_mat), 50),
#'                        tau_mu_log = TRUE, tau_delta_log = TRUE, optim_method = "Nelder-Mead",
#'                        set.init.params = TRUE)
#' simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list, sumstat_se_list = SE_list,
#'                        r_mat_list = rep(list(r_mat), 50),
#'                        tau_mu_log = TRUE, tau_delta_log = TRUE, optim_method = "L-BFGS-B")
#'
#' observed_data <- simplemodel_sim(gamma = 0.7, tau_mu = 0.5, tau_delta = 0, SE_list = SE_list,
#'                                  vars = 50, pops = 3, r_mat = r_mat)
#' sumstat_beta_list <- apply(observed_data$beta_matrix, MARGIN = 1, function(x) {return(x)},
#'                            simplify = FALSE)
#' simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list, sumstat_se_list = SE_list,
#'                        r_mat_list = rep(list(r_mat), 50),
#'                        is.fixed = c(FALSE, FALSE, TRUE), fix.params = c(NA, NA, 0),
#'                        tau_mu_log = TRUE, tau_delta_log = FALSE, optim_method = "Nelder-Mead")
#' simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list, sumstat_se_list = SE_list,
#'                        r_mat_list = rep(list(r_mat), 50),
#'                        is.fixed = c(FALSE, FALSE, TRUE), fix.params = c(NA, NA, 0),
#'                        tau_mu_log = TRUE, tau_delta_log = FALSE, optim_method = "L-BFGS-B")
simple_loglik_optimize <- function(sumstat_beta_list, sumstat_se_list, is_overlap = FALSE, r_mat_list= NA, is.fixed = c(FALSE, FALSE, FALSE), fix.params = c(NA, NA, NA), set.init.params = FALSE, tau_mu_log = FALSE, tau_delta_log = FALSE, optim_method = c("Nelder-Mead", "L-BFGS-B", "Brent"), select_zscore = NA) {

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
  if (set.init.params) {
    #If tau_delta is set to fixed, set initial parameters estimating tau_delta (use_tau_delta = FALSE); otherwise, estimate using tau_delta
    use_tau_delta <- !is.fixed[3]
    init.params <- set_initial_params(sumstat_beta_list, sumstat_se_list, tau_mu_log, tau_delta_log, use_tau_delta = use_tau_delta)
    # names(init.params) <- c("gamma", "tau_mu", "tau_delta")
  } else {
    init.params <- c(gamma = 0, tau_mu = 0, tau_delta = 0)
  }
  init.params <- init.params[!is.fixed]

  #Brent's method needs a different function to work
  if (optim_method == "Brent"){
    free_names <- names(init.params)
    if (length(free_names) > 1) {
      stop("Brent's method won't work with >1 dimension!")
    }
    free_index <- which(c("gamma", "tau_mu", "tau_delta") == free_names)
    optimize_fn <- function(x) {
      full <- fix.params
      full[free_index] <- x
      gamma <- full[1]
      tau_mu <- full[2]
      tau_delta <- full[3]
      return(-simple_loglik_full(sumstat_beta_list = sumstat_beta_list,
                                 sumstat_se_list = sumstat_se_list,
                                 gamma = gamma,
                                 is_overlap = is_overlap,
                                 r_mat_list = r_mat_list,
                                 tau_mu = tau_mu,
                                 tau_delta = tau_delta,
                                 tau_mu_log = tau_mu_log,
                                 tau_delta_log = tau_delta_log,
                                 select_zscore = select_zscore))
    }
  } else {
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
                                 tau_delta_log = tau_delta_log,
                                 select_zscore = select_zscore))
    }
  }

  #minimize the negative log-likelihood using Brent's method, Nelder-Mead or L-BFGS-B
  if (optim_method == "Brent") {#for (bounded) one-dimensional optimization using Brent's method
    if ((is.fixed[2] == FALSE & tau_mu_log == FALSE) ||
        (is.fixed[3] == FALSE & tau_delta_log == FALSE)) { #optimizing a non-transformed variance component, bounds should be positive
      stats::optim(par = init.params[free_names],
                   fn = optimize_fn,
                   method = optim_method,
                   lower = 0, upper = 10,
                   hessian = TRUE)
    } else { #optimizing gamma or a log-transformed variance component, search space centered around zero
      stats::optim(par = init.params[free_names],
                   fn = optimize_fn,
                   method = optim_method,
                   lower = -16, upper = 16,
                   hessian = TRUE)
    }
  } else {
    if (optim_method == "L-BFGS-B") { #L-BFGS-B needs bounds
      stats::optim(par = init.params,
                   fn = optim_fn,
                   method = optim_method,
                   lower = c(-Inf, -16, -16)[!is.fixed], upper = c(Inf, 10, 10)[!is.fixed],
                   hessian = TRUE)
    } else { #Nelder-Mead is unbounded and is the default method if a "wrong" optimization method is specified
      stats::optim(par = init.params,
                   fn = optim_fn,
                   method = optim_method,
                   hessian = TRUE)
    }
  }
}
