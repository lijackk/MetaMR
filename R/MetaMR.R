#' MetaMR using a simple random-effects meta-analysis model
#'
#' This function estimates an exposure-outcome causal effect, \eqn{\gamma}, in the target population, using exposure and outcome summary statistics in the target population as well as exposure summary statistics from \eqn{K - 1} auxiliary populations using a simple hierarchical model. This function inputs a list of GWAS effect size estimates across all \eqn{K+1} studies for a set of \eqn{n} variants selected as instrumental variables, a list of standard errors for these variants, and a user-specified residual correlation matrix reflecting sample overlap between studies, then maximizes a likelihood that models exposure summary statistics for each variant \eqn{j} under a random effects meta-analysis framework based on central values \eqn{\mu_j} and population-specific deviations \eqn{\delta_{kj}}, where \eqn{\mu_j} are i.i.d. normal with variance \eqn{\tau_{\mu}} and \eqn{\delta_{kj}} are also i.i.d. normal with variance \eqn{\tau_{\delta}}. Note that this function does NOT incorporate variant selection, potential horizontal pleiotropy in the outcome, or alternative distributions for \eqn{\mu} and \eqn{\delta} values.
#'
#' @param sumstat_beta_list a list of vectors of GWAS effect size estimates, each with length K+1: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations. Each vector in the list represents summary statistics for one of N variants.
#' @param sumstat_se_list a list of vectors of the standard errors for GWAS effect size estimates in sumstat_beta_list
#' @param is_overlap A boolean describing whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat_list a list of matrices of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap.
#' @param check_zeros whether to check if either variance parameter (tau_mu) or (tau_delta) was estimated to be exactly zero, which can destabilize the algorithm. If either estimated variance parameter falls below a certain threshold, we rerun optimization fixing that component at exactly zero.
#' @param zero_thresh the threshold specified in check_zeros. Is 1e-8 by default.
#' @param tau_mu_zero Whether tau_mu should be fixed at zero when performing MetaMR.
#' @param tau_delta_zero Whether tau_delta should be fixed at zero when performing MetaMR.
#' @param set.init.params Whether the optim() function should use the observed data to set initial parameters for optimization.
#' @param select_zscore if the instrument was selected because its association with the exposure in the target exceeded some Z-score threshold, what threshold this is. If no selection was performed to obtain the summary statistics, put NA.
#'
#' @returns a list that includes point estimates for the exposure-outcome causal effect \eqn{\gamma} and variance terms \eqn{\tau_{\mu}} and \eqn{\tau_{\delta}} estimated using maximum likelihood, the covariance matrix of these point estimates (calculated using the inverse of the Hessian), and the results of a likelihood ratio test that tests against null hypothesis \eqn{\gamma = 0} (including the test statistic and p-value)
#' @export
#'
#' @examples
#' SE_list <- rep(list(c(0.2, 0.1, 0.05, 0.05)), 50)
#' r_mat <- diag(4)
#' r_mat[1,2] <- 0.2
#' r_mat[2,1] <- 0.2
#' observed_data <- simplemodel_sim(gamma = 0.7, tau_mu = 0.05, tau_delta = 0.02, SE_list = SE_list,
#'                                  vars = 50, pops = 3, r_mat = r_mat)
#' sumstat_beta_list <- apply(observed_data$beta_matrix, MARGIN = 1, function(x) {return(x)},
#'                            simplify = FALSE)
#' MetaMR_simplemodel(sumstat_beta_list = sumstat_beta_list, sumstat_se_list = SE_list,
#'                    r_mat_list = rep(list(r_mat), 50))
#'
MetaMR_simplemodel <- function(sumstat_beta_list, sumstat_se_list, is_overlap = FALSE, r_mat_list= NA, check_zeros = TRUE, zero_thresh = 1e-8, set.init.params = FALSE, tau_mu_zero = FALSE, tau_delta_zero = FALSE, select_zscore = NA) {
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

  #Keeping track of which parameters to be fixed or not fixed throughout
  is.fixed <- c(FALSE, tau_mu_zero, tau_delta_zero)
  fix.params <- c(NA, ifelse(tau_mu_zero, 0, NA), ifelse(tau_delta_zero, 0, NA))
  tau_mu_log <- ifelse(tau_mu_zero, FALSE, TRUE)
  tau_delta_log <-ifelse(tau_delta_zero, FALSE, TRUE)

  #if K = 1, then variance parameters tau_mu and tau_delta are not identifiable.
    #If both variance parameters are not fixed at zero, we force tau_delta to zero.
  if (k == 1) {
    if (!tau_mu_zero & !tau_delta_zero) {
      print("No auxiliary populations, cannot estimate both variance parameters - tau_delta set to zero!")
      tau_delta_zero <- TRUE
      tau_delta_log <- FALSE
      is.fixed[3] <- TRUE
      fix.params[3] <- 0
    }
    #Optimizing the full likelihood
    MetaMR_point_est <- simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                               sumstat_se_list = sumstat_se_list,
                                               is_overlap = is_overlap,
                                               r_mat_list = r_mat_list,
                                               is.fixed = is.fixed,
                                               fix.params = fix.params,
                                               tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                               optim_method = ifelse(sum(is.fixed) >= 2, "Brent", "Nelder-Mead"),
                                               set.init.params = set.init.params,
                                               select_zscore = select_zscore)

    #Optimizing the constrained likelihood under the null hypothesis
    MetaMR_null_loglik <- simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                                 sumstat_se_list = sumstat_se_list,
                                                 is_overlap = is_overlap,
                                                 r_mat_list = r_mat_list,
                                                 is.fixed = c(TRUE, is.fixed[2:3]),
                                                 fix.params = c(0, fix.params[2:3]),
                                                 tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                                 optim_method = ifelse(sum(is.fixed[2:3]) >= 1, "Brent", "Nelder-Mead"),
                                                 set.init.params = set.init.params,
                                                 select_zscore = select_zscore)

  } else {
    #Optimizing the full likelihood
    MetaMR_point_est <- simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                               sumstat_se_list = sumstat_se_list,
                                               is_overlap = is_overlap,
                                               r_mat_list = r_mat_list,
                                               is.fixed = is.fixed,
                                               fix.params = fix.params,
                                               tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                               optim_method = ifelse(sum(is.fixed) >= 2, "Brent", "Nelder-Mead"),
                                               set.init.params = set.init.params,
                                               select_zscore = select_zscore)

    #Optimizing the constrained likelihood under the null hypothesis
    MetaMR_null_loglik <- simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                                 sumstat_se_list = sumstat_se_list,
                                                 is_overlap = is_overlap,
                                                 r_mat_list = r_mat_list,
                                                 is.fixed = c(TRUE, is.fixed[2:3]),
                                                 fix.params = c(0, fix.params[2:3]),
                                                 tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                                 optim_method = ifelse(sum(is.fixed[2:3]) >= 1, "Brent", "Nelder-Mead"),
                                                 set.init.params = set.init.params,
                                                 select_zscore = select_zscore)
  }

  if (check_zeros) { #checks whether tau_mu or tau_delta were optimized to be nearly zero, then reruns MetaMR if they were.

    #Which variance parameters were estimated
    tm_indicator <- "tau_mu" %in% names(MetaMR_point_est$par)
    td_indicator <- "tau_delta" %in% names(MetaMR_point_est$par)

    #If a variance parameter was estimated, was it "close enough" to zero?
    taumu_nearzero <- ifelse(!tm_indicator, FALSE, ifelse(exp(MetaMR_point_est$par["tau_mu"]) < zero_thresh, TRUE, FALSE))
    taudelta_nearzero <- ifelse(!td_indicator, FALSE, ifelse(exp(MetaMR_point_est$par["tau_delta"]) < zero_thresh, TRUE, FALSE))

    if (!taumu_nearzero & !taudelta_nearzero) { #neither parameter was estimated near zero; there is nothing to do
      print("no variance components estimated near zero by MetaMR.")
    } else {
      if (taumu_nearzero) {
        print("tau_mu estimated near zero; rerunning estimation fixed at exactly zero.")
        tau_mu_zero <- TRUE
        tau_mu_log <- FALSE
        is.fixed[2] <- TRUE
        fix.params[2] <- 0
      }

      if (taudelta_nearzero) {
        print("tau_delta estimated near zero; rerunning estimation fixed at exactly zero.")
        tau_delta_zero <- TRUE
        tau_delta_log <- FALSE
        is.fixed[3] <- TRUE
        fix.params[3] <- 0
      }
        #Optimizing the full likelihood with the new constraints
        MetaMR_point_est <- simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                                   sumstat_se_list = sumstat_se_list,
                                                   is_overlap = is_overlap,
                                                   r_mat_list = r_mat_list,
                                                   is.fixed = is.fixed,
                                                   fix.params = fix.params,
                                                   tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                                   optim_method = ifelse(sum(is.fixed) >= 2, "Brent", "Nelder-Mead"),
                                                   set.init.params = set.init.params,
                                                   select_zscore = select_zscore)

        #Optimizing the constrained likelihood under the null hypothesis
        MetaMR_null_loglik <- simple_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                                     sumstat_se_list = sumstat_se_list,
                                                     is_overlap = is_overlap,
                                                     r_mat_list = r_mat_list,
                                                     is.fixed = c(TRUE, is.fixed[2:3]),
                                                     fix.params = c(0, fix.params[2:3]),
                                                     tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                                     optim_method = ifelse(sum(is.fixed[2:3]) >= 1, "Brent", "Nelder-Mead"),
                                                     set.init.params = set.init.params,
                                                     select_zscore = select_zscore)
    }
  }

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

#' MetaMR using a simple random-effects meta-analysis model on a list of pre-selected variants
#'
#' Like MetaMR_simplemodel, but now assumes that the variants under consideration have passed some selection method. Currently supports variants being selected based on the p-value of exposure association in a single population, their minimum p-value of exposure association across all populations, or the p-value of Fisher's method when applied to the exposure summary statistics for that variant. This method currently does NOT check whether variants being considered fulfills variant selection, but this may be included in the future.
#'
#' @param sumstat_beta_list a list of vectors of GWAS effect size estimates, each with length K+1: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations. Each vector in the list represents summary statistics for one of N variants.
#' @param sumstat_se_list a list of vectors of the standard errors for GWAS effect size estimates in sumstat_beta_list
#' @param is_overlap A boolean describing whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat_list a list of matrices of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap.
#' @param check_zeros whether to check if either variance parameter (tau_mu) or (tau_delta) was estimated to be exactly zero, which can destabilize the algorithm. If either estimated variance parameter falls below a certain threshold, we rerun optimization fixing that component at exactly zero.
#' @param zero_thresh the threshold specified in check_zeros. Is 1e-8 by default.
#' @param set.init.params Whether the optim() function should use the observed data to set initial parameters for optimization.
#' @param tau_mu_zero Whether tau_mu should be fixed at zero when performing MetaMR.
#' @param tau_delta_zero Whether tau_delta should be fixed at zero when performing MetaMR.
#' @param select_method The method by which variants were selected: "single_exposure" = a single population's exposure association p-value, "minp" = minimum exposure association p-value across all populations, "fisher" = Fisher's method
#' @param single_exp_pop If select_method = "single_exposure", what exposure population to apply it to
#' @param select_pthresh The p-value threshold for variant selection
#' @param resample_number If select_method = "fisher", the number of resamples performed to estimate the probability of a variant being selected under Fisher's method.
#' @param resample_seed The seed used for resampling
#' @param standardized A boolean indicating whether summary statistics are standardized; i.e. whether they have the same set of standard errors. If they do, we only need to calculate the probability of selection once across all variants.
#' @param verbose A boolean indicating whether to print out current values of parameters at every iteration of optim(). Use for troubleshooting purposes.
#'
#' @returns a list that includes point estimates for the exposure-outcome causal effect \eqn{\gamma} and variance terms \eqn{\tau_{\mu}} and \eqn{\tau_{\delta}} estimated using maximum likelihood, the covariance matrix of these point estimates (calculated using the inverse of the Hessian), and the results of a likelihood ratio test that tests against null hypothesis \eqn{\gamma = 0} (including the test statistic and p-value)
#' @export
#'
#' @examples
#' nvar_causal <- 200
#' gamma <- 0.1
#' h2 <- 0.001*nvar_causal
#' persnp_h2 <- h2 / nvar_causal
#' tau_mu <- 0.9*persnp_h2
#' tau_delta <- persnp_h2 - tau_mu
#' N_k <- c(100000, 5000, 50000)
#' sigma_Xk <- 1/sqrt(N_k)
#' special_conditions <- "selection_allvars"
#' SE_list_causal <- rep(list(sigma_Xk), nvar_causal)
#' nrep <- 200
#' seed <- 2025
#' observed_data_causal <- simplemodel_sim(gamma = gamma, tau_mu = tau_mu,
#'                                        tau_delta = tau_delta,
#'                                        SE_list = SE_list_causal,
#'                                        vars = nvar_causal,
#'                                        pops = 2,
#'                                        seed = seed)
#' sumstat_beta_list_causal <- apply(observed_data_causal$beta_matrix, MARGIN = 1,
#'                                   function(x){return(x)}, simplify = FALSE)
#' crit.p <- 0.05/nvar_causal
#' select_minp <- which(abs(observed_data_causal$beta_matrix[,2]/observed_data_causal$se_matrix[,2]) >
#'                      abs(qnorm(crit.p/2)) |
#' abs(observed_data_causal$beta_matrix[,3]/observed_data_causal$se_matrix[,3]) >
#' abs(qnorm(crit.p/2)))
#' sumstat_beta_list_causal <- apply(observed_data_causal$beta_matrix, MARGIN = 1,
#'                                   function(x) {return(x)}, simplify = FALSE)
#' minp_metamr <- MetaMR_selectionmodel(sumstat_beta_list = sumstat_beta_list_causal[select_minp],
#'                                      sumstat_se_list = SE_list_causal[select_minp],
#'                                      is_overlap = FALSE, r_mat_list= NA,
#'                                      check_zeros = TRUE, zero_thresh = 1e-8,
#'                                      set.init.params = FALSE,
#'                                      tau_mu_zero = FALSE, tau_delta_zero = FALSE,
#'                                      select_method = "minp",
#'                                      single_exp_pop = 1, select_pthresh = 0.05,
#'                                      resample_number = 1e5, resample_seed = 2025,
#'                                      standardized = TRUE, verbose = FALSE)
MetaMR_selectionmodel <- function(sumstat_beta_list, sumstat_se_list, is_overlap = FALSE, r_mat_list= NA, check_zeros = TRUE, zero_thresh = 1e-8, set.init.params = FALSE, tau_mu_zero = FALSE, tau_delta_zero = FALSE, select_method = c("single_exposure", "minp", "fisher"), single_exp_pop = 1, select_pthresh = 0.05, resample_number = 1e5, resample_seed = 2025, standardized = FALSE, verbose = FALSE) {

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

  #Keeping track of which parameters to be fixed or not fixed throughout
  is.fixed <- c(FALSE, tau_mu_zero, tau_delta_zero)
  fix.params <- c(NA, ifelse(tau_mu_zero, 0, NA), ifelse(tau_delta_zero, 0, NA))
  tau_mu_log <- ifelse(tau_mu_zero, FALSE, TRUE)
  tau_delta_log <-ifelse(tau_delta_zero, FALSE, TRUE)

  #if K = 1, then variance parameters tau_mu and tau_delta are not identifiable.
  #If both variance parameters are not fixed at zero, we force tau_delta to zero.
  if (k == 1) {
    if (!tau_mu_zero & !tau_delta_zero) {
      print("No auxiliary populations, cannot estimate both variance parameters - tau_delta set to zero!")
      tau_delta_zero <- TRUE
      tau_delta_log <- FALSE
      is.fixed[3] <- TRUE
      fix.params[3] <- 0
    }
    #Optimizing the full likelihood
    MetaMR_point_est <- selection_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                               sumstat_se_list = sumstat_se_list,
                                               is_overlap = is_overlap,
                                               r_mat_list = r_mat_list,
                                               is.fixed = is.fixed,
                                               fix.params = fix.params,
                                               tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                               optim_method = ifelse(sum(is.fixed) >= 2, "Brent", "Nelder-Mead"),
                                               set.init.params = set.init.params,
                                               select_method = select_method,
                                               single_exp_pop = single_exp_pop,
                                               select_pthresh = select_pthresh,
                                               resample_number = resample_number,
                                               resample_seed = resample_seed,
                                               standardized = standardized,
                                               verbose = verbose)

    #Optimizing the constrained likelihood under the null hypothesis
    MetaMR_null_loglik <- selection_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                                 sumstat_se_list = sumstat_se_list,
                                                 is_overlap = is_overlap,
                                                 r_mat_list = r_mat_list,
                                                 is.fixed = c(TRUE, is.fixed[2:3]),
                                                 fix.params = c(0, fix.params[2:3]),
                                                 tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                                 optim_method = ifelse(sum(is.fixed[2:3]) >= 1, "Brent", "Nelder-Mead"),
                                                 set.init.params = set.init.params,
                                                 select_method = select_method,
                                                 single_exp_pop = single_exp_pop,
                                                 select_pthresh = select_pthresh,
                                                 resample_number = resample_number,
                                                 resample_seed = resample_seed,
                                                 standardized = standardized,
                                                 verbose = verbose)

  } else {
    #Optimizing the full likelihood
    MetaMR_point_est <- selection_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                               sumstat_se_list = sumstat_se_list,
                                               is_overlap = is_overlap,
                                               r_mat_list = r_mat_list,
                                               is.fixed = is.fixed,
                                               fix.params = fix.params,
                                               tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                               optim_method = ifelse(sum(is.fixed) >= 2, "Brent", "Nelder-Mead"),
                                               set.init.params = set.init.params,
                                               select_method = select_method,
                                               single_exp_pop = single_exp_pop,
                                               select_pthresh = select_pthresh,
                                               resample_number = resample_number,
                                               resample_seed = resample_seed,
                                               standardized = standardized,
                                               verbose = verbose)

    #Optimizing the constrained likelihood under the null hypothesis
    MetaMR_null_loglik <- selection_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                                 sumstat_se_list = sumstat_se_list,
                                                 is_overlap = is_overlap,
                                                 r_mat_list = r_mat_list,
                                                 is.fixed = c(TRUE, is.fixed[2:3]),
                                                 fix.params = c(0, fix.params[2:3]),
                                                 tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                                 optim_method = ifelse(sum(is.fixed[2:3]) >= 1, "Brent", "Nelder-Mead"),
                                                 set.init.params = set.init.params,
                                                 select_method = select_method,
                                                 single_exp_pop = single_exp_pop,
                                                 select_pthresh = select_pthresh,
                                                 resample_number = resample_number,
                                                 resample_seed = resample_seed,
                                                 standardized = standardized,
                                                 verbose = verbose)
  }

  if (check_zeros) { #checks whether tau_mu or tau_delta were optimized to be nearly zero, then reruns MetaMR if they were.

    #Which variance parameters were estimated
    tm_indicator <- "tau_mu" %in% names(MetaMR_point_est$par)
    td_indicator <- "tau_delta" %in% names(MetaMR_point_est$par)

    #If a variance parameter was estimated, was it "close enough" to zero?
    taumu_nearzero <- ifelse(!tm_indicator, FALSE, ifelse(exp(MetaMR_point_est$par["tau_mu"]) < zero_thresh, TRUE, FALSE))
    taudelta_nearzero <- ifelse(!td_indicator, FALSE, ifelse(exp(MetaMR_point_est$par["tau_delta"]) < zero_thresh, TRUE, FALSE))

    if (!taumu_nearzero & !taudelta_nearzero) { #neither parameter was estimated near zero; there is nothing to do
      print("no variance components estimated near zero by MetaMR.")
    } else {
      if (taumu_nearzero) {
        print("tau_mu estimated near zero; rerunning estimation fixed at exactly zero.")
        tau_mu_zero <- TRUE
        tau_mu_log <- FALSE
        is.fixed[2] <- TRUE
        fix.params[2] <- 0
      }

      if (taudelta_nearzero) {
        print("tau_delta estimated near zero; rerunning estimation fixed at exactly zero.")
        tau_delta_zero <- TRUE
        tau_delta_log <- FALSE
        is.fixed[3] <- TRUE
        fix.params[3] <- 0
      }
      #Optimizing the full likelihood with the new constraints
      MetaMR_point_est <- selection_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                                 sumstat_se_list = sumstat_se_list,
                                                 is_overlap = is_overlap,
                                                 r_mat_list = r_mat_list,
                                                 is.fixed = is.fixed,
                                                 fix.params = fix.params,
                                                 tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                                 optim_method = ifelse(sum(is.fixed) >= 2, "Brent", "Nelder-Mead"),
                                                 set.init.params = set.init.params,
                                                 select_method = select_method,
                                                 single_exp_pop = single_exp_pop,
                                                 select_pthresh = select_pthresh,
                                                 resample_number = resample_number,
                                                 resample_seed = resample_seed,
                                                 standardized = standardized,
                                                 verbose = verbose)

      #Optimizing the constrained likelihood under the null hypothesis
      MetaMR_null_loglik <- selection_loglik_optimize(sumstat_beta_list = sumstat_beta_list,
                                                   sumstat_se_list = sumstat_se_list,
                                                   is_overlap = is_overlap,
                                                   r_mat_list = r_mat_list,
                                                   is.fixed = c(TRUE, is.fixed[2:3]),
                                                   fix.params = c(0, fix.params[2:3]),
                                                   tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                                   optim_method = ifelse(sum(is.fixed[2:3]) >= 1, "Brent", "Nelder-Mead"),
                                                   set.init.params = set.init.params,
                                                   select_method = select_method,
                                                   single_exp_pop = single_exp_pop,
                                                   select_pthresh = select_pthresh,
                                                   resample_number = resample_number,
                                                   resample_seed = resample_seed,
                                                   standardized = standardized,
                                                   verbose = verbose)
    }
  }

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
