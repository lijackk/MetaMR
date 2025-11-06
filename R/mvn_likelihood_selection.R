#' Log-likelihood function under simple model with variant selection (single-variant)
#'
#' A conditional log-likelihood function for our simple model that accounts for variants being pre-selected using either a single population's exposure association p-value, minimum exposure association p-value across all populations, or Fisher's method.
#'
#' @param sumstat_beta the GWAS effect size estimates for a single variant across K+1 studies: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations.
#' @param sumstat_se the standard errors of the effect size estimates in sumstat_beta
#' @param gamma The exposure-outcome causal effect parameter
#' @param is_overlap a boolean whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat a matrix of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap.
#' @param tau_mu The variance of the normal distribution underlying mu. It should be log-transformed unless set to be exactly zero.
#' @param tau_delta The variance of the normal distributions underlying delta. It should be log-transformed unless set to be exactly zero.
#' @param tau_mu_log Whether tau_mu represents its log-transformed value or not
#' @param tau_delta_log Whether tau_delta represents its log-transformed value or not
#' @param select_method The method by which variants were selected: "single_exposure" = a single population's exposure association p-value, "minp" = minimum exposure association p-value across all populations, "fisher" = Fisher's method
#' @param single_exp_pop If select_method = "single_exposure", what exposure population to apply it to
#' @param select_pthresh The p-value threshold for variant selection
#' @param resample_number If select_method = "fisher", the number of resamples performed to estimate the probability of a variant being selected under Fisher's method.
#' @param resample_seed The seed used for resampling
#'
#' @returns A log-likelihood for the summary statistics of a single (selected) variant given some parameters gamma, mu, tau
#' @export
selection_loglik_single <- function(sumstat_beta,
                                    sumstat_se,
                                    gamma,
                                    is_overlap = FALSE,
                                    r_mat = NA,
                                    tau_mu,
                                    tau_delta,
                                    tau_mu_log = FALSE,
                                    tau_delta_log = FALSE,
                                    select_method = c("single_exposure", "minp", "fisher"),
                                    single_exp_pop = 1,
                                    select_pthresh = 0.05,
                                    resample_number = 1e5,
                                    resample_seed = 2025) {
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

  #calculate the full likelihood first
  unconditional_loglik <- simple_loglik_single(sumstat_beta,
                                               sumstat_se,
                                               gamma,
                                               is_overlap,
                                               r_mat,
                                               tau_mu,
                                               tau_delta,
                                               tau_mu_log,
                                               tau_delta_log)

  #Calculating the probability of selection

  #first, we need to undo the transformation of tau_mu and tau_delta, if it exists
  unt_tau_mu <- ifelse(tau_mu_log == TRUE, exp(tau_mu), tau_mu)
  unt_tau_delta <- ifelse(tau_delta_log == TRUE, exp(tau_delta), tau_delta)

  #we also need a diagonal matrix placeholder for r_mat, if it doesn't exist
  if(identical(r_mat, NA)) {
    r_mat <- diag(1, length(sumstat_beta))
  }

  if (select_method == "single_exposure") { #the probability of a variant being selected based on a single population strength of exposure association
    select_zscore <- abs(stats::qnorm(select_pthresh / 2))
    log_pselect <- log(2 * stats::pnorm(abs(select_zscore) * sumstat_se[1 + single_exp_pop]/sqrt(unt_tau_mu + unt_tau_delta + sumstat_se[1 + single_exp_pop]^2), lower.tail = FALSE))
  }

  if (select_method == "minp") { #the probability of a variant being selected based on the minimum p-value of exposure association across all populations
    log_pselect <- log(numint_prob(unt_tau_mu, unt_tau_delta, r_mat = r_mat[-1, -1], SE_vector = sumstat_se[-1], select_pthresh))
  }

  if (select_method == "fisher") { #the probability of a variant being selected based on Fisher's method of combining p-values
    log_pselect <- log(resample_prob(resample_number, resample_seed, unt_tau_mu, unt_tau_delta, r_mat = r_mat[-1, -1], SE_vector = sumstat_se[-1], select_pthresh))
  }

  #conditional log-likelihood = unconditional log-likelihood minus log probability of selection
  return(unconditional_loglik - log_pselect)
}

#' Calculating the probability of a variant being selected using resampling (Fisher's method)
#'
#' @param resample_number the number of resamples performed to estimate the probability of a variant being selected under Fisher's method.
#' @param resample_seed The seed used for resampling
#' @param tau_mu The variance of the normal distribution underlying mu.
#' @param tau_delta The variance of the normal distributions underlying delta.
#' @param r_mat the matrix from r_mat in selection_loglik_single, but for only exposure summary statistics (i.e. first row and column removed)
#' @param SE_vector the vector from sumstat_se in selection_loglik_single, but for only exposure summary statistics (i.e. first entry removed)
#' @param select_pthresh The p-value threshold for variant selection
#'
#' @returns The probability that a variant was selected under a Fisher's method selection strategy
#' @export
resample_prob <- function(resample_number, resample_seed, tau_mu, tau_delta, r_mat, SE_vector, select_pthresh) {
  pops <- length(SE_vector)
  if (!is.na(resample_seed)) {
    set.seed(resample_seed)
  }

  #Empirically assessing the probability of selection - having a seed
  #because the number of resamples can get very large, I'll just use a pure for loop for everything

  # The covariance matrix stays the same for all loops
  SE_matrix <- r_mat * outer(SE_vector, SE_vector)
  select_indicator <- vector()

  #Establishing thresholds for either selection method once so we don't have to do it multiple times
  fisher.crit <- stats::qchisq(select_pthresh, df = 2*pops, lower.tail = FALSE)


  for (i in 1:resample_number) {
    #Mean vector
    mu <- stats::rnorm(n = 1, mean = 0, sd = sqrt(tau_mu))
    mean_vector <- mu + stats::rnorm(n = pops, mean = 0, sd = sqrt(tau_delta))

    #Generating random draws of observed exposure association data based on our parameters
    sumstat_vector <- mvtnorm::rmvnorm(n = 1, mean = mean_vector, sigma = SE_matrix)
    zscores <- sumstat_vector/SE_vector

    logpvals <- log(stats::pnorm(-abs(zscores)) * 2)
    select_indicator[i] <- -2 * sum(logpvals) > fisher.crit
  }

  select_prob <- ifelse(sum(select_indicator) > 0, sum(select_indicator)/resample_number, 0.5/resample_number)
  return(select_prob)
}

#' Calculating the probability of a variant being selected using numerical integration (multi-exposure minimum p-value)
#'
#' @param tau_mu The variance of the normal distribution underlying mu.
#' @param tau_delta The variance of the normal distributions underlying delta.
#' @param r_mat the matrix from r_mat in selection_loglik_single, but for only exposure summary statistics (i.e. first row and column removed)
#' @param SE_vector the vector from sumstat_se in selection_loglik_single, but for only exposure summary statistics (i.e. first entry removed)
#' @param select_pthresh The p-value threshold for variant selection
#'
#' @returns The probability that a variant was selected under a minimum exposure association p-value selection strategy
#' @export
numint_prob <- function(tau_mu, tau_delta, r_mat, SE_vector, select_pthresh) {
  select_zscore <- abs(stats::qnorm(select_pthresh/2))
  int_lb <- -select_zscore*SE_vector
  int_ub <- select_zscore*SE_vector

  pops <- length(SE_vector)

  sigma_Ej <- r_mat * outer(SE_vector, SE_vector)
  sigma_G <- diag(tau_delta, pops) + matrix(tau_mu, nrow = pops, ncol = pops)
  sigma_j <- sigma_Ej + sigma_G

  #the integral function is the PDF of MVN
  p_cube <- mvtnorm::pmvnorm(lower = int_lb, upper = int_ub, mean = rep(0, pops), sigma = sigma_j, keepAttr = FALSE)

  return(1 - p_cube)
}

#' Log-likelihood function under simple model under variant selection (multi-variant)
#'
#' A conditional log-likelihood function for our simple model that accounts for variants being pre-selected using either a single population's exposure association p-value, minimum exposure association p-value across all populations, or Fisher's method.
#' @param sumstat_beta_list a list of vectors of GWAS effect size estimates with length K+1: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations. Each vector in the list represents summary statistics for one of many variants.
#' @param sumstat_se_list a list of vectors of the standard errors for GWAS effect size estimates in sumstat_beta_list
#' @param gamma The exposure-outcome causal effect parameter
#' @param is_overlap a boolean describing whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat_list a list of matrices of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap
#' @param tau_mu The variance of the normal distribution underlying mu. It should be log-transformed unless set to be exactly zero.
#' @param tau_delta The variance of the normal distributions underlying delta. It should be log-transformed unless set to be exactly zero.
#' @param tau_mu_log Whether tau_mu represents its log-transformed value or not.
#' @param tau_delta_log Whether tau_delta represents its log-transformed value or not.
#' @param select_method The method by which variants were selected: "single_exposure" = a single population's exposure association p-value, "minp" = minimum exposure association p-value across all populations, "fisher" = Fisher's method
#' @param single_exp_pop If select_method = "single_exposure", what exposure population to apply it to
#' @param select_pthresh The p-value threshold for variant selection
#' @param resample_number If select_method = "fisher", the number of resamples performed to estimate the probability of a variant being selected under Fisher's method.
#' @param resample_seed The seed used for resampling
#' @param standardized A boolean indicating whether summary statistics are standardized; i.e. whether they have the same set of standard errors. If they do, we only need to calculate the probability of selection once across all variants.
#'
#' @returns A conditional log-likelihood for the summary statistics of a set of (pre-selected) variants given some parameters gamma, mu, tau
#' @export
#'
selection_loglik_full <- function(sumstat_beta_list,
                                  sumstat_se_list,
                                  gamma,
                                  is_overlap = FALSE,
                                  r_mat_list = NA,
                                  tau_mu,
                                  tau_delta,
                                  tau_mu_log = FALSE,
                                  tau_delta_log = FALSE,
                                  select_method = c("single_exposure", "minp", "fisher"),
                                  single_exp_pop = 1,
                                  select_pthresh = 0.05,
                                  resample_number = 1e5,
                                  resample_seed = 2025,
                                  standardized = FALSE) {
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

  if (standardized) { #"standardized" simply means that every variant has the same set of standard errors, meaning that we only need to perform resampling/numerical integration once

    unconditional_loglik <- simple_loglik_full(sumstat_beta = sumstat_beta_list,
                                               sumstat_se = sumstat_se_list,
                                               is_overlap = TRUE,
                                               r_mat = r_mat_list,
                                               gamma = gamma,
                                               tau_mu = tau_mu,
                                               tau_delta = tau_delta,
                                               tau_mu_log = tau_mu_log,
                                               tau_delta_log = tau_delta_log)

    #If standardized, the standard errors provided for the first variant should be the same across all variants.

    #First, we untransform tau_mu and tau_delta
    unt_tau_mu <- ifelse(tau_mu_log == TRUE, exp(tau_mu), tau_mu)
    unt_tau_delta <- ifelse(tau_delta_log == TRUE, exp(tau_delta), tau_delta)

    #We then resample and calculate probabilities based on the first variant
    #For Fisher's method, we need to use resampling
    if (select_method == "fisher") {
      log_pselect <- log(resample_prob(resample_number, resample_seed, unt_tau_mu, unt_tau_delta, r_mat = r_mat_list[[1]][-1, -1], SE_vector = sumstat_se_list[[1]][-1], select_pthresh))
      return(unconditional_loglik - log_pselect*n)
    } else {
      if (select_method == "minp") {
        log_pselect <- log(numint_prob(unt_tau_mu, unt_tau_delta, r_mat = r_mat_list[[1]][-1, -1], SE_vector = sumstat_se_list[[1]][-1], select_pthresh))
        return(unconditional_loglik - log_pselect*n)
      } else {#if not Fisher's method or minimum p, default to single-exposure
        sumstat_se <- sumstat_se_list[[1]]
        select_zscore <- abs(stats::qnorm(select_pthresh / 2))
        log_pselect <- log(2 * stats::pnorm(abs(select_zscore) * sumstat_se[1 + single_exp_pop]/sqrt(unt_tau_mu + unt_tau_delta + sumstat_se[1 + single_exp_pop]^2), lower.tail = FALSE))
        return(unconditional_loglik - log_pselect*n)
      }
    }
  } else { #otherwise, we need to go the full way
    loglik_variants <- vector()
    for (i in 1:n) {
      loglik_variants[i] <- selection_loglik_single(sumstat_beta = sumstat_beta_list[[i]],
                                                    sumstat_se = sumstat_se_list[[i]],
                                                    is_overlap = TRUE,
                                                    r_mat = r_mat_list[[i]],
                                                    gamma = gamma,
                                                    tau_mu = tau_mu,
                                                    tau_delta = tau_delta,
                                                    tau_mu_log = tau_mu_log,
                                                    tau_delta_log = tau_delta_log,
                                                    select_method = select_method,
                                                    single_exp_pop = single_exp_pop,
                                                    select_pthresh = select_pthresh,
                                                    resample_number = resample_number,
                                                    resample_seed = resample_seed)
    }
    return(sum(loglik_variants))
  }
}

#' Optimizing the conditional log-likelihood under variant selection for a list of selected variants.
#'
#' @param sumstat_beta_list a list of vectors of GWAS effect size estimates with length K+1: first, for the outcome in the target population, then the exposure in the target population, then the exposures across K-1 auxiliary populations. Each vector in the list represents summary statistics for one of many variants.
#' @param sumstat_se_list a list of vectors of the standard errors for GWAS effect size estimates in sumstat_beta_list
#' @param is_overlap a boolean describing whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat_list a list of matrices of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap.
#' @param is.fixed which parameters of interested should be fixed during optimization. Follows order (gamma, tau_mu, tau_delta).
#' @param fix.params if a parameter is set as fixed in is.fixed, what it should be. Follows order (gamma, tau_mu, tau_delta). Our function will check to make sure is.fixed and fix.params are compatible with each other.
#' @param set.init.params Whether to manually set initial parameters from the initial data. If not, the initial parameters used are (0, 0, 0).
#' @param tau_mu_log whether tau_mu is log-transformed. It should be log-transformed unless set to be exactly zero.
#' @param tau_delta_log whether tau_delta is log-transformed. It should be log-transformed unless set to be exactly zero.
#' @param optim_method what method to use for the optim function (Nelder-Mead (default), Brent (for 1-dimensional optimization), or L-BFGS-B (not recommended))
#' @param select_method The method by which variants were selected: "single_exposure" = a single population's exposure association p-value, "minp" = minimum exposure association p-value across all populations, "fisher" = Fisher's method
#' @param single_exp_pop If select_method = "single_exposure", what exposure population to apply it to
#' @param select_pthresh The p-value threshold for variant selection
#' @param resample_number If select_method = "fisher", the number of resamples performed to estimate the probability of a variant being selected under Fisher's method.
#' @param resample_seed The seed used for resampling
#' @param standardized A boolean indicating whether summary statistics are standardized; i.e. whether they have the same set of standard errors. If they do, we only need to calculate the probability of selection once across all variants.
#' @param verbose A boolean indicating whether to print out current values of parameters at every iteration of optim(). Use for troubleshooting purposes.
#'
#' @returns the output from optim(), detailing the optimized values for gamma, tau_mu and tau_delta, information about convergence, and the negative log-likelihood
#' @export
selection_loglik_optimize <- function(sumstat_beta_list, sumstat_se_list,
                                      is_overlap = FALSE, r_mat_list= NA,
                                      is.fixed = c(FALSE, FALSE, FALSE), fix.params = c(NA, NA, NA), set.init.params = FALSE,
                                      tau_mu_log = FALSE, tau_delta_log = FALSE,
                                      optim_method = c("Nelder-Mead", "L-BFGS-B", "Brent"),
                                      select_method = c("single_exposure", "minp", "fisher"),
                                      single_exp_pop = 1,
                                      select_pthresh = 0.05,
                                      resample_number = 1e5,
                                      resample_seed = 2025,
                                      standardized = FALSE,
                                      verbose = FALSE) {
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
      return(-selection_loglik_full(sumstat_beta_list,
                                    sumstat_se_list,
                                    gamma,
                                    is_overlap = is_overlap,
                                    r_mat_list = r_mat_list,
                                    tau_mu,
                                    tau_delta,
                                    tau_mu_log = tau_mu_log,
                                    tau_delta_log = tau_delta_log,
                                    select_method = select_method,
                                    single_exp_pop = single_exp_pop,
                                    select_pthresh = select_pthresh,
                                    resample_number = resample_number,
                                    resample_seed = resample_seed,
                                    standardized = standardized))
    }
  } else {
    optim_fn <- function(params) {
      gamma <- ifelse("gamma" %in% names(init.params), params["gamma"], fix.params[1])
      tau_mu <- ifelse("tau_mu" %in% names(init.params), params["tau_mu"], fix.params[2])
      tau_delta <- ifelse("tau_delta" %in% names(init.params), params["tau_delta"], fix.params[3])

      #verbose prints out current values of parameters at every iteration. Use for troubleshooting purposes.
      if (verbose == TRUE) {
        print(paste("current values of gamma/tau_mu/tau_delta:", gamma, tau_mu, tau_delta))
      }

      return(-selection_loglik_full(sumstat_beta_list,
                                    sumstat_se_list,
                                    gamma,
                                    is_overlap = is_overlap,
                                    r_mat_list = r_mat_list,
                                    tau_mu,
                                    tau_delta,
                                    tau_mu_log = tau_mu_log,
                                    tau_delta_log = tau_delta_log,
                                    select_method = select_method,
                                    single_exp_pop = single_exp_pop,
                                    select_pthresh = select_pthresh,
                                    resample_number = resample_number,
                                    resample_seed = resample_seed,
                                    standardized = standardized))
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
      #the choice of lower bound depends on whether the tau parameters are log-transformed - upper bound remains unchanged.
      taumu_bound <- ifelse(tau_mu_log, -16, 1e-16)
      taudelta_bound <- ifelse(tau_delta_log, -16, 1e-16)
      stats::optim(par = init.params,
                   fn = optim_fn,
                   method = optim_method,
                   lower = c(-Inf, taumu_bound, taudelta_bound)[!is.fixed], upper = c(Inf, 10, 10)[!is.fixed],
                   hessian = TRUE)
    } else { #Nelder-Mead is unbounded and is the default method if a "wrong" optimization method is specified
      stats::optim(par = init.params,
                   fn = optim_fn,
                   method = optim_method,
                   hessian = TRUE)
    }
  }
}
