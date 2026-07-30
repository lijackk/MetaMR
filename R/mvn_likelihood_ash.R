#' Log-likelihood function incorporating variant selection and ASH pleiotropy (single-variant)
#'
#' A conditional log-likelihood function that models pleiotropy using the adaptive shrinkage (ASH) prior distribution proposed by Stephens et al (2015), accounting for variants being pre-selected using either a single population's exposure association p-value, minimum exposure association p-value across all populations, or Fisher's method. The ASH distribution is a mixture of normals with mean zero, and fixed variance components.
#'
#' @param sumstat_beta the GWAS effect size estimates for a single variant across K+1 studies: first, for the outcome in the target population, then across K exposure studies.
#' @param sumstat_se the standard errors of the effect size estimates in sumstat_beta
#' @param is_overlap a boolean whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat a matrix of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap.
#' @param kernel_matrix a K+1 x K+1 correlation matrix that captures correlations in exposure effect size heterogeneity across all populations (including the target population).
#' @param gamma The exposure-outcome causal effect
#' @param tau_mu The variance of the normal distribution underlying mu. May be log-transformed.
#' @param tau_delta The variance of the normal distributions underlying delta. May be log-transformed.
#' @param ash_pi A vector of length c representing the mixture probabilities of each of the components of the ASH distribution we use. By default, we assume c = 7 and that probabilities are equally distributed across all mixture variance parameters (1/7)
#' @param sigma_c A vector of length c representing the variance of each of the components of the ASH distribution we used. By default, we specify c = 7 with components 0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, and 1.
#' @param tau_mu_log Whether tau_mu is log-transformed
#' @param tau_delta_log Whether tau_delta is log-transformed
#' @param select_method The method by which variants were selected: "single_exposure" = a single population's exposure association p-value, "minp" = minimum exposure association p-value across all populations, "fisher" = Fisher's method.
#' @param single_exp_pop If select_method = "single_exposure", what exposure study it applies to. By default, we assume it is the first.
#' @param select_pthresh The p-value threshold for variant selection
#' @param random_seed For Fisher's method of selection, what random seed to use for sampling of observations
#' @param mc_iter  For Fisher's method of selection, how many random samples to take
#'
#' @returns A log-likelihood for the summary statistics of a single (selected) variant given parameters gamma, mu, tau, ash_pi
#' @export
#'
metamrash_loglik.single <- function(sumstat_beta, sumstat_se, is_overlap = FALSE, r_mat = NA, kernel_matrix = NA,
                                    gamma, tau_mu, tau_delta,
                                    ash_pi = rep(1/7, 7), sigma_c = c(0, 10^seq(-10, 0, by = 2)),
                                    tau_mu_log = FALSE, tau_delta_log = FALSE,
                                    select_method = c("single_exposure", "minp", "fisher", "none"), single_exp_pop = 1, select_pthresh = 0.05,
                                    random_seed = 2026, mc_iter = 250) {

  #Initializing parameters
  if (tau_mu_log) {
    tau_mu <- exp(tau_mu)
  }

  if (tau_delta_log) {
    tau_delta <- exp(tau_delta)
  }
  K <- length(sumstat_beta) - 1

  #Specifying a kernel matrix if none is given
  if (identical(kernel_matrix, NA)) {
    kernel_matrix <- diag(K + 1)
  }

  if (dim(kernel_matrix)[1] != K + 1 & dim(kernel_matrix)[2] != K + 1) {
    stop("Kernel matrix dimensions specified incorrectly!")
  }

  #Basic checks for compatibility of the ash model of horizontal pleiotropy
  if (length(ash_pi) != length(sigma_c)) {
    stop("incompatible lengths for ash pi and sigma_c components!")
  }

  if (abs(1 - sum(ash_pi)) > 1e-5) {
    stop("ash pi components do not add up!")
  }

  #Setting up the base covariance matrix of the conditional log-likelihood
  if (is_overlap == FALSE & identical(r_mat, NA)) { #the default is a standard diagonal matrix
    r_mat <- diag(K+1)
  }

  m <- c(gamma, rep(1, K))
  M <- diag(m)

  sigma_Ej <- r_mat * outer(sumstat_se, sumstat_se)
  sigma_j <- sigma_Ej + tau_mu * outer(m, m) + tau_delta * crossprod(crossprod(kernel_matrix, M), M)

  #Calculating the unconditional mixture log-likelihood
  C <- length(sigma_c)
  sigma_jc <- list()

  for (comp in 1:C) {
    sigma_jc[[comp]] <- sigma_j
    sigma_jc[[comp]][1,1] <- sigma_jc[[comp]][1,1] + sigma_c[comp]
  }

  sigma_jc_det <- lapply(sigma_jc, function(X){det(X)})
  sigma_jc_inv <- lapply(sigma_jc, function(X){solve(X)})

  mixcomp <- vector()
  for (comp in 1:C) {
    mixcomp[comp] <- log(ash_pi[comp]) - 1/2 * (log(sigma_jc_det[[comp]]) + crossprod(crossprod(sigma_jc_inv[[comp]], sumstat_beta), sumstat_beta) + (K+1)*log(2*pi))
  }
  unconditional_loglik <- matrixStats::logSumExp(mixcomp)

  #Calculating (or estimating) the probability of selection
  unt_tau_mu <- ifelse(tau_mu_log == TRUE, exp(tau_mu), tau_mu)
  unt_tau_delta <- ifelse(tau_delta_log == TRUE, exp(tau_delta), tau_delta)


  if (select_method == "single_exposure") { #the probability of a variant being selected based on a single population strength of exposure association
    select_zscore <- abs(stats::qnorm(select_pthresh / 2))
    log_pselect <- log(2 * stats::pnorm(abs(select_zscore) * sumstat_se[1 + single_exp_pop]/sqrt(unt_tau_mu + unt_tau_delta + sumstat_se[1 + single_exp_pop]^2), lower.tail = FALSE))
  }

  if (select_method == "minp") { #the probability of a variant being selected based on the minimum p-value of exposure association across all populations
    log_pselect <- log(numint_prob_kernel(unt_tau_mu, unt_tau_delta, r_mat = r_mat[-1, -1], SE_vector = sumstat_se[-1], select_pthresh, kernel_matrix = kernel_matrix[-1, -1]))
  }

  if (select_method == "fisher") { #the probability of a variant being selected based on Fisher's method of combining p-values
    log_pselect <- log(resample_prob_kernel(mc_iter, random_seed, unt_tau_mu, unt_tau_delta, r_mat = r_mat[-1, -1], SE_vector = sumstat_se[-1], select_pthresh, kernel_matrix = kernel_matrix[-1, -1]))
  } else {
    log_pselect <- 0
  }

  return(unconditional_loglik - log_pselect)
}

#' Calculating the probability of a variant being selected by Fisher's method of meta-analysis using Monte-Carlo sampling
#'
#' @param resample_number the number of resamples performed to estimate the probability of a variant being selected under Fisher's method.
#' @param resample_seed The seed used for resampling
#' @param tau_mu The variance of the normal distribution underlying mu.
#' @param tau_delta The variance of the normal distributions underlying delta.
#' @param r_mat the matrix from r_mat in metamrash_loglik.single, but for only exposure summary statistics (i.e. first row and column removed)
#' @param SE_vector the vector from sumstat_se in metamrash_loglik.single, but for only exposure summary statistics (i.e. first entry removed)
#' @param select_pthresh The p-value threshold for variant selection
#' @param kernel_matrix the kernel matrix from kernel_matrix in metamrash_loglik.single, but only for the K exposure studies (i.e. first row and column removed)
#'
#' @returns The probability that a variant was selected under Fisher's method for variant selection
#' @export
resample_prob_kernel <- function(resample_number, resample_seed, tau_mu, tau_delta, r_mat, SE_vector, select_pthresh, kernel_matrix) {
  pops <- length(SE_vector)

  if (!is.na(resample_seed)) {
    set.seed(resample_seed)
  }

  #Empirically assessing the probability of selection

  # The covariance matrix stays the same for all loops
  sigma_Ej <- r_mat * outer(SE_vector, SE_vector) + tau_mu * matrix(1, nrow = pops, ncol = pops) + tau_delta * kernel_matrix

  #Critical value for selection based on Fisher's method
  fisher.crit <- stats::qchisq(select_pthresh, df = 2*pops, lower.tail = FALSE)

  #Generating exposure association betas based on their underlying marginal distribution
  sumstat_matrix <- mvtnorm::rmvnorm(n = resample_number, mean = rep(0, pops), sigma = sigma_Ej)
  zscore_matrix <- apply(sumstat_matrix, MARGIN = 1, FUN = function(X){X/SE_vector})
  logpval_matrix <- log(stats::pnorm(-abs(zscore_matrix)) * 2)
  select_indicator <- -2 * colSums(logpval_matrix) > fisher.crit

  select_prob <- ifelse(sum(select_indicator) > 0, sum(select_indicator)/resample_number, 0.5/resample_number)
  return(select_prob)
}

#' Probability of variant selection under the minimum p-value strategy, using numerical integration
#'
#' @param tau_mu The variance of the normal distribution underlying mu.
#' @param tau_delta The variance of the normal distributions underlying delta.
#' @param r_mat the matrix from r_mat in metamrash_loglik.single, but for only exposure summary statistics (i.e. first row and column removed)
#' @param SE_vector the vector from sumstat_se in metamrash_loglik.single, but for only exposure summary statistics (i.e. first entry removed)
#' @param select_pthresh The p-value threshold for variant selection
#' @param kernel_matrix the kernel matrix from kernel_matrix in metamrash_loglik.single, but only for the K exposure studies (i.e. first row and column removed)
#'
#' @returns The probability that a variant was selected under the minimum p-value strategy for variant selection
#' @export
numint_prob_kernel <- function(tau_mu, tau_delta, r_mat, SE_vector, select_pthresh, kernel_matrix) {
  select_zscore <- abs(stats::qnorm(select_pthresh/2))
  int_lb <- -select_zscore*SE_vector
  int_ub <- select_zscore*SE_vector

  pops <- length(SE_vector)

  m <- rep(1, pops)
  M <- diag(m)

  sigma_Ej <- r_mat * outer(SE_vector, SE_vector)
  #sigma_G <- diag(tau_delta, pops) + matrix(tau_mu, nrow = pops, ncol = pops)
  sigma_j <- sigma_Ej + tau_mu * outer(m, m) + tau_delta * crossprod(crossprod(kernel_matrix, M), M)

  #the integral function is the PDF of MVN
  p_cube <- mvtnorm::pmvnorm(lower = int_lb, upper = int_ub, mean = rep(0, pops), sigma = sigma_j, keepAttr = FALSE)

  return(1 - p_cube)
}

#' Log-likelihood function incorporating variant selection and ASH pleiotropy (multi-variant)
#'
#' A conditional log-likelihood function that models pleiotropy using the adaptive shrinkage (ASH) prior distribution proposed by Stephens et al (2015), accounting for variants being pre-selected using either a single population's exposure association p-value, minimum exposure association p-value across all populations, or Fisher's method. The ASH distribution is a mixture of normals with mean zero, and fixed variance components.
#'
#' @param sumstat_beta_list a list of vectors of GWAS effect size estimates with length K+1: first, for the outcome in the target population, then for exposures across K studies. Each vector in the list represents summary statistics for one of J total variants.
#' @param sumstat_se_list a list of vectors of the standard errors for GWAS effect size estimates in sumstat_beta_list
#' @param is_overlap a boolean describing whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat_list a list of J matrices of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap
#' @param kernel_matrix a K+1 x K+1 correlation matrix that captures correlations in exposure effect size heterogeneity across all populations (including the target population).
#' @param gamma The exposure-outcome causal effect
#' @param tau_mu The variance of the normal distribution underlying mu. May be log-transformed.
#' @param tau_delta The variance of the normal distributions underlying delta. May be log-transformed.
#' @param ash_pi A vector of length c representing the mixture probabilities of each of the components of the ASH distribution we use. By default, we assume c = 7 and that probabilities are equally distributed across all mixture variance parameters (1/7).
#' @param sigma_c A vector of length c representing the variance of each of the components of the ASH distribution we used. By default, we specify c = 7 with components 0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, and 1.
#' @param tau_mu_log Whether tau_mu is log-transformed.
#' @param tau_delta_log Whether tau_delta is log-transformed.
#' @param select_method The method by which variants were selected: "single_exposure" = a single population's exposure association p-value, "minp" = minimum exposure association p-value across all populations, "fisher" = Fisher's method.
#' @param single_exp_pop If select_method = "single_exposure", what exposure study it applies to. By default, we assume it is the first.
#' @param select_pthresh The p-value threshold for variant selection
#' @param random_seed For Fisher's method of selection, what random seed to use for sampling of observations
#' @param mc_iter For Fisher's method of selection, how many random samples to take
#' @param standardized A boolean indicating whether summary statistics are standardized; i.e. whether they have the same set of standard errors. If they do, we only need to calculate the probability of selection once across all J variants.
#'
#' @returns A log-likelihood for the summary statistics across multiple (selected) variants given parameters gamma, mu, tau, ash_pi
#' @export
metamrash_loglik.full <- function(sumstat_beta_list, sumstat_se_list, is_overlap = FALSE, r_mat_list = NA, kernel_matrix = NA,
                                  gamma, tau_mu, tau_delta,
                                  ash_pi = rep(1/7, 7), sigma_c = c(0, 10^seq(-10, 0, by = 2)),
                                  tau_mu_log = FALSE, tau_delta_log = FALSE,
                                  select_method = c("single_exposure", "minp", "fisher", "none"), single_exp_pop = 1, select_pthresh = 0.05,
                                  random_seed = 2026, mc_iter = 250, standardized = FALSE) {

  #Basic checks for compatibility of arguments
  if (length(ash_pi) != length(sigma_c)) {
    stop("incompatible lengths for ash pi and sigma_c components!")
  }

  if (abs(1 - sum(ash_pi)) > 1e-5) {
    stop("ash pi components do not add up!")
  }

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

  #If variants are standardized, we only need to calculate the probability of selection once per log-likelihood calculation
  if (standardized) {
    for (i in 1:n) {
      loglik_variants[i] <- metamrash_loglik.single(sumstat_beta = sumstat_beta_list[[i]],
                                                    sumstat_se = sumstat_se_list[[i]],
                                                    is_overlap = TRUE,
                                                    r_mat = r_mat_list[[i]],
                                                    kernel_matrix = kernel_matrix,
                                                    gamma = gamma, tau_mu = tau_mu, tau_delta = tau_delta,
                                                    ash_pi = ash_pi, sigma_c = sigma_c,
                                                    tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                                    select_method = "none", single_exp_pop = single_exp_pop, select_pthresh = select_pthresh,
                                                    random_seed = random_seed, mc_iter = mc_iter)
    }
    unconditional_loglik <- sum(loglik_variants)
    unt_tau_mu <- ifelse(tau_mu_log == TRUE, exp(tau_mu), tau_mu)
    unt_tau_delta <- ifelse(tau_delta_log == TRUE, exp(tau_delta), tau_delta)
    if (select_method == "fisher") {
      log_pselect <- log(resample_prob_kernel(mc_iter, random_seed, unt_tau_mu, unt_tau_delta, r_mat = r_mat_list[[1]][-1, -1], SE_vector = sumstat_se_list[[1]][-1], select_pthresh,
                                              kernel_matrix = kernel_matrix[-1, -1]))
      return(unconditional_loglik - log_pselect*n)
    } else {
      if (select_method == "minp") {
        log_pselect <- log(numint_prob_kernel(unt_tau_mu, unt_tau_delta, r_mat = r_mat_list[[1]][-1, -1], SE_vector = sumstat_se_list[[1]][-1], select_pthresh, kernel_matrix = kernel_matrix[-1, -1]))
        return(unconditional_loglik - log_pselect*n)
      } else {
        if (select_method == "single_exposure") {
          sumstat_se <- sumstat_se_list[[1]]
          select_zscore <- abs(stats::qnorm(select_pthresh / 2))
          log_pselect <- log(2 * stats::pnorm(abs(select_zscore) * sumstat_se[1 + single_exp_pop]/sqrt(unt_tau_mu + unt_tau_delta + sumstat_se[1 + single_exp_pop]^2), lower.tail = FALSE))
          return(unconditional_loglik - log_pselect*n)
        } else {#if not Fisher's method, minimum p, or single-exposure, then default to assuming no selection was performed.
          return(unconditional_loglik)
        }
      }
    }
  } else { #otherwise, we calculate a probability of selection for each variant.
    for (i in 1:n) {
      loglik_variants[i] <- metamrash_loglik.single(sumstat_beta = sumstat_beta_list[[i]],
                                                    sumstat_se = sumstat_se_list[[i]],
                                                    is_overlap = TRUE,
                                                    r_mat = r_mat_list[[i]],
                                                    kernel_matrix = kernel_matrix,
                                                    gamma = gamma, tau_mu = tau_mu, tau_delta = tau_delta,
                                                    ash_pi = ash_pi, sigma_c = sigma_c,
                                                    tau_mu_log = tau_mu_log, tau_delta_log = tau_delta_log,
                                                    select_method = select_method, single_exp_pop = single_exp_pop, select_pthresh = select_pthresh,
                                                    random_seed = random_seed, mc_iter = mc_iter)
    }
    return(sum(loglik_variants))
  }
}

#' Optimizing the conditional log-likelihood under variant selection and ASH pleiotropy for a list of selected variants
#' Due to the large number of parameters that need to be optimized (gamma, tau_mu, tau_delta, and especially ash_pi), direct numerical optimization via optim() is unreliable. Instead, we opt for an EM algorithm.
#'
#' @param sumstat_beta_list a list of vectors of GWAS effect size estimates with length K+1: first, for the outcome in the target population, then for exposures across K studies. Each vector in the list represents summary statistics for one of J total variants.
#' @param sumstat_se_list a list of vectors of the standard errors for GWAS effect size estimates in sumstat_beta_list
#' @param is_overlap a boolean describing whether there is any sample overlap between any of the K+1 studies used for summary statistics. Usually this is true due to overlap between the outcome and exposure GWAS in the target population, but we assume no overlap by default.
#' @param r_mat_list a list of J matrices of estimated (residual) summary statistics correlations due to sample overlap between each set of GWAS summary statistics. Our method does not provide a matrix by default, but specifies a simple diagonal matrix in the event of no sample overlap
#' @param kernel_matrix a K+1 x K+1 correlation matrix that captures correlations in exposure effect size heterogeneity across all populations (including the target population).
#' @param sigma_c A vector of length c representing the variance of each of the components of the ASH distribution we used. By default, we specify c = 7 with components 0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, and 1.
#' @param closest_exposure_pop A single number (from 1 to K) indicating which exposure study is sampled from a population best matches the target population. Used for variant selection.
#' @param select_flag A boolean that indicates whether the EM algorithm should account for selection for not. Ignoring selection allows for a closed-form (but possibly inaccurate) algorithm, while accounting for selection needs to be done numerically.
#' @param select_method The method by which variants were selected: "single_exposure" = a single population's exposure association p-value, "minp" = minimum exposure association p-value across all populations, "fisher" = Fisher's method.
#' @param matching_exp_pop How many exposure studies are sampled from the same population as the target. The EM algorithm changes slightly depending on whether we have population-matched exposures
#' @param single_exp_pop If select_method = "single_exposure", what exposure study it applies to. By default, we assume it is the first.
#' @param select_pthresh The p-value threshold for variant selection
#' @param random_seed For Fisher's method of selection, what random seed to use for sampling of observations
#' @param mc_iter For Fisher's method of selection, how many random samples to take
#' @param fixed_gamma A single number indicating the value that gamma should be fixed at through the EM algorithm. Set to NA by default to indicate that it should not be fixed.
#' @param fixed_tau_mu A single number indicating the value that tau_mu should be fixed at through the EM algorithm. Set to NA by default to indicate that it should not be fixed.
#' @param fixed_tau_delta A single number indicating the value that tau_delta should be fixed at through the EM algorithm. Set to NA by default to indicate that it should not be fixed.
#' @param fixed_ash_pi A vector of length c indicating the values that ash_pi should be fixed at. Set to NA by default to indicate that it should not be fixed.
#' @param em_iter The maximum number of iterations to perform for the EM algorithm.
#' @param em_tol The tolerance threshold for stopping the EM algorithm early. If the Euclidean distance between successive EM iterations follows below this threshold, the algorithm ends early.
#' @param standardized A boolean indicating whether summary statistics are standardized; i.e. whether they have the same set of standard errors. If they do, we only need to calculate the probability of selection once across all J variants.
#' @param verbose If true, prints out the current values of each parameter, as well as the observed data log-likelihood every user-defined number of iterations (defined in verbose_iteration).
#' @param ashpi_zero_thresh The threshold at which a mixture parameter in the ash model is automatically set to zero. By default, set to 0.001.
#' @param verbose_iteration The interval at which to print out parameter updates and iterations. By default, set to 50.
#'
#' @returns A list containing the final parameter values estimated by the EM algorithm, as well as successive iterations. params = the final quantities for gamma/tau_mu/tau_delta/ash_pi, final_loglik = the final observed data log-likelihood at params, params_list = a list of values that each parameter took across successive iterations, loglik_vector = the observed data log-likelihood across all iterations.
#' @export
metamrash_em_select <- function(sumstat_beta_list, sumstat_se_list,
                                is_overlap = FALSE, r_mat_list = NA, kernel_matrix = NA,
                                sigma_c = c(10^seq(-10, 0, by = 2)), closest_exposure_pop = 1,
                                select_flag = FALSE, select_method = c("single_exposure", "minp", "fisher", "none"),
                                matching_exp_pop = "none", single_exp_pop = 1, select_pthresh = 0.05,
                                random_seed = 2026, mc_iter = 250,
                                fixed_gamma = NA, fixed_tau_mu = NA, fixed_tau_delta = NA, fixed_ash_pi = NA,
                                em_iter = 100, em_tol = 1e-5, standardized = FALSE, verbose = FALSE,
                                ashpi_zero_thresh = 1e-3, verbose_iteration = 50) {
  #=Basic initial checks=
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

  pops <- length(sumstat_beta_list[[1]]) - 1
  J <- length(sumstat_beta_list)
  #If no residual correlation matrix is given due to a lack of overlap, generate diag(k) matrices by default
  if (is_overlap == FALSE & identical(r_mat_list, NA)) { #the default is a standard diagonal matrix for all matrices
    diag_mat <- list(diag(pops + 1))
    r_mat_list <- rep(diag_mat, J)
  }

  #=Initializing the parameter vector=
  beta_matrix <- matrix(unlist(sumstat_beta_list), nrow = length(sumstat_beta_list), byrow = TRUE)
  se_matrix <- matrix(unlist(sumstat_se_list), nrow = length(sumstat_beta_list), byrow = TRUE)

  #gamma = simple IVW regression of outcome/exposure from the closest population, or fixed at zero
  if (is.na(fixed_gamma)) {
    init_gamma <- unname(stats::lm(beta_matrix[,1] ~ beta_matrix[,closest_exposure_pop + 1] - 1, weights = 1/se_matrix[,closest_exposure_pop + 1]^2)$coefficients)
  } else {
    init_gamma <- fixed_gamma
  }

  #tau_mu and tau_delta = estimated from the covariance matrix of exposure summary statistics
  beta_x_vars <- stats::var(beta_matrix[,-1])
  if (length(beta_x_vars) == 1) {
    init_tau_mu <- beta_x_vars
    init_tau_delta <- 0
  } else {
    if (identical(unname(kernel_matrix[-1, -1]), diag(pops))) { #if the kernel matrix for exposures is just diag(pops), we can't run regression to get our tau_mu and tau_delta estimates
      init_tau_mu <- unname(mean(c(beta_x_vars[upper.tri(beta_x_vars, diag = FALSE)], beta_x_vars[lower.tri(beta_x_vars, diag = FALSE)])))
      if (init_tau_mu <= 0) {init_tau_mu <- 1e-16}
      se_x_squared <- se_matrix[1,-1]^2
      tau_delta_estimates <- unname(diag(beta_x_vars) - init_tau_mu - se_x_squared)
      init_tau_delta <- mean(tau_delta_estimates)
      if (init_tau_delta <= 0) {init_tau_delta <- 1e-16}
    } else {
      #tau_mu and tau_delta are estimated using linear regression of exposure covariance terms against kernel matrix terms
      reg.y <- as.vector(beta_x_vars[upper.tri(beta_x_vars, diag = FALSE)])
      reg.x <- as.vector(kernel_matrix[-1, -1][upper.tri(kernel_matrix[-1, -1], diag = FALSE)])

      covar_vs_kernel <- unname(stats::lm(reg.y ~ reg.x)$coefficients)

      init_tau_mu <- covar_vs_kernel[1]
      init_tau_delta <- covar_vs_kernel[2]

      if (init_tau_mu <= 0) {init_tau_mu <- 1e-16}
      if (init_tau_delta <= 0) {init_tau_delta <- 1e-16}
    }
  }

  if (!is.na(fixed_tau_mu)) {
    init_tau_mu <- fixed_tau_mu
  }
  if (!is.na(fixed_tau_delta)) {
    init_tau_delta <- fixed_tau_delta
  }

  #Initial ash_pi: mixture probabilities
  if (identical(fixed_ash_pi, NA)) {
    C <- length(sigma_c)
    init_ash_pi <- (C:1)/sum(1:C)
  } else {
    C <- length(sigma_c)
    init_ash_pi <- fixed_ash_pi
  }

  starting_loglik <- metamrash_loglik.full(sumstat_beta_list, sumstat_se_list, is_overlap = is_overlap, r_mat_list = r_mat_list, kernel_matrix = kernel_matrix,
                                           gamma = init_gamma, tau_mu = init_tau_mu, tau_delta = init_tau_delta,
                                           ash_pi = init_ash_pi, sigma_c = sigma_c,
                                           tau_mu_log = FALSE, tau_delta_log = FALSE,
                                           select_method = select_method, single_exp_pop = single_exp_pop, select_pthresh = select_pthresh,
                                           random_seed = random_seed, mc_iter = mc_iter, standardized = standardized)

  if (verbose) {
    print(paste("parameter vectors initialized at", paste(round(c(init_gamma, init_tau_mu, init_tau_delta, init_ash_pi), 6), collapse = "/"), "with log-likelihood", starting_loglik))
  }

  current_gamma <- init_gamma
  current_tau_mu <- init_tau_mu
  current_tau_delta <- init_tau_delta
  current_ash_pi <- init_ash_pi

  param_list <- list()
  loglik_vector <- vector()

  #if there are population-matched exposure GWAS, moving them to the front WLOG
  if (identical(matching_exp_pop, "none")) {
    print("No exposure populations match, proceeding without exposure-swapping and kernel matrix pruning")
  } else {
    #some checks about matching_exp_pop
    if (is.numeric(matching_exp_pop)) {
      stop("matching_exp_pop is not numeric!")
    }

    if (any(matching_exp_pop > pops | matching_exp_pop < 1 | matching_exp_pop %% 1 != 0)) {
      stop("invalid index specified in matching_exp_pop!")
    }

    mcount <- length(matching_exp_pop)
    print(paste("Reindexing population-matched exposure GWAS", paste(matching_exp_pop, collapse = " and "), "as the first", mcount, "exposure GWAS:"))

    nonmatching <- (1:pops)[-matching_exp_pop]

    #Swapping entries of sumstat_beta_list and sumstat_se_list
    sumstat_beta_list <- lapply(sumstat_beta_list, function(x){x[c(1, matching_exp_pop + 1, nonmatching + 1)]})
    sumstat_se_list <- lapply(sumstat_se_list, function(x){x[c(1, matching_exp_pop + 1, nonmatching + 1)]})

    #Swapping the kernel matrix
    kernel_matrix <- kernel_matrix[c(1, matching_exp_pop + 1, nonmatching + 1), c(1, matching_exp_pop + 1, nonmatching + 1)]

    #Swapping r_mat_list
    r_mat_list <- lapply(r_mat_list, function(x){x[c(1, matching_exp_pop + 1, nonmatching + 1), c(1, matching_exp_pop + 1, nonmatching + 1)]})

    #Swapping single_exp_pop
    repo_vector <- c(matching_exp_pop, nonmatching) #where all the exposure populations end up after the swap
    single_exp_pop <- which(repo_vector == single_exp_pop)
  }

  #calculating sigma_ej and its inverses, which remains constant between E-M iterations
  sigma_ej <- list()
  sigma_ej_inv <- list()
  for (j in 1:J) {
    sigma_ej[[j]] <- r_mat_list[[j]] * outer(sumstat_se_list[[j]], sumstat_se_list[[j]])
    sigma_ej_inv[[j]] <- solve(sigma_ej[[j]])
  }

  #same for the inverse of the kernel
  #if there are any population-matched exposures, remove them.
  if (identical(matching_exp_pop, "none")) {
    inv_kernel <- solve(kernel_matrix)
  } else {
    print(paste("Due to the presence of population-matched exposures, the kernel matrix as written is non-invertible; removing the first", mcount, "exposures:"))
    first_m_exp <- 2:(mcount+1)
    inv_kernel <- solve(kernel_matrix[-first_m_exp, -first_m_exp])
  }

  if (!is.na(fixed_gamma) & !is.na(fixed_tau_mu) & !is.na(fixed_tau_delta) & !identical(fixed_ash_pi, NA) ) {
    stop("At least one of gamma/tau_mu/tau_delta/ash_pi needs to be set as not fixed for the algorithm to run!")
  }

  for (iter in 1:em_iter) {
    #=E-step: calculate per-variant conditional distributions of c_j=
    pcj <- matrix(0, nrow = J, ncol = C)

    m <- c(current_gamma, rep(1, pops))
    M <- diag(m)

    for (j in 1:J) {
      #Setting up the base covariance matrix
      sigma_Ej <- r_mat_list[[j]] * outer(sumstat_se_list[[j]], sumstat_se_list[[j]])
      sigma_j <- sigma_Ej + current_tau_mu * outer(m, m) + current_tau_delta * crossprod(crossprod(kernel_matrix, M), M)

      #Calculating the non-logged MVN likelihoods for each component
      component_mvn <- vector()
      for (comp in 1:C) {
        sigma_jc <- sigma_j
        sigma_jc[1,1] <- sigma_jc[1,1] + sigma_c[comp]
        component_mvn[comp] <- mvtnorm::dmvnorm(sumstat_beta_list[[j]], sigma = sigma_jc)
      }

      pcj[j,] <- current_ash_pi * component_mvn/sum(current_ash_pi * component_mvn)
    }

    #=E-step: calculate posterior mean and covariance of mu_j, delta_j, alpha_j conditional on c_j=
    Ezj <- list() #a list of J x C vectors
    Covzj <- list() #a list of J x C matrices

    if (identical(matching_exp_pop, "none")) {
      L0 <- matrix(0, nrow = pops + 1, ncol = pops + 3)
      L0[,1] <- 1
      L0[1,1] <- L0[1,2] <- current_gamma
      L0[1, pops + 3] <- 1
      L0[2:(pops + 1), 3:(pops+2)] <- diag(pops)
    } else {
      L0 <- matrix(0, nrow = pops + 1, ncol = pops + 3 - mcount)
      L0[,1] <- 1
      L0[1,1] <- L0[1,2] <- current_gamma
      L0[1, pops + 3 - mcount] <- 1
      L0[2:(mcount+1), 2] <- 1
      L0[(mcount+2):(pops+1), 3:(pops + 2 - mcount)] <- diag(pops - mcount)
    }

    #the design matrix if a component exists with exact zeroes.
    last_col <- ncol(L0)
    L0_0 <- L0[, -last_col]

    Vzc0_inverse_list <- list()

    if (identical(matching_exp_pop, "none")) {
      for (comp in 1:C) {
        if (current_ash_pi[comp] == 0) {
          next
        }
        if (sigma_c[comp] == 0) {
          Vzc0_inverse <- matrix(0, nrow = pops + 2, ncol = pops + 2)
          Vzc0_inverse[1,1] <- 1/current_tau_mu
          Vzc0_inverse[2:(pops + 2),2:(pops + 2)] <- 1/current_tau_delta * inv_kernel
          Vzc0_inverse_list[[comp]] <- Vzc0_inverse
        } else {
          Vzc0_inverse <- matrix(0, nrow = pops + 3, ncol = pops + 3)
          Vzc0_inverse[1,1] <- 1/current_tau_mu
          Vzc0_inverse[2:(pops + 2),2:(pops + 2)] <- 1/current_tau_delta * inv_kernel
          Vzc0_inverse[(pops + 3), (pops + 3)] <- 1/sigma_c[comp]
          Vzc0_inverse_list[[comp]] <- Vzc0_inverse
        }
      }
    } else {
      for (comp in 1:C) {
        if (current_ash_pi[comp] == 0) {
          next
        }
        if (sigma_c[comp] == 0) {
          Vzc0_inverse <- matrix(0, nrow = pops + 2 - mcount, ncol = pops + 2 - mcount)
          Vzc0_inverse[1,1] <- 1/current_tau_mu
          Vzc0_inverse[2:(pops + 2 - mcount), 2:(pops + 2 - mcount)] <- 1/current_tau_delta * inv_kernel
          Vzc0_inverse_list[[comp]] <- Vzc0_inverse
        } else {
          Vzc0_inverse <- matrix(0, nrow = pops + 3 - mcount, ncol = pops + 3 - mcount)
          Vzc0_inverse[1,1] <- 1/current_tau_mu
          Vzc0_inverse[2:(pops + 2 - mcount), 2:(pops + 2 - mcount)] <- 1/current_tau_delta * inv_kernel
          Vzc0_inverse[(pops + 3 - mcount), (pops + 3 - mcount)] <- 1/sigma_c[comp]
          Vzc0_inverse_list[[comp]] <- Vzc0_inverse
        }
      }
    }

    for (j in 1:J) {
      Covzj_j <- list()
      Ezj_j <- list()
      for (comp in 1:C) {
        if (current_ash_pi[comp] == 0) {
          next
        }
        if (sigma_c[comp] == 0) {
          Covzj_j[[comp]] <- solve(t(L0_0) %*% sigma_ej_inv[[j]] %*% L0_0 + Vzc0_inverse_list[[comp]])
          Ezj_j[[comp]] <- Covzj_j[[comp]] %*% (t(L0_0) %*% sigma_ej_inv[[j]] %*% sumstat_beta_list[[j]])
        } else {
          Covzj_j[[comp]] <- solve(t(L0) %*% sigma_ej_inv[[j]] %*% L0 + Vzc0_inverse_list[[comp]])
          Ezj_j[[comp]] <- Covzj_j[[comp]] %*% (t(L0) %*% sigma_ej_inv[[j]] %*% sumstat_beta_list[[j]])
        }
      }
      Covzj[[j]] <- Covzj_j
      Ezj[[j]] <- Ezj_j
    }

    #=M-step: updating ash_pi.=
    if (identical(fixed_ash_pi, NA)) {
      current_ash_pi <- colSums(pcj)/J
      #as a slight time-save, we can get rid of components that are clearly shrinking to zero
      current_ash_pi <- ifelse(current_ash_pi < ashpi_zero_thresh, 0, current_ash_pi)
      current_ash_pi <- current_ash_pi/sum(current_ash_pi)
    }


    #=M-step: updating gamma=
    if (is.na(fixed_gamma)) {
      Bjc <- 0
      Ajc <- 0
      for (j in 1:J) {
        for (comp in 1:C) {
          if (current_ash_pi[comp] == 0) {
            next
          }
          s_mat_j <- sigma_ej_inv[[j]]
          mu_vec_jc <- Ezj[[j]][[comp]]
          sigma_mat_jc <- Covzj[[j]][[comp]]
          sumstat_beta <- sumstat_beta_list[[j]]

          #the derivative of L^T %*% \Sigma_{Ej}^{-1} %*% L with respect to gamma, but without the gamma components
          if (identical(matching_exp_pop, "none")) {
            dlsl_nogamma <- matrix(0, nrow = pops + 3, ncol = pops + 3)
            dlsl_nogamma[1,1] <- 2 * sum(s_mat_j[1, 2:(pops + 1)])
            dlsl_nogamma[1,2] <- dlsl_nogamma[2,1] <- sum(s_mat_j[1, 2:(pops + 1)])
            dlsl_nogamma[1, 3:(pops + 2)] <- dlsl_nogamma[2, 3:(pops + 2)] <- dlsl_nogamma[3:(pops + 2), 1] <- dlsl_nogamma[3:(pops + 2), 2] <-
              s_mat_j[1,2:(pops + 1)]
            dlsl_nogamma[1, pops + 3] <- dlsl_nogamma[2, pops + 3] <- dlsl_nogamma[pops + 3, 1] <- dlsl_nogamma[pops + 3, 2] <- s_mat_j[1,1]
          } else {
            dlsl_nogamma <- matrix(0, nrow = pops + 3 - mcount, ncol = pops + 3 - mcount)
            dlsl_nogamma[1,1] <- 2 * sum(s_mat_j[1, 2:(pops + 1)])
            dlsl_nogamma[1,2] <- dlsl_nogamma[2,1] <- sum(s_mat_j[1, 2:(pops + 1)]) + sum(s_mat_j[1, 2:(mcount + 1)])
            dlsl_nogamma[1, 3:(pops+2-mcount)] <- dlsl_nogamma[2, 3:(pops+2-mcount)] <- dlsl_nogamma[3:(pops+2-mcount), 1] <- dlsl_nogamma[3:(pops+2-mcount), 2] <-
              s_mat_j[1,(mcount+2):(pops + 1)]
            dlsl_nogamma[1, pops+3-mcount] <- dlsl_nogamma[2, pops+3-mcount] <- dlsl_nogamma[pops+3-mcount, 1] <- dlsl_nogamma[pops+3-mcount, 2] <- s_mat_j[1,1]
          }

          #if sigma_c = 0 then "horizontal pleiotropy" doesn't exist and we chop off the last row and column of dlsl (which corresponds specifically to the contribution from pleiotropy)
          if (sigma_c[comp] == 0) {dlsl_nogamma <- dlsl_nogamma[-nrow(dlsl_nogamma), -ncol(dlsl_nogamma)]}

          #component 1: -1/2 * (dez_lsl_ez)
          dez_lsl_ez_A <- -s_mat_j[1,1] * (mu_vec_jc[1] + mu_vec_jc[2])^2
          #dez_lsl_ez_B <- -(mu_vec_jc[1] + mu_vec_jc[2]) * (mu_vec_jc[1] * sum(s_mat_j[2:(pops + 1), 1]) + sum(s_mat_j[2:(pops + 1), 1] * mu_vec_jc[3:(pops + 2)]) + s_mat_j[1,1] * mu_vec_jc[pops + 3])
          dez_lsl_ez_B <- -1/2 * crossprod(crossprod(dlsl_nogamma, mu_vec_jc), mu_vec_jc)

          #component 2: dbsl_ez
          dbsl_ez_B <- (s_mat_j[1,1] * sumstat_beta[1] + sum(sumstat_beta[2:(pops + 1)] * s_mat_j[2:(pops + 1), 1])) * (mu_vec_jc[1] + mu_vec_jc[2])

          #component 3: -1/2 * (dtr_lsl_cov)
          dtr_lsl_cov_A <- -(s_mat_j[1,1]) * (sigma_mat_jc[1,1] + 2 * sigma_mat_jc[2,1] + sigma_mat_jc[2,2])
          dtr_lsl_cov_B <- -1/2 * sum(dlsl_nogamma * sigma_mat_jc)

          Ajc <- Ajc + pcj[j, comp] * (dez_lsl_ez_A + dtr_lsl_cov_A)
          Bjc <- Bjc + pcj[j, comp] * (dez_lsl_ez_B + dbsl_ez_B + dtr_lsl_cov_B)
        }
      }
      current_gamma <- Bjc/-Ajc
    }

    #=M-step: updating tau_mu and tau_delta=
    if (is.na(fixed_tau_mu) | is.na(fixed_tau_delta)) {
      numerator_mu <- 0
      numerator_delta <- 0

      if (identical(matching_exp_pop, "none")) {
        mcount <- 0
      }

      for (j in 1:J) {
        for (comp in 1:C) {
          if (current_ash_pi[comp] == 0) {
            next
          }
          mu_vec_jc <- Ezj[[j]][[comp]]
          sigma_mat_jc <- Covzj[[j]][[comp]]
          mu_summand <- mu_vec_jc[1]^2 + sigma_mat_jc[1,1]
          delta_summand <- crossprod(crossprod(inv_kernel, mu_vec_jc[2:(pops+2-mcount)]), mu_vec_jc[2:(pops+2-mcount)]) + sum(inv_kernel * sigma_mat_jc[2:(pops+2-mcount), 2:(pops+2-mcount)])
          numerator_mu <- numerator_mu + mu_summand * pcj[j, comp]
          numerator_delta <- numerator_delta + delta_summand * pcj[j, comp]
        }
      }
    }

    if (select_flag) { #coordinate ascent on log-scale
      lps <- function(tau_mu, tau_delta) { #the sum of log probabilities of selection across all variants
        lpsj <- rep(0, J)
        for (j in 1:J) {
          if (standardized == TRUE & j > 1) {
            lpsj[j] <- lpsj[1]
          } else {
            select_zscore <- abs(stats::qnorm(select_pthresh / 2))
            if (select_method == "single_exposure") {
              sumstat_se <- sumstat_se_list[[j]]
              lpsj[j] <- log(2 * stats::pnorm(abs(select_zscore) * sumstat_se[1 + single_exp_pop]/sqrt(tau_mu + tau_delta + sumstat_se[1 + single_exp_pop]^2), lower.tail = FALSE))
            }
            if (select_method == "minp") {
              lpsj[j] <- log(numint_prob_kernel(tau_mu, tau_delta, r_mat = r_mat_list[[j]][-1, -1], SE_vector = sumstat_se_list[[j]][-1], select_pthresh, kernel_matrix = kernel_matrix[-1, -1]))
            }
            if (select_method == "fisher") {
              lpsj[j] <- log(resample_prob_kernel(mc_iter, random_seed, tau_mu, tau_delta, r_mat = r_mat_list[[j]][-1, -1], SE_vector = sumstat_se_list[[j]][-1], select_pthresh, kernel_matrix = kernel_matrix[-1, -1]))
            }
          }
        }
        return(sum(lpsj))
      }

      Qtau <- function(X) { #X = c(log(tm), log(td))
        tau_mu <- exp(X[1])
        tau_delta <- exp(X[2])
        set.seed(random_seed)

        tm_only <- -1/2 * (numerator_mu/tau_mu + J * log(tau_mu))
        td_only <- -1/2 * (numerator_delta/tau_delta + J * (pops+1-mcount) * log(tau_delta))
        lselect <- lps(tau_mu, tau_delta)

        return((tm_only + td_only - lselect)[1,1])
      }

      #optimizing only tau_mu
      if (is.na(fixed_tau_mu)) {
        Qtau_mu <- function(tm) {
          X <- c(tm, log(current_tau_delta))
          return(-Qtau(X))
        }
        log_ctm <- stats::optimize(Qtau_mu, interval = c(-32, 8))$minimum
        current_tau_mu <- exp(log_ctm)
      }

      #optimizing only tau_delta
      if (is.na(fixed_tau_delta)) {
        Qtau_delta <- function(td) {
          X <- c(log(current_tau_mu), td)
          return(-Qtau(X))
        }
        log_ctd <- stats::optimize(Qtau_delta, interval = c(-32, 8))$minimum
        current_tau_delta <- exp(log_ctd)
      }
    } else { #direct optimization
      if (is.na(fixed_tau_mu)) {
        current_tau_mu <- numerator_mu/J
      }

      if (is.na(fixed_tau_delta)) {
        current_tau_delta <- numerator_delta[1,1]/J/(pops+1-mcount)
      }
    }

    current_loglik <- metamrash_loglik.full(sumstat_beta_list, sumstat_se_list, is_overlap = is_overlap, r_mat_list = r_mat_list, kernel_matrix = kernel_matrix,
                                            gamma = current_gamma, tau_mu = current_tau_mu, tau_delta = current_tau_delta,
                                            ash_pi = current_ash_pi, sigma_c = sigma_c,
                                            tau_mu_log = FALSE, tau_delta_log = FALSE,
                                            select_method = select_method, single_exp_pop = single_exp_pop, select_pthresh = select_pthresh,
                                            random_seed = random_seed, mc_iter = mc_iter, standardized = standardized)

    if (verbose) {
      if (iter %% verbose_iteration == 0) {
        print(paste("At iteration", iter, "parameter vectors are now at", paste(round(c(current_gamma, current_tau_mu, current_tau_delta, current_ash_pi), 6), collapse = "/"), "with log-likelihood", current_loglik))
      }
    }

    param_list[[iter]] <- c(current_gamma, current_tau_mu, current_tau_delta, current_ash_pi)
    loglik_vector[iter] <- current_loglik
    if (iter > 1) {
      if (sum((param_list[[iter]] - param_list[[iter - 1]])^2) < em_tol) {print("Tolerance of EM algorithm reached, terminating early!"); break}
    }
  }

  output_list <- list(params = c(gamma = current_gamma, tau_mu = current_tau_mu, tau_delta = current_tau_delta, ash_pi = current_ash_pi), final_loglik = current_loglik, param_list = param_list, loglik_vector = loglik_vector)
  return(output_list)
}
