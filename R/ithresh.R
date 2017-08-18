# =========================== ithresh ===========================
#
#' Threshold selection in the i.i.d. case (peaks over threshold)
#'
#' Produces a diagnostic plot to assist in the selection of an extreme value
#' threshold in the case where the data can be treated as independent and
#' identically distributed (i.i.d.) observations.  For example, it could be
#' that these observations are the cluster maxima resulting from the
#' declustering of time series data.  The plot is based on the predictive
#' ability of the models fitted using each of a user-supplied set of
#' thresholds, assessed using leave-one-out cross-validation.
#' These models are based on a Generalized Pareto (GP) distribution for
#' threshold excesses and a binomial model for the probability of threshold
#' exceedance.  See
#' \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#' for details.
#'
#' @param data  A numeric vector of observations.
#' @param u_vec A numeric vector. A vector of \emph{training} thresholds
#'   at which inferences are made from the GP model.  Any duplicated values
#'   will be removed. These could be set at sample quantiles of \code{data}
#'   using \code{\link[stats]{quantile}}.
#' @param n_v A numeric scalar.
#'   The \code{n_v} largest values in \code{u_vec} will be used as
#'   validation thresholds to quantify the predictive performance of the
#'   GP models fitted at the thresholds in \code{u_vec}.
#' @param use_rcpp A logical scalar.  If TRUE (the default) use the
#'   revdbayes function \code{\link[revdbayes]{rpost_rcpp}} for posterior
#'   simulation.  Otherwise, we use \code{\link[revdbayes]{rpost}}.
#' @param ... Further (optional) arguments to be passed to the
#'   \code{\link[revdbayes]{revdbayes}} function
#'   \code{\link[revdbayes]{rpost_rcpp}} (or \code{\link[revdbayes]{rpost}}).
#'   In particular:
#' \itemize{
#'   \item {\code{n}} {The size of the posterior sample used to perform
#'     predictive inference.  Default: \code{n = 1000}.}
#'   \item {\code{prior}} {A prior for the GP parameters, set using
#'     \code{\link[revdbayes]{set_prior}}.  Default: \code{prior = "flat"}
#'     with \code{min_xi = -1}.}
#'   \item {\code{h_prior}} {A list of further arguments (hyperparameters)
#'     for the GP prior specified in \code{prior}.}
#'   \item {\code{bin_prior}} {A prior for the threshold exceedance
#'     probability \eqn{p}, set using \code{\link[revdbayes]{set_bin_prior}}.
#'     Default: \code{prior = "jeffreys"}.}
#'   \item {\code{h_bin_prior}} {A list of further arguments (hyperparameters)
#'     for the binomial prior specified in \code{bin_prior}.}
#' }
#' @details For a given threshold in \code{u_vec}:
#' \itemize{
#'   \item {the number of values in \code{data} that exceed the threshold,
#'     and the amounts (the \emph{threshold excesses}) by which these value
#'     exceed the threshold are calculated;}
#'   \item {\code{\link[revdbayes]{rpost_rcpp}}
#'     (or \code{\link[revdbayes]{rpost}}) is used to sample from the posterior
#'     distributions of the parameters of a GP model for the threshold
#'     excesses and a binomial model for the probability of threshold
#'     exceedance;}
#'   \item {the ability of this binomial-GP model to predict data
#'     thresholded at the threshold(s) specified by \code{n_v} is
#'     assessed using leave-one-out cross-validation (the measure of
#'     this is given in equation (7) of
#'     \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}.}
#' }
#'   See \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#'   and the introductory threshr vignette for further details and examples.
#' @return An object (list) of class \code{"thresh"}, containing the components
#'   \itemize{
#'     \item{\code{pred_perf}:} A numeric matrix with \code{length(u_vec)}
#'     rows and \code{n_v} columns.  Each column contains the values of
#'     the measure of predictive performance
#'     \item{\code{u_vec}:} The argument \code{u_vec} to \code{ithresh}.
#'     \item{\code{v_vec}:} A numeric vector.  The validation thresholds
#'       implied by the argument \code{n_v} to \code{ithresh}.
#'     \item{\code{u_ps}:} A numeric vector. The approximate levels of the
#'       sample quantiles to which the values in \code{u_vec} correspond,
#'       i.e. the approximate percentage of the data the lie at or below
#'       each element in \code{u_vec}.
#'     \item{\code{v_ps}:} A numeric vector.  The values in \code{u_ps}
#'       that correspond to the validation thresholds.
#'     \item{\code{data}:} The argument \code{data} to \code{ithresh}
#'       detailed above, with any missing values removed.
#'   }
#' @seealso \code{\link{plot.thresh}} for the S3 plot method for objects of
#'   class \code{thresh}.
#' @seealso \code{\link[revdbayes]{rpost}} in the
#'   \code{\link[revdbayes]{revdbayes}} package for details of the arguments
#'   that can be passed to
#'   \code{\link[revdbayes]{rpost_rcpp}}/\code{\link[revdbayes]{rpost}}.
#' @seealso \code{\link[revdbayes]{set_prior}} in the
#'   \code{\link[revdbayes]{revdbayes}} package for details of how to set a
#'   prior distribution for GP parameters.
#' @seealso \code{\link[revdbayes]{set_bin_prior}} in the
#'   \code{\link[revdbayes]{revdbayes}} package for details of how to set a
#'   prior distribution for the exceedance probability \eqn{p}.
#' @seealso \code{\link[stats]{quantile}}.
#' @examples
#'
#' # Gulf of Mexico significant wave heights
#' u_vec <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
#' gom_cv <- ithresh(data = gom, u_vec = u_vec, n_v = 4)
#' plot(gom_cv, lwd = 2)
#'
#' # North Sea significant wave heights
#' u_vec <- quantile(ns, probs = seq(0, 0.95, by = 0.05))
#' ns_cv <- ithresh(data = ns, u_vec = u_vec, n_v = 3,
#'   prior = "mdi", h_prior = list(a = 0.6))
#' plot(ns_cv, lwd = 2)
#' @references Northrop, P.J. and Attalides, N. (2016) Posterior propriety in
#'   Bayesian extreme value analyses using reference priors
#'   \emph{Statistica Sinica}, \strong{26}(2), 721--743
#'   \url{http://dx.doi.org/10.5705/ss.2014.034}.
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{http://dx.doi.org/10.1111/rssc.12159}
#' @export
ithresh <- function(data, u_vec, n_v = 1, use_rcpp = TRUE, ...) {
  # Remove missing values from data
  data <- as.numeric(na.omit(data))
  # Put thresholds in ascending order and remove any repeated values.
  u_vec <- unique(sort(u_vec))
  n_u <- length(u_vec)
  # Check that the highest threshold is lower than the largest observation.
  if (u_vec[n_u] >= max(data)) {
    stop("max(u_vec) must be less than max(data)")
  }
  if (n_v > n_u) {
    n_v <- n_u
    warning("n_v has been set to length(u_vec)")
  }
  v_vec <- u_vec[(n_u - n_v + 1):n_u]
  temp <- cv_fn(data = data, u_vec = u_vec, v_vec = v_vec, n_u = n_u,
                n_v = n_v, use_rcpp = use_rcpp, ...)
  if (is.null(names(u_vec))) {
    temp$u_ps <- round(100 * sapply(u_vec, function(x) mean(data < x)))
  } else {
    temp$u_ps <- as.numeric(substr(names(u_vec), 1, nchar(names(u_vec),
                                                     type = "c") - 1))
  }
  if (is.null(names(v_vec))) {
    temp$v_ps <- round(100 * sapply(v_vec, function(x) mean(data < x)))
  } else {
    temp$v_ps <- as.numeric(substr(names(v_vec), 1, nchar(names(v_vec),
                                                     type = "c") - 1))
  }
  temp$data <- data
  class(temp) <- "thresh"
  return(temp)
}

# =========================== cv_fn ===========================

cv_fn <- function(data, u_vec, v_vec, n_u, n_v, use_rcpp, ...) {
  # Extract arguments for passing to revdbayes function rpost -----------------
  # GP prior.
  cv_control <- list(...)
  if (is.null(cv_control$prior)) {
    cv_control$prior <- "flat"
  }
  if (is.null(cv_control$h_prior$min_xi)) {
    cv_control$h_prior$min_xi <- -1
  }
  for_set_prior <- c(list(prior = cv_control$prior, model = "gp"),
                     cv_control$h_prior)
  gp_prior <- do.call(revdbayes::set_prior, for_set_prior)
  cv_control$prior <- NULL
  cv_control$h_prior <- NULL
  # Binomial prior.
  if (is.null(cv_control$bin_prior)) {
    cv_control$bin_prior <- "jeffreys"
  }
  for_set_bin_prior <- c(list(prior = cv_control$bin_prior),
                         cv_control$h_bin_prior)
  bin_prior <- do.call(revdbayes::set_bin_prior, for_set_bin_prior)
  cv_control$bin_prior <- NULL
  cv_control$h_bin_prior <- NULL
  # Size of samples from posterior distributions.
  if (is.null(cv_control$n)) {
    n <- 1000
  } else {
    n <- cv_control$n
  }
  cv_control$n <- NULL
  # Check for arguments passed in ... that are not formal arguments of
  # rpost_rcpp/rpost or ru_rcpp/ru.
  in_rpost <- names(cv_control) %in% methods::formalArgs(revdbayes::rpost_rcpp)
  in_ru <- names(cv_control) %in% methods::formalArgs(rust::ru_rcpp)
  rogue_names <- !in_rpost & !in_ru
  rogue_args <- cv_control[rogue_names]
  if (length(rogue_args) > 0) {
    warning("The following arguments have been ignored", immediate. = TRUE)
    print(rogue_args)
    cv_control[rogue_names] <- NULL
  }
  # Set up quantities for passing to revdbayes::rpost()
  # Locate the largest observation(s) in the data.
  j_max <- which(data == max(data))
  data_max <- data[j_max]
  n_max <- length(data_max)
  # Remove (one of the) the largest observation(s)
  data_rm <- data[-j_max[1]]
  n_rm <- length(data_rm)
  #
  # A matrix to store the measures of predictive performance.
  pred_perf <- matrix(NA, nrow = n_u, ncol = n_v)
  pred_perf_rcpp <- matrix(NA, nrow = n_u, ncol = n_v)
  # Set the revdbayes posterior simulation function.
  if (use_rcpp) {
    gp_postsim <- revdbayes::rpost_rcpp
  } else {
    gp_postsim <- revdbayes::rpost
  }
  #
  for_post <- c(list(n = n, model = "bingp", prior = gp_prior,
                   bin_prior = bin_prior), cv_control)
  # Loop over the training thresholds -------------------------------------
  for (i in 1:n_u){
    u <- u_vec[i]
    # Simulate from (full) bin-GP posterior.
    temp <- do.call(gp_postsim, c(for_post, list(data = data, thresh = u)))
    # Simulate from the bin-GP posterior after removal of the maximum value.
    temp_rm <- do.call(gp_postsim, c(for_post, list(data = data_rm, thresh = u)))
    # Combine binomial and GP posterior simulated values.
    theta <- cbind(temp$bin_sim_vals, temp$sim_vals)
    theta_rm <- cbind(temp_rm$bin_sim_vals, temp_rm$sim_vals)
    # Carry out leave-one-out Bayesian cross-validation -------------------
    which_v <- which(v_vec >= u)
    v_vals <- v_vec[which_v]           # select the validation thresholds
    pred_perf[i, which_v] <- bloocv(z = data, theta = theta,
                                    theta_rm = theta_rm, u1 = u,
                                    u2_vec = v_vals, z_max = data_max,
                                    z_rm = data_rm, n = n)
  }
  temp <- list(pred_perf = pred_perf, u_vec = u_vec, v_vec = v_vec)
  return(temp)

}

# =========================== new_bloocv ===========================

bloocv <- function(z, theta, theta_rm, u1, u2_vec, z_max, z_rm, n){
  #
  n1 <- sum(z_rm <= u1)      # number of values that do not exceed u1
  z_gt_u1 <- z_rm[z_rm > u1] # data greater than u1 (not including maximum)
  n2 <- length(u2_vec)       # number of validation thresholds considered
  n_max <- length(z_max)     # cardinality of max(z)
  #
  # Posterior samples of (p,sigma,xi) at training threshold u1 ...
  p1 <- theta[, 1]           # p.u from posterior sample
  sigma_1 <- theta[, 2]
  xi <- theta[, 3]
  xi_near_0 <- abs(xi) < 1e-6
  pow1 <- -1 / xi
  pow2 <- -1 + pow1
  mp1 <- 1 / (1 - p1)
  # Probability of exceeding validation threshold u2, for each u2
  p_gt_u2 <- function(x) {
    xx <- (x - u1) / sigma_1
    ifelse(xi_near_0, exp(-xx + xi * xx ^ 2 / 2), (1 + xi * xx) ^ pow1)
  }
  templ <- numeric(n)
  p_2_vec <- p1 * vapply(u2_vec, p_gt_u2, templ)
  # n.sim by (1+length(u2_vec)) matrix
  # bin-GP density at each observation greater than u1
  f_gt_u1 <- function(x) {
    xx <- (x - u1) / sigma_1
    ifelse(xi_near_0, exp(-xx + xi * xx * (xx - 2) / 2), (1 + xi * xx) ^ pow2)
  }
  fgp_r <- vapply(z_gt_u1, f_gt_u1, templ) / sigma_1
  #
  temp <- p1 * fgp_r
  p2 <- p_2_vec
  #
  #----------------------- Function pred_u2 ---------------------------#
  pred_u2 <- function(k) {
    # Contribution for each z that is less than or equal to u1
    t1 <- 1
    if (n1 > 0) {
      t1 <- mean((1 - p2[, k]) * mp1) / mean(mp1)
    }
    # Contributions from zs that are in (u1,u2]
    t2 <- 1
    if (min(z_gt_u1) <= u2_vec[k]) {
      w2 <- which(z_gt_u1 <= u2_vec[k])
      fr2 <- temp[, w2, drop = FALSE]
      mp2 <- 1 / fr2
      t2 <- colMeans((1 - p2[, k]) * mp2) / colMeans(mp2)
    }
    # Contributions from zs that are greater than u2
    t3 <- 1
    if (max(z_gt_u1) > u2_vec[k]) {
      w3 <- which(z_gt_u1 > u2_vec[k])
      fr3 <- temp[, w3, drop = FALSE]
      t3 <- 1 / colMeans(1 / fr3)
    }
    #
    return(n1 * log(t1) + sum(log(t2)) + sum(log(t3)))
  }
  #--------------------- End of function pred_u2 ----------------------#
  #
  t123 <- vapply(1:n2, pred_u2, 0)
  #
  # Contribution from z_max
  p1 <- theta_rm[, 1]
  sigma_1 <- theta_rm[, 2]
  xi <- theta_rm[, 3]
  xi_near_0 <- abs(xi) < 1e-6
  pow2 <- -(1 + 1 / xi)
  # bin-GP density at max(z)
  xx <- (z_max - u1) / sigma_1
  lin_max <- 1 + xi * xx
  t4 <- ifelse(lin_max > 0,
               ifelse(xi_near_0, exp(-xx + xi * xx * (xx - 2) / 2),
                      lin_max ^ pow2),
               0)
  t4 <- mean(p1 * t4 / sigma_1)
  #
  return(t123 + n_max * log(t4))
}

