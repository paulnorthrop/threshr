# ================================= bcthresh ==================================
#
#' Threshold selection in the i.i.d. case, assisted by Box-Cox transformation
#'
#' Uses the function \code{\link{ithresh}} to assist in the selection of
#' an extreme value threshold after transforming the data with a Box-Cox
#' transformation.  The general aim is to seek a Box-Cox transformation
#' parameter \eqn{\lambda} for which we can set a relatively low threshold,
#' that is, a threshold that has a larger number of threshold excesses than
#' for other values of \eqn{\lambda}.  For
#'
#' @inheritParams ithresh
#' @param probs A numeric vector of non-exceedance probabilities in [0, 1).
#'   \code{\link[stats]{quantile}}, specifically
#'   \code{u_vec = stats::quantile(data, probs)}, is used to set a vector
#'   \code{u_vec} of \emph{training} thresholds for use on the untransformed
#'   (\eqn{\lambda = 1}) data scale.
#' @param lambda A numeric vector containing values of the Box-Cox
#'   transformation parameter \eqn{\lambda}.
#' @details See \code{ithresh} and/or
#'   \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#'   for details of the threshold selection algorithm that is applied
#'   for a given value of \eqn{\lambda}.  The measure of predictive
#'   performance is calculated on the scale of the raw (\eqn{\lambda = 1})
#'   data, in order that values resulting from different values of
#'   \eqn{\lambda} can be compared.
#' @return An object (list) of class \code{"ithresh"}, containing the
#'   components
#'   \itemize{
#'     \item{\code{pred_perf}:} A numeric array with dimensions
#'     \code{length(probs)} (the number of training thresholds),
#'     \code{length(n_v)} (the number of validation thresholds)
#'     and \code{length(lambda)} (the number of value of \eqn{lambda}),
#'     containing values of the measures predictive performance.
#'     \code{pred_perf[ , , k]} is an \code{length(probs)} by
#'     \code{length(n_v)} matrix resulting from the analysis that uses
#'     \eqn{\lambda} = \code{lambda[k]}.
#'     Entries corresponding to cases where the training threshold is above
#'     the validation threshold will be \code{NA}.
#'     \item{\code{u_vec}:} The argument \code{u_vec} to \code{ithresh}.
#'     \item{\code{v_vec}:} A numeric vector.  The validation thresholds
#'       implied by the argument \code{n_v} to \code{ithresh}.
#'     \item{\code{u_ps}:} A numeric vector. The approximate levels of the
#'       sample quantiles to which the values in \code{u_vec} correspond,
#'       i.e. the approximate percentage of the data the lie at or below
#'       each element in \code{u_vec}.
#'     \item{\code{v_ps}:} A numeric vector.  The values in \code{u_ps}
#'       that correspond to the validation thresholds.
#'     \item{\code{sim_vals}:} A numeric matrix with 4 columns and
#'       \code{n} x \code{length(u_vec)} rows.  The \eqn{j}th block of
#'       \code{n} rows contains in columns 1-3 the posterior samples of
#'       the threshold exceedance probability, the GP scale
#'       parameter and the GP shape parameter respectively, and in
#'       column 4 the value of \eqn{j}.
#'     \item{\code{n}:} A numeric scalar.  The value of \code{n}.
#'     \item{\code{npy}:} A numeric scalar.  The value of \code{npy}.
#'     \item{\code{data}:} The argument \code{data} to \code{ithresh}
#'       detailed above, with any missing values removed.
#'     \item{\code{use_rcpp}:} A logical scalar indicating whether
#'       \code{\link[revdbayes]{rpost_rcpp}} (\code{use_rcpp = TRUE}) or
#'       \code{\link[revdbayes]{rpost}} (\code{use_rcpp = FALSE})
#'       was used for posterior simulation.
#'     \item{\code{for_post}:} A list containing arguments with which
#'       \code{\link[revdbayes]{rpost_rcpp}}
#'       (or \code{\link[revdbayes]{rpost}}) was called, including
#'       any user-supplied arguments to these functions.
#'   }
#' @seealso \code{\link{ithresh}} for threshold selection in the i.i.d. case.
#' @seealso \code{\link[revdbayes]{rpost}} in the
#'   \code{\link[revdbayes]{revdbayes}} package for details of the arguments
#'   that can be passed to
#'   \code{\link[revdbayes]{rpost_rcpp}}/\code{\link[revdbayes]{rpost}}.
#' @seealso \code{\link[revdbayes]{set_prior}} and
#'   \code{\link[revdbayes]{set_bin_prior}} in the
#'   \code{\link[revdbayes]{revdbayes}} package for details of how to set a
#'   prior distributions for GP parameters and for the exceedance probability
#'   \eqn{p}.
#' @seealso \code{\link[stats]{quantile}}.
#' @examples
#'
#' library(revdbayes)
#' probs <- seq(0.1, 0.9, 0.4)
#' lambda <- seq(-1/2, 1, 0.5)
#' y <- rexp(1000)
#' x <- exp(y)
#' pjn <- bcthresh(data = x, probs = probs, lambda = lambda,
#'                 prior = "flatflat", bin_prior = "haldane")
#'
#' library(revdbayes)
#' probs <- seq(0.1, 0.9, 0.4)
#' lambda <- seq(1, 3, 0.5)
#' gom_lambda <- bcthresh(data = gom, probs = probs, lambda = lambda,
#'                        prior = "flatflat", bin_prior = "haldane")
#'
#' library(revdbayes)
#' probs <- seq(0.1, 0.9, 0.4)
#' lambda <- seq(0, 2, 0.5)
#' ns_lambda <- bcthresh(data = ns, probs = probs, lambda = lambda,
#'                       prior = "flatflat", bin_prior = "haldane", n_v = 2)
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{http://dx.doi.org/10.1111/rssc.12159}
#' @export
bcthresh <- function(data, probs, lambda, ..., n_v = 1, npy = NULL,
                      use_rcpp = TRUE) {
  # Store npy (if it has been supplied)
  if (!is.null(attr(data, "npy"))) {
    return_npy <- attr(data, "npy")
  }
  if (!is.null(npy) & !is.null(attr(data, "npy"))) {
    warning(paste("Two values of npy supplied.  attr(data, ''npy'') = ",
                  return_npy, " will be returned.", sep=""),
            immediate. = TRUE)
  } else if (!is.null(npy)) {
    return_npy <- npy
  } else if (is.null(attr(data, "npy"))) {
    return_npy <- NULL
  }
  # Remove missing values from data
  data <- as.numeric(stats::na.omit(data))
  # Check that probs contains values in [0, 1)
  if (any(probs < 0 || probs >= 1)) {
    stop("''probs'' must contain values in [0, 1)")
  }
  # Set the thresholds on the untransformed data scale
  # Put the non-exceedance probabilites in increasing order, remove repeats
  probs <- unique(sort(probs))
  u_vec <- stats::quantile(data, probs)
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
  # Loop over the values in lambda
  # Save the raw_data
  raw_data <- data
  n_lambda <- length(lambda)
  store_pred_perf <- array(dim = c(n_u, n_v, n_lambda))
  for (i in 1:n_lambda) {
    print(lambda[i])
    # Transform the data and the thresholds
#    bc_data <- bc_vec(raw_data, lambda = lambda[i])
#    bc_u_vec <- bc_vec(u_vec, lambda = lambda[i])
#    bc_v_vec <- bc_vec(v_vec, lambda = lambda[i])
    lngm <- mean(log(raw_data))
    lngm <- 0
    bc_data <- bc_gm(raw_data, lambda = lambda[i], lngm = lngm)
    bc_u_vec <- bc_gm(u_vec, lambda = lambda[i], lngm = lngm)
    bc_v_vec <- bc_gm(v_vec, lambda = lambda[i], lngm = lngm)
#    print(summary(bc_data))
#    print(bc_u_vec)
#    print(bc_v_vec)
#    stop()
    temp <- bccv_fn(data = bc_data, u_vec = bc_u_vec, v_vec = bc_v_vec,
                    n_u = n_u, n_v = n_v, use_rcpp = use_rcpp,
                    raw_data = raw_data, lambda = lambda[i], lngm = lngm, ...)
    store_pred_perf[ , , i] <- temp$pred_perf
  }
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
  if (!is.null(return_npy)) {
    temp$npy <- return_npy
  }
  temp$data <- data
  temp$pred_perf <- store_pred_perf
  class(temp) <- "ithresh"
  return(temp)
}

# =========================== cv_fn ===========================

bccv_fn <- function(data, u_vec, v_vec, n_u, n_v, use_rcpp, raw_data, lambda,
                    lngm, ...) {
  # Extract arguments for passing to revdbayes function rpost -----------------
  # GP prior.
  cv_control <- list(...)
  # If prior is a (user-supplied) R function then set use_rcpp = FALSE
  if (is.function(cv_control$prior)) {
    use_rcpp <- FALSE
  }
  if (is.null(cv_control$prior)) {
    cv_control$prior <- "mdi"
  }
  if (is.null(cv_control$h_prior$min_xi)) {
    cv_control$h_prior$min_xi <- -1
  }
  if (is.character(cv_control$prior)) {
    if (cv_control$prior == "mdi" & is.null(cv_control$h_prior$a)) {
      cv_control$h_prior$a <- 0.6
    }
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
#  # If the "flatflat" prior is used then the posterior mode is equal to the
#  # MLE.  This may cause a false positive convergence warning, because we
#  # start too close to the mode.  To avoid this we calculate the MLE, perturb
#  # it slightly and pass this to revdbayes::rpost() or revdbayes::rpost_rcpp()
#  # as the initial estimate of the posterior mode.
#  print(cv_control)
#  if (gp_prior$prior = "flatflat" && is.null(cv_control$init_ests)) {
#    cv_control$init_ests <- gp_mle(data - u)
#  }
#  print(cv_control)
  for_post <- c(list(n = n, model = "bingp", prior = gp_prior,
                     bin_prior = bin_prior), cv_control)
  # A matrix to store the posterior samples for each threshold
  sim_vals <- matrix(NA, ncol = 4, nrow = n * n_u)
  # Loop over the training thresholds -------------------------------------
  for (i in 1:n_u){
    u <- u_vec[i]
    # Simulation from posterior distributions.
    #
    # If an error occurs (this can sometimes happen if there are few excesses)
    # then try trans = "BC".  This tends to be slower but the Box-Cox
    # transformation towards normality can improve stability of the
    # ratio-of-uniforms boxing optimizations.  If trans = "BC" was used already
    # then try trans = "none".
    #
    # Simulate from (full) bin-GP posterior.
    temp <- try(do.call(gp_postsim, c(for_post, list(data = data,
                                                     thresh = u))),
                silent = TRUE)
    if (inherits(temp, "try-error")) {
      if (is.null(for_post$trans) || for_post$trans == "none") {
        for_post$trans <- "BC"
      } else {
        for_post$trans <- "none"
      }
      temp <- do.call(gp_postsim, c(for_post, list(data = data, thresh = u)))
    }
    # Simulate from the bin-GP posterior after removal of the maximum value
    temp_rm <- try(do.call(gp_postsim, c(for_post, list(data = data_rm,
                                                        thresh = u))),
                   silent = TRUE)
    if (inherits(temp_rm, "try-error")) {
      if (is.null(for_post$trans) || for_post$trans == "none") {
        for_post$trans <- "BC"
      } else {
        for_post$trans <- "none"
      }
      temp_rm <- do.call(gp_postsim, c(for_post, list(data = data_rm,
                                                      thresh = u)))
    }
    # Combine binomial and GP posterior simulated values.
    theta <- cbind(temp$bin_sim_vals, temp$sim_vals)
    theta_rm <- cbind(temp_rm$bin_sim_vals, temp_rm$sim_vals)
    # Carry out leave-one-out Bayesian cross-validation -------------------
    which_v <- which(v_vec >= u)
    v_vals <- v_vec[which_v]           # select the validation thresholds
    pred_perf[i, which_v] <- bcbloocv(z = data, theta = theta,
                                      theta_rm = theta_rm, u1 = u,
                                      u2_vec = v_vals, z_max = data_max,
                                      z_rm = data_rm, n = n,
                                      raw_rm = raw_data[-j_max],
                                      raw_z_max = raw_data[j_max],
                                      lambda = lambda, lngm = lngm)
    # Save the posterior samples values
    which_rows <- (1 + (i - 1) * n):(i * n)
    sim_vals[which_rows, ] <- cbind(theta, i)
  }
  temp <- list(pred_perf = pred_perf, u_vec = u_vec, v_vec = v_vec,
               sim_vals = sim_vals, n = n, for_post = for_post,
               use_rcpp = use_rcpp)
  return(temp)

}

# =========================== new_bloocv ===========================

bcbloocv <- function(z, theta, theta_rm, u1, u2_vec, z_max, z_rm, n,
                     raw_rm, raw_z_max, lambda, lngm){
  gm <- exp(lngm)
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
#      # Use the Box-Cox Jacobian to make all predictions on the raw data scale
#      if (lambda == 0) {
#        print(exp(-z_gt_u1[w3] / gm) * gm)
#        print((raw_rm_gt_u1[w3] / gm) ^ (lambda - 1))
#        stop()
#        t3 <- t3 * (raw_rm_gt_u1[w3] / gm) ^ (lambda - 1)
#        t3 <- t3 * exp(-z_gt_u1[w3] / gm) * gm
#      }
#      if (lambda != 0) {
#        t3 <- t3 * (lambda * gm ^ (lambda - 1) * z_gt_u1[w3] + 1) ^
#          (1 - 1 / lambda) / gm ^ (lambda - 1)
#      }
      t3 <- t3 * (raw_rm_gt_u1[w3] / gm) ^ (lambda - 1)
    }
    #
    return(n1 * log(t1) + sum(log(t2)) + sum(log(t3)))
  }
  #--------------------- End of function pred_u2 ----------------------#
  #
  raw_rm_gt_u1 <- raw_rm[z_rm > u1]
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
#  # Use the Box-Cox Jacobian to make all predictions on the raw data scale
#  if (lambda==0) {
#    t4 <- t4 * exp(-z_max / gm) * gm
#  }
#  if (lambda!=0) {
#    t4 <- t4 * (lambda * gm ^ (lambda - 1) * z_max + 1) ^ (1 - 1 / lambda) /
#      gm ^ (lambda - 1)
#  }
  t4 <- t4 * raw_z_max ^ (lambda - 1)
  #
  return(t123 + n_max * log(t4))
}

