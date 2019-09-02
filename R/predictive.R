# ============================== predict.ithresh ==============================

#' Predictive inference for the largest value observed in N years.
#'
#' \code{predict} method for class \code{"ithresh"}.  Predictive inferences can
#' either be based on a \emph{single training threshold} or using a weighted
#' average of inferences over \emph{multiple training thresholds}.  A single
#' threshold may chosen to be the best performing threshold, as judged by the
#' measure of predictive performance calculated by \code{\link{ithresh}} or
#' chosen by the user.  The weights used in the latter case are based on the
#' measures of predictive performance and prior probabilities assigned to the
#' training thresholds.  By default all thresholds are given the same
#' prior probability but the user can specify their own prior.
#'
#' @param object An object of class \code{"ithresh"}, a result of a call to
#'   \code{\link{ithresh}} or \code{\link{choose_lambda}}, or an object of
#'   class "bcthresh", a result of a call to \code{\link{bcthresh}}.
#' @param npy A numeric scalar. The mean number of observations per year
#'   of data, after excluding any missing values, i.e. the number of
#'   non-missing observations divided by total number of years of non-missing
#'   data.
#' @param n_years A numeric vector. Value(s) of N.  If \code{which_u = "all"}
#'   then \code{n_years} must have length one.
#' @param which_u Either a character scalar or a numeric scalar.
#'   If \code{which_u} is a character scalar it must be either "best" or "all".
#'
#'   If \code{which_u = "best"} then the single threshold achieving the largest
#'   measure of predictive performance in \code{object$pred_perf}, based
#'   on the validation threshold selected using \code{which_v}, is used to
#'   perform prediction.  See \code{\link{summary.ithresh}} to print the
#'   best thresholds for each validation threshold.
#'
#'   If \code{which_u = "all"} then \emph{all} the thresholds are used to
#'   perform prediction.  The inferences from each threshold are weighted
#'   according to the \emph{posterior threshold weights} given in
#'   equation (15) of
#'   \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#'   based on the prior probabilities of thresholds in \code{u_prior}
#'   and column \code{which_v} of the measures of predictive performance in
#'   \code{object$pred_perf}.
#'
#'   Otherwise, \code{which_u} is a numeric scalar that indicates which
#'   element of \code{object$u_vec} the user wishes to select as a single
#'   threshold on which to base prediction, that is, \code{which_u} must
#'   be an integer in {1, ..., \code{length(object$u_vec)}}.
#' @param which_v A numeric scalar. Indicates which element of
#'   \code{object$v_vec} is used in selecting a single threshold
#'   (if \code{which_u = "best"}) or weighting the inferences from
#'   all thresholds (if \code{which_u = "all"}).
#'   Note: the default, \code{which_v = 1} gives the \emph{lowest} of the
#'   validation thresholds in \code{object$v_vec}.
#' @param u_prior  A numeric vector.  Prior probabilities for the training
#'   thresholds in \code{object$u_vec}.  Only used if \code{which_u = "all"}.
#'
#'   Only the first
#'   \code{length(object$u_vec) - length(object$v_vec) + which_v}
#'   elements of \code{u_prior} are used.  This is because only training
#'   thresholds up to and including \code{object$v_vec[which_v]} are relevant.
#'   \code{u_prior} must have length \code{length(object$u_vec)} or
#'   \code{length(object$u_vec) - length(object$v_vec) + which_v}.
#'
#'   If \code{u_prior} is not supplied then all (relevant) training thresholds
#'   are given equal prior probability.
#'   \code{u_prior} is normalized to have sum equal to 1 inside
#'   \code{predict.ithresh}.
#' @param type A character vector.
#'   Passed to \code{\link[revdbayes]{predict.evpost}}.
#'   Indicates which type of inference is required:
#' \itemize{
#'   \item "p" for the predictive distribution function,
#'   \item "d" for the predictive density function,
#'   \item "q" for the predictive quantile function,
#'   \item "i" for predictive intervals (see \code{...} to set \code{level}),
#'   \item "r" for random generation from the predictive distribution.
#' }
#'   If \code{which_u = "all"} then only \code{type = "p"} or \code{type = "d"}
#'   are allowed.
#' @param hpd A logical scalar.  The argument \code{hpd} of
#'   \code{\link[revdbayes]{predict.evpost}}. Only relevant if
#'   \code{type = "i"}.
#' @param x A numeric vector.  The argument \code{x} of
#'   \code{\link[revdbayes]{predict.evpost}}.  In the current context this
#'   must be a vector (not a matrix, as suggested by the documentation of
#'   \code{\link[revdbayes]{predict.evpost}}).  If \code{x} is not supplied
#'   then a default value is set within
#'   \code{\link[revdbayes]{predict.evpost}}.
#' @param ... Additional arguments to be passed to
#'   \code{\link[revdbayes]{predict.evpost}}.  In particular:
#'   \code{level}, which can be used to set (multiple) levels
#'   for predictive intervals when \code{type = "i"};
#'   \code{lower_tail} (relevant when \code{type = "p"} or \code{"q"}) and
#'   \code{log} (relevant when \code{type = "d"}).
#' @details The function \code{\link[revdbayes]{predict.evpost}} is used to
#'   perform predictive based on the binomial-GP posterior sample generated
#'   using a given training threshold.  For mathematical details of the
#'   single threshold and multiple threshold cases see Sections 2.3 and 3 of
#'   \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#'   respectively.
#' @return An list object of class "ithreshpred" with a similar structure to
#'   an object of class "evpred" returned from
#'   \code{\link[revdbayes]{predict.evpost}} is returned \emph{invisibly}.
#'   In addition, the object contains
#'   \code{u_vec = object$u_vec} and \code{v_vec = object$v_vec},
#'   \code{which_v} and the index \code{best_u} in
#'   \code{u_vec = object$u_vec} of the best training threshold based on
#'   the value of \code{which_v}.
#'   It also contains the value of the Box-Cox transformation parameter
#'   \code{lambda}.  This will always be equal to 1 if \code{object} was
#'   returned from \code{ithresh}.
#'
#'   If \code{which_u = "all"} then
#' \itemize{
#'   \item the list also contains the \emph{posterior threshold weights}
#'     in component \code{post_thresh_wts}
#'   \item the component \code{y} is a matrix with \code{length{x}} rows
#'     and 1 + \code{length(object$u_vec) - length(object$v_vec) + which_v}
#'     columns.  Column 1 contains the estimated predictive distribution
#'     function (\code{type = "p"}) or density function (\code{type = "d"})
#'     obtained using a weighted average of the inferences over different
#'     training thresholds.  The other columns contain the estimated
#'     functions for each of the training thresholds in \code{u_vec}.
#' }
#'
#' @examples
#' # Note:
#' #'  In the examples below validation thresholds rather higher than is
#' #   advisable have been used, with far fewer excesses than the minimum of
#' #   50 suggested by Jonathan and Ewans (2013).
#'
#' # Gulf of Mexico significant wave heights, default priors.
#' u_vec_gom <- quantile(gom, probs = seq(0, 0.9, by = 0.05))
#' gom_cv <- ithresh(data = gom, u_vec = u_vec_gom, n_v = 3)
#'
#' # Note: gom_cv$npy contains the correct value of npy (it was set in the
#' #       call to ithresh, via attr(gom, "npy").
#' #       If object$npy doesn't exist then the argument npy must be supplied
#' #       in the call to predict().
#'
#' ### Best training threshold based on the lowest validation threshold
#'
#' # Predictive distribution function
#' best_p <- predict(gom_cv, n_years = c(100, 1000))
#' plot(best_p)
#'
#' # Predictive density function
#' best_d <- predict(gom_cv, type = "d", n_years = c(100, 1000))
#' plot(best_d)
#'
#' # Predictive intervals
#' best_i <- predict(gom_cv, n_years = c(100, 1000), type = "i", hpd = TRUE,
#'                   level = c(95, 99))
#' plot(best_i, which_int = "both")
#'
#' # See which threshold was used
#' summary(gom_cv)
#'
#' ### All thresholds plus weighted average of inferences over all thresholds
#'
#' # Predictive distribution function
#' all_p <- predict(gom_cv, which_u = "all")
#' plot(all_p)
#'
#' # Predictive density function
#' all_d <- predict(gom_cv, which_u = "all", type = "d")
#' plot(all_d)
#' @seealso \code{\link{ithresh}} for threshold selection in the i.i.d. case
#'   based on leave-one-out cross-validation.
#' @seealso \code{\link{plot.ithreshpred}} for the S3 plot method for objects
#'   of class \code{ithreshpred}.
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{https://doi.org/10.1111/rssc.12159}
#' @export
predict.ithresh <- function(object, npy = NULL, n_years = 100,
                            which_u = c("best", "all"), which_v = 1L,
                            u_prior = rep(1, length(object$u_vec)),
                            type = c("p", "d", "q", "i", "r"), hpd = FALSE,
                            x = NULL, ...) {
  if (!inherits(object, "ithresh")) {
    stop("object must be of class ''ithresh''")
  }
  # Check that which_u and n_years are compatible
  if (which_u[1] == "all" & length(n_years) > 1){
    stop("If which = \"all\" then n_years must have length 1")
  }
  # From which function was object returned?
  if (is.null(object$lambda)) {
    fn_object <- "ithresh"
  } else {
    fn_object <- "choose_lambda"
  }
  # If object was returned from choose_lambda() then we use object$lambda
  if (fn_object == "choose_lambda") {
    lambda <- object$lambda
  }
  # Numbers of training and validation thresholds
  n_u <- length(object$u_vec)
  n_v <- length(object$v_vec)
  # Check inputs
  type <- match.arg(type)
  if (is.character(which_u)) {
    which_u <- match.arg(which_u)
  } else if (is.numeric(which_u)) {
    if (!(which_u %in% 1:n_u)) {
      stop("'which_u' must be in 1:length(object$u_vec)")
    }
  } else {
    stop("'which_u' must be ''best'' or ''all'' or a numeric scalar")
  }
  if (which_u == "all" & type != "p" & type != "d") {
    stop("If 'which_u = ''all'' then 'type' must be ''p'' or ''d''")
  }
  if (!is.numeric(n_years)) {
    stop("'n_years' must be numeric")
  }
  if (!is.numeric(which_v) || !(which_v %in% 1:n_v)) {
    stop("'which_v' must be in 1:length(object$v_vec)")
  }
  # Check that npy has been supplied
  if (!is.null(object$npy) & !is.null(npy)) {
    warning(paste("Two values of npy supplied.  The value npy = ", npy,
                  " from the current call has been used.", sep=""))
  }
  if (!is.null(object$npy) & is.null(npy)) {
    npy <- object$npy
  }
  if (is.null(object$npy) & is.null(npy)) {
    stop("npy must be given, here or in call to ithresh.")
  }
  if (!is.numeric(npy)) {
    stop("'npy' must be numeric")
  }
  # Store best use so that we can return it later
  return_best_u <- which.max(object$pred_perf[, which_v])
  # Select the user's option based on which_u -----------
  if (which_u == "best" || is.numeric(which_u)) {
    # Create a list object of class "evpost" for revdbayes::predict.evpost().
    evpost_obj <- list()
    class(evpost_obj) <- "evpost"
    n <- object$n
    evpost_obj$model <- "bingp"
    # Best threshold
    if (which_u == "best") {
      best_u <- which.max(object$pred_perf[, which_v])
    } else {
      best_u <- which_u
    }
    # Extract posterior sample for this threshold
    which_rows <- (1 + (best_u - 1) * n):(best_u * n)
    evpost_obj$bin_sim_vals <- object$sim_vals[which_rows, 1]
    evpost_obj$sim_vals <- object$sim_vals[which_rows, 2:3]
    evpost_obj$thresh <- object$u_vec[best_u]
    for_predict_evpost <- list(object = evpost_obj, n_years = n_years,
                               npy = npy, type = type, hpd = hpd, x = x, ...)
    # predict is revdbayes:::predict.evpost
    ret_obj <- do.call(revdbayes:::predict.evpost, for_predict_evpost)
  }
  if (which_u == "all") {
    # All thresholds
    ptw <- post_thresh_weights(x = object, which_v = which_v,
                               u_prior = u_prior)
    n_t <- length(ptw)
    # Create a list object of class "evpost" for revdbayes::predict.evpost().
    evpost_obj <- list()
    class(evpost_obj) <- "evpost"
    n <- object$n
    evpost_obj$model <- "bingp"
    # We need to call predict.evpost() with the same x for each threshold
    # Base this on the default x for the 'best' threshold.
    best_u <- which.max(object$pred_perf[, which_v])
    # Extract posterior sample for this threshold
    which_rows <- (1 + (best_u - 1) * n):(best_u * n)
    evpost_obj$bin_sim_vals <- object$sim_vals[which_rows, 1]
    evpost_obj$sim_vals <- object$sim_vals[which_rows, 2:3]
    evpost_obj$thresh <- object$u_vec[best_u]
    for_predict_evpost <- list(object = evpost_obj, n_years = n_years,
                               npy = npy, type = type, hpd = hpd, x = x, ...)
    # predict is revdbayes:::predict.evpost
    ret_obj <- do.call(revdbayes:::predict.evpost, for_predict_evpost)
    x_vals <- ret_obj$x
    others <- (1:n_t)[-best_u]
    # Create matrix: y_val
    y_val <- matrix(NA, ncol = 1 + n_t, nrow = length(x_vals))
    y_val[, 1 + best_u] <- ret_obj$y
    for (i in others) {
      # Extract posterior sample for this threshold
      which_rows <- (1 + (i - 1) * n):(i * n)
      evpost_obj$bin_sim_vals <- object$sim_vals[which_rows, 1]
      evpost_obj$sim_vals <- object$sim_vals[which_rows, 2:3]
      evpost_obj$thresh <- object$u_vec[i]
      for_predict_evpost <- list(object = evpost_obj, n_years = n_years,
                                 npy = npy, type = type, hpd = hpd, x = x_vals,
                                 ...)
      # predict is revdbayes:::predict.evpost
      temp <- do.call(revdbayes:::predict.evpost, for_predict_evpost)
      y_val[, 1 + i] <- temp$y
    }
    # Row-wise mean weighted by the posterior threshold weights
    y_val[, 1] <- apply(y_val[, -1, drop = FALSE], 1, stats::weighted.mean,
                        w = ptw)
    ret_obj$y <- y_val
    ret_obj$post_thresh_wts <- ptw
  }
  ret_obj$which_u <- which_u
  ret_obj$u_vec <- object$u_vec
  ret_obj$which_v <- which_v
  ret_obj$v_vec <- object$v_vec
  ret_obj$best_u <- return_best_u
  # If object was returned from choose_lambda() then transform
  # back to the original scale
  if (fn_object == "choose_lambda") {
    # Add the value of lambda
    ret_obj$lambda <- lambda
    if (ret_obj$type == "p") {
      ret_obj$x <- inv_bc(ret_obj$x, lambda)
    }
    if (ret_obj$type == "d") {
      # Transform the ordinates
      ret_obj$x <- inv_bc(ret_obj$x, lambda)
      # Transformed density, using Jacobian based on untransformed ordinates
      # Code so this works when which_u = "all" and so ret_obj$y is a matrix
      ret_obj$y <- ret_obj$y * as.vector(ret_obj$x ^ (lambda - 1))
    }
    if (ret_obj$type == "q" || ret_obj$type == "r") {
      ret_obj$y <- inv_bc(ret_obj$y, lambda)
    }
    if (ret_obj$type == "i") {
      temp <- apply(ret_obj$long[, c("lower", "upper"), drop = FALSE], 2,
                   inv_bc, lambda = lambda)
      ret_obj$long[, c("lower", "upper")] <- temp
      if (hpd) {
        temp <- apply(ret_obj$short[, c("lower", "upper"), drop = FALSE], 2,
                      inv_bc, lambda = lambda)
        ret_obj$short[, c("lower", "upper")] <- temp
      }
    }
  }
  class(ret_obj) <- "ithreshpred"
  return(invisible(ret_obj))
}

post_thresh_weights <- function(x, which_v = 1, u_prior = NULL) {
  #
  # Calculate the posterior threshold weights using equation (15) of
  # Northrop, Attalides and Jonathan (2017).
  #
  # Args:
  #   x       : A object of class "ithresh" returned by ithresh.
  #   which_v : A numeric scalar.  Indicates which validation threshold,
  #             i.e. which element of x$v_vec and hence which column of
  #             x$pred_perf containing the measures of predictive performance,
  #             is used.
  #   u_prior : A numeric vector.  Prior probabilities for the thresholds
  #             in u_vec.  If this is NULL then it is set to a vector of
  #             length length(x$u_vec) with each element equal to
  #             1/\code{length(x$u_vec)}.
  # Returns:
  #
  # Use only the validation thresholds in column which_v.
  t_vals <- as.numeric(stats::na.omit(x$pred_perf[, which_v]))
  # Check that u_prior is suitable:
  # a numeric vector of length n_u or n_t with no negative entries
  if (!is.numeric(u_prior)) {
    stop("'u_prior' must be a numeric vector")
  }
  n_t <- length(t_vals)
  n_u <- length(x$u_vec)
  if (length(u_prior) != n_u & length(u_prior) != n_t) {
    mess1 <- "'u_prior' must have length \n length(object$u_vec) or"
    mess2 <- "\n length(object$u_vec) - length(object$v_vec) + which_v"
    stop(mess1, mess2)
  }
  if (any(u_prior < 0)) {
    stop("'u_prior' must no have any negative entries")
  }
  # Shoof t_vals so that it is centred on zero (ignoring any -Inf values)
  t_vals <- t_vals - mean(t_vals[is.finite(t_vals)])
  # Truncate u_prior to have length equal to that of t_vals.
  u_prior <- u_prior[1:length(t_vals)]
  # Normalize the remaining weights.
  u_prior <- u_prior / sum(u_prior)
  # Calculate posterior threshold weights.
  ptw <- exp(t_vals) * u_prior
  ptw <- ptw / sum(ptw)
  return(ptw)
}

# ============================== predict.bcthresh =============================

#' Predictive inference for the largest value observed in N years.
#'
#' \code{predict} method for class \code{"bcthresh"}.  Calls
#' \code{\link{predict.ithresh}} for each value of \eqn{\lambda} in
#' \code{lambda}.
#' @param object An object of class \code{"bcthresh"}, a result of a call to
#'   \code{\link{bcthresh}}.
#' @param lambda A numeric vector containing values of the Box-Cox
#'   transformation parameter \eqn{\lambda} to use when performing predictive
#'   inference.  These must be contained in \code{object$lambda}.
#' @param ... Additional arguments to be passed to
#'   \code{\link{predict.ithresh}}.
#' @details If \code{which_u = "all"} then \code{n_years} must have length 1.
#' @examples
#' # Set a prior: flat for GP parameters, Haldane for P(exceedance)
#' prior_args <- list(prior = "flatflat", bin_prior = "haldane",
#'                    h_prior = list(min_xi = -Inf))
#'
#' ## Gulf of Mexico significant wave heights ------------------
#'
#' gprobs <- seq(0.1, 0.8, 0.1)
#' glambda <- seq(1, 3, 0.5)
#' gom_args <- list(data = gom, probs = gprobs, lambda = glambda)
#' gom_lambda <- do.call(bcthresh, c(gom_args, prior_args))
#' dgom <- predict(gom_lambda, type = "d")
#' @export
predict.bcthresh <- function(object, lambda = object$lambda, ...) {
  if (!inherits(object, "bcthresh")) {
    stop("object must be of class ''bcthresh''")
  }
  # Extract arguments that the user want to pass to predict.ithresh
  user_args <- list(...)
  if (is.null(user_args$type)){
    user_args$type <- "p"
  }
  # Check that which_u and n_years are compatible
  if (!is.null(user_args$which_u) & !is.null(user_args$n_years)) {
    if (user_args$which_u == "all" & length(user_args$n_years) > 1){
      stop("If which = \"all\" then n_years must have length 1")
    }
  }
  # Create an object of class "ithresh" for lambda = lambda[1]
  temp <- choose_lambda(object, lambda = lambda[1])
  # If type = "p" or "d" then we need to set values at which to evaluate
  # the cdf or pdf that are equally-spaced on the original (lambda = 1) scale
  # To do this we first make a call to predict.ithresh using x_dum = 2.
  # This estimates the 0.1% and 99% quantiles on the original scale, using
  # we set values to be equally-spaced between these quantiles on this scale
  # and transform them to the transformed scale.
  # However, if the user has supplied x (in user_args$x) then we don't do this.
  # We also respect user_args$x_num, should this exist.
  if (user_args$type %in% c("p", "d") && is.null(user_args$x)) {
    tempx <- predict(temp, x_num = 2, ...)
    if (is.null(user_args$x_num)) {
      x_num <- 100
    } else {
      x_num <- user_args$x_num
    }
    x_vals <- seq(tempx$x[1, 1], tempx$x[2, 1], len = x_num)
    x_vals <- bc(x_vals, lambda[1])
    ret_obj <- predict(temp, x = x_vals, ...)
  } else {
    ret_obj <- predict(temp, ...)
  }
  # Create arrays in which to store the relevant results
  # The 3rd dimension is lambda
  if (user_args$type %in% c("p", "d", "q")) {
    ret_obj$x <- array(ret_obj$x, dim = c(nrow(ret_obj$x), ncol(ret_obj$x),
                                          length(lambda)))
    ret_obj$y <- array(ret_obj$y, dim = c(nrow(ret_obj$y), ncol(ret_obj$y),
                                          length(lambda)))
  } else if (user_args$type == "i") {
    ret_obj$long <- array(ret_obj$long, dim = c(nrow(ret_obj$long),
                                                ncol(ret_obj$long),
                                                length(lambda)))
    ret_obj$short <- array(ret_obj$short, dim = c(nrow(ret_obj$short),
                                                  ncol(ret_obj$short),
                                                  length(lambda)))
  } else if (user_args$type == "r") {
    ret_obj$y <- array(ret_obj$y, dim = c(nrow(ret_obj$y), ncol(ret_obj$y),
                                          length(lambda)))
  }
  # If there are other lambdas then call predict.ithresh() again
  if (length(lambda) > 1) {
    for (i in 2:length(lambda)) {
      temp <- choose_lambda(object, lambda = lambda[i])
      if (user_args$type %in% c("p", "d") && is.null(user_args$x)) {
        tempx <- predict(temp, x_num = 2, ...)
        if (is.null(user_args$x_num)) {
          x_num <- 100
        } else {
          x_num <- user_args$x_num
        }
        x_vals <- seq(tempx$x[1, 1], tempx$x[2, 1], len = x_num)
        x_vals <- bc(x_vals, lambda[i])
        temp <- predict(temp, x = x_vals, ...)
      } else {
        temp <- predict(temp, ...)
      }
      if (user_args$type %in% c("p", "d", "q")) {
        ret_obj$x[, , i] <- temp$x
        ret_obj$y[, , i] <- temp$y
      } else if (user_args$type == "i") {
        ret_obj$long[, , i] <- temp$long
        ret_obj$short[, , i] <- temp$short
      } else if (user_args$type == "r") {
        ret_obj$y[, , i] <- temp$y
      }
      # Add the best threshold
      ret_obj$best_u <- c(ret_obj$best_u, temp$best_u)
    }
  }
  ret_obj$lambda <- lambda
  class(ret_obj) <- "bcthreshpred"
  return(invisible(ret_obj))
}
