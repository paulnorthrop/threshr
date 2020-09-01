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
