# =================================== rbc =====================================

#' Simulation with (Inverse) Box-Cox Transformation
#'
#' Random generation for a randm variable that has the following property:
#' a Box-Cox transformation of the variable has a specified target
#' distribution.
#' @param n An integer scalar.  The sample size required.
#' @param lambda A numeric scalar.  The Box-Cox transformation parameter.
#' @param sim_fn An R function that simulates from the target probability
#'   distribution.
#' @param ... Further arguments to be passed to \code{sim_fn}.
#' @details A sample of size \code{n} is simulated using the function supplied
#'   in \code{sim_fn}.  Then the inverse of the Box-Cox transformation is
#'   applied to each value in this sample.
#' @return A numeric vector of length \code{n}.
#' @examples
#' x <- rbc(n = 100, lambda = 2)
#' x <- rbc(n = 100, lambda = 2, sim_fn = revdbayes::rgp)
#' x <- rbc(n = 100, lambda = 2, sim_fn = revdbayes::rgp, shape = 0.1)
#' @export
rbc <- function(n = 1, lambda = 1, sim_fn = stats::rexp, ...) {
  return(inv_bc(sim_fn(n, ...), lambda = lambda))
}

# ============================== bc_sim_study =================================

# Change arguments to lists: rbc_args and bcthresh_args?
# Less ambiguity but clunkier.

#' Simulation study
#'
#' Add description
#' @param sims An integer scalar.  The number of simulated datasets.
#' @param rbc_args A named list of arguments to be passed to \code{rbc},
#'   including any arguments to be passed to \code{sim_fn} via \code{...}
#'   in \code{rbc}.
#' @param bcthresh_args A named list of arguments to be passed to
#'   \code{bcthresh}.  In particular: \code{probs, lambda, n, prior}.
#'   The defaults for \code{probs, lambda} and \code{n} are
#'   \code{probs = seq(0, 0.9, 0.1)}, \code{lambda = 1} and \code{n = 1000}.
#'   By default \code{n_v = 1}.  If \code{n_v} is supplied then this value is
#'   ignored, without warning.  Similarly, if \code{data} is supplied then
#'   this is also ignored, without warning.
#' @details Add details.
#' @return An object (a list) of class \code{"bc_sim_study"} containing the
#'   following components:
#'     \item{pred_perf }{An array with dimension
#'       (\code{length(probs)}, \code{sims}, \code{length(lambda)).}}
#'     \item{bcthresh_args }{The arguments passed to \code{\link{bcthresh}}.}
#'     \item{rbc_args }{The arguments passed to \code{\link{rbc}}.}
#'     \item{call }{The call to \code{bc_sim_study}.}
#' @examples
#' # Set a prior: flat for GP parameters, Haldane for P(exceedance)
#' prior_args <- list(prior = "flatflat", bin_prior = "haldane",
#'                    h_prior = list(min_xi = -Inf))
#' bcthresh_args <- c(list(lambda = c(0.5, 1, 1.5)), prior_args)
#' rbc_args <- list(n = 1000, lambda = 1)
#' res <- bc_sim_study(2, rbc_args = rbc_args, bcthresh_args = bcthresh_args)
#' plot(res)
#' @export
bc_sim_study <- function(sims, rbc_args, bcthresh_args) {
  # Record the call for later use
  Call <- match.call()
#  rbc_args <- c(list(n = m, lambda = true_lambda, sim_fn = sim_fn),
#                sim_fn_args)
#  # Set arguments to and bcthresh()
#  bcthresh_args <- list(...)
  if (!is.null(bcthresh_args$n_v)) {
    bcthresh_args$n_v <- 1
  }
  if (is.null(bcthresh_args$probs)) {
    bcthresh_args$probs <- seq(0, 0.9, 0.1)
  }
  if (is.null(bcthresh_args$lambda)) {
    bcthresh_args$lambda <- 1
  }
  if (is.null(bcthresh_args$n)) {
    bcthresh_args$n <- 1000
  }
  n_lambda <- length(bcthresh_args$lambda)
  # Set up an array in which to store the results
  res_array <- array(dim = c(length(bcthresh_args$prob), sims, n_lambda))
  for (i in 1:sims) {
    # Simulate a new dataset
    x <- do.call(rbc, rbc_args)
    # Overwrite data with the newly-simulated data
    bcthresh_args$data <- x
    # Call bcthresh() using the new dataset
    res <- do.call(bcthresh, bcthresh_args)
    # Store the lambda-specific results
    for (j in 1:n_lambda) {
      res_array[, i, j] <- res$pred_perf[, , j]
    }
  }
  bcthresh_args$data <- NULL
  temp <- list(pred_perf = res_array, rbc_args = rbc_args,
               bcthresh_args = bcthresh_args, call = Call)
  class(temp) <- "bc_sim_study"
  return(temp)
}

# ============================== plot.bc_sim_study ============================

#' Plot method for objects of class "bc_sim_study"
#'
#' \code{plot} method for class "bc_sim_study". Plots the summaries of
#' predictive performance against quantile of training threshold for different
#' values of \eqn{\lambda}: those stored in \code{x$lambda}.
#'
#' @param x an object inheriting from class "bc_sim_study", a result of a call to
#'   \code{\link{bc_sim_study}}.
#' @param type A numeric vector.  Determines the statistic(s) used to
#'   summarise measures of performance across different simulated datasets,
#'   for each each combination of \eqn{\lambda} and threshold level.
#'   If \code{type = 0} then the mean is used.  Otherwise, sample quantiles
#'   are used and \code{type} is a numeric vector of values in (0, 1) that
#'   is passed as the argument \code{probs} to \code{\link[stats]{quantile}}.
#'   For example, if \code{type = 0.5} then the sample median is used.
#' @param which_lambdas A numeric vector.  Specifies which values of
#'   \eqn{\lambda}, that is, the components of \code{x$bcthresh_args$lambda},
#'   to include in the plot.  The default is to use all these values.
#' @param legend_pos The position of the legend specified using the argument
#'   \code{x} in \code{\link[graphics]{legend}}.
#' @param ... Additional arguments to be passed to
#'   \code{\link[graphics]{matplot}} and/or \code{\link[graphics]{legend}}.
#'   Lines that relate to a common value of \eqn{\lambda} will share
#'   the same value of \code{lty}, \code{lwd} and \code{col}.
#'   The default setting plots solid lines (\code{lty = 1}) of width 2
#'   (\code{lwd = 2}) using \code{col = 1:length(x$bcthresh_args$lambda)}.
#' @details If \code{type = 0} then it is possible that a curve on the plot
#'   may be incomplete.  This indicates that, for a particular \eqn{\lambda} and
#'   threshold level, a measure of predictive performance is \code{-Inf} for
#'   at least one simulated dataset.  It occurs when an observation in the
#'   data lies above the estimated upper end point of the predictive
#'   distribution produced when this observation is removed.
#' @return Nothing.
#' @section Examples:
#' See the examples in \code{\link{bc_sim_study}}.
#' @export
plot.bc_sim_study <- function(x, type = 0,
                              which_lambdas = 1:length(x$bcthresh_args$lambda),
                              legend_pos = "bottom", ...) {
  # Choose the values of lambda
  x$pred_perf <- x$pred_perf[, , which_lambdas, drop = FALSE]
  lambda <- x$bcthresh_args$lambda[which_lambdas]
  n_lambda <- length(lambda)
  n_probs <- length(x$bcthresh_args$probs)
  if (length(type) == 1 && type == 0) {
    my_ylab <- "mean CV performance"
    ymat <- apply(x$pred_perf, MARGIN = c(1, 3), mean)
  } else if (length(type) == 1 && type == 0.5) {
    my_ylab <- "median CV performance"
    ymat <- apply(x$pred_perf, MARGIN = c(1, 3), stats::median)
  } else {
    my_ylab <- "quantile(s) of CV performance"
    my_fn <- function(type) {
      return(apply(x$pred_perf, MARGIN = c(1, 3), stats::quantile,
                   probs = type))
    }
    # This returns an array in which the 3rd dimension indexes probs
    ymat <- vapply(type, my_fn, matrix(0, nrow = n_probs, ncol = n_lambda))
    # Collapse to a matrix for plotting
    dim(ymat) <- c(n_probs, n_lambda * length(type))
  }
  my_lty <- 1
  my_col <- 1:n_lambda
  my_matplot <- function(x, y, ..., lty = my_lty, col = my_col, lwd = 2,
                         xlab = "quantile of training threshold / %",
                         ylab = my_ylab, type = "l") {
    lty <- rep(rep_len(lty, n_lambda), times = length(type))
    col <- rep(rep_len(col, n_lambda), times = length(type))
    lwd <- rep(rep_len(lwd, n_lambda), times = length(type))
    graphics::matplot(x, y, ..., lty = lty, col = col, lwd = lwd, xlab = xlab,
                      ylab = ylab, type = type)
  }
  my_legend <- function(..., x = legend_pos,
                        legend = paste0("lambda = ", lambda), lty = my_lty,
                        col = my_col, lwd = 2) {
    lty <- rep_len(lty, n_lambda)
    col <- rep_len(col, n_lambda)
    lwd <- rep_len(lwd, n_lambda)
    graphics::legend(x = x, legend = legend, lty = lty, col = col, lwd = lwd)
  }
  # Apply the averaging function over dimensions 1 and 3 of x$pred_perf
  # i.e. average over simulations (dimension 2) for for different threshold
  # levels (dimension 1) and values of lambda (dimension 3)
  my_matplot(x$bcthresh_args$probs, ymat, ...)
  my_legend(...)
  return(invisible())
}
