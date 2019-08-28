# ================================= rboxcox ===================================

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
#' x <- rboxcox(n = 100, lambda = 2)
#' x <- rboxcox(n = 100, lambda = 2, sim_fn = revdbayes::rgp)
#' x <- rboxcox(n = 100, lambda = 2, sim_fn = revdbayes::rgp, shape = 0.1)
#' @export
rboxcox <- function(n = 1, lambda = 1, sim_fn = stats::rexp, ...) {
  return(inv_bc(sim_fn(n, ...), lambda = lambda))
}

# ============================== bc_sim_study =================================

# Change arguments to lists: rboxcox_args and bcthresh_args?
# Less ambiguity but clunkier.

#' Simulation study
#'
#' Add description
#' @param sims An integer scalar.  The number of simulated datasets.
#' @param m,true_lambda,sim_fn Arguments to be passed to \code{\link{rboxcox}}:
#'   \code{m} is passed to \code{\link{rboxcox}} as \code{n} (that is, the
#'   size of the simulated sample),
#'   \code{true_lambda} is passed to \code{\link{rboxcox}} as \code{lambda}.
#' @param sim_fn_args A list of arguments to \code{sim_fn} to be passed
#'   to \code{\link{rboxcox}}.
#' @param ... Further arguments to be passed to \code{bcthresh}.
#'   In particular: \code{probs, lambda, n, prior}.
#'   The defaults for \code{probs, lambda} and \code{n} are
#'   \code{probs = seq(0, 0.9, 0.1)}, \code{lambda = 1} and \code{n = 1000}.
#' @details Add details.
#'   By default \code{n_v = 1}.  If \code{n_v} is supplied via \code{...} then
#'   this value is ignored, without warning.
#' @return An object (a list) of class \code{"bc_sim_study"} containing the
#'   following components:
#'     \item{pred_perf }{An array with dimension
#'       (\code{length(probs)}, \code{sims}, \code{length(lambda)).}}
#' @examples
#' # Set Haldane and flatflat priors
#' # Set a prior: flat for GP parameters, Haldane for P(exceedance)
#' prior_args <- list(prior = "flatflat", bin_prior = "haldane",
#'                    h_prior = list(min_xi = -Inf))
#' bcthresh_args <- c(list(lambda = c(0.5, 1, 1.5)), prior_args)
#' rboxcox_args <- list(m = 1000, true_lambda = 1)
#' res <- bc_sim_study(2, lambda = c(0.5, 1, 1.5))
#' plot(res)
#' @export
bc_sim_study <- function(sims, m = 1000, true_lambda = 1, sim_fn = stats::rexp,
                         sim_fn_args = list(), ...) {
  # Record the call for later use
  Call <- match.call()
  # Set arguments to rboxcox()
  rboxcox_args <- c(list(n = m, lambda = true_lambda, sim_fn = sim_fn),
                    sim_fn_args)
  # Set arguments to and bcthresh()
  bcthresh_args <- list(...)
  if (!is.null(bcthresh_args$data)) {
    stop("''data'' cannot be supplied in ...: data are simulated")
  }
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
    x <- do.call(rboxcox, rboxcox_args)
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
  temp <- list(pred_perf = res_array, rboxcox_args = rboxcox_args,
               bcthresh_args = bcthresh_args, call = Call)
  class(temp) <- "bc_sim_study"
  return(temp)
}

# ============================== plot.bc_sim_study ============================

#' Plot method for objects of class "bc_sim_study"
#'
#' \code{plot} method for class "bc_sim_study".
#'
#' @param x an object inheriting from class "bc_sim_study", a result of a call to
#'   \code{\link{bc_sim_study}}.
#' @param sum Type of summary of the CV measures.
#' @param probs If \code{sum = "quantiles"} then \code{probs} provides
#'   the argument to \code{\link[stats]{quantile}}.
#' @param which_lambdas A numeric vector.  Specifies which values of
#'   \eqn{\lambda}, that is, the components of \code{x$lambda}, to use in the
#'   plot.  The default is to use all these values.
#' @param legend_pos The position of the legend specified using the argument
#'   \code{x} in \code{\link[graphics]{legend}}.
#' @param ... Additional arguments to be passed to
#'   \code{\link[graphics]{matplot}}.
#' @details Plots the measure of predictive performance against quantile
#'   of training threshold for different values of \eqn{\lambda}: those
#'   stored in \code{x$lambda}.
#' @return Nothing.
#' @section Examples:
#' See the examples in \code{\link{bc_sim_study}}.
#' @export
plot.bc_sim_study <- function(x, sum = c("mean", "median", "quantiles"),
                              probs = 0.5,
                              which_lambdas = 1:length(x$bcthresh_args$lambda),
                              legend_pos = "bottom", ...) {
  sum <- match.arg(sum)
  # Choose the values of lambda
  x$pred_perf <- x$pred_perf[, , which_lambdas, drop = FALSE]
  lambda <- x$bcthresh_args$lambda[which_lambdas]
  n_lambda <- length(lambda)
  n_probs <- length(x$bcthresh_args$probs)
  if (sum == "mean") {
    my_ylab <- "mean CV performance"
    ymat <- apply(x$pred_perf, MARGIN = c(1, 3), mean)
    print(dim(ymat))
  } else if (sum == "median") {
    my_ylab <- "median CV performance"
    ymat <- apply(x$pred_perf, MARGIN = c(1, 3), stats::median)
  } else {
    my_ylab <- "quantile(s) of CV performance"
    my_fn <- function(probs) {
      return(apply(x$pred_perf, MARGIN = c(1, 3), stats::quantile,
                   probs = probs))
    }
    # This returns an array in which the 3rd dimension indexes probs
    ymat <- vapply(probs, my_fn, matrix(0, nrow = n_probs, ncol = n_lambda))
    print(ymat)
    dim(ymat) <- c(n_probs, n_lambda * length(probs))
    print(ymat)
  }
  #
  my_lty <- rep(1:n_lambda, times = length(probs))
  my_col <- rep(1:n_lambda, times = length(probs))
  my_matplot <- function(x, y, ..., lty = my_lty, col = my_col,
                         xlab = "quantile of training threshold / %",
                         ylab = my_ylab, type = "l") {
    graphics::matplot(x, y, ..., lty = lty, col = col, xlab = xlab,
                      ylab = ylab, type = type)
  }
  # Apply the averaging function over dimensions 1 and 3 of x$pred_perf
  # i.e. average over simulations (dimension 2) for for different threshold
  # levels (dimension 1) and values of lambda (dimension 3)
  my_matplot(x$bcthresh_args$probs, ymat, ...)
  graphics::legend(legend_pos, legend = paste0("lambda = ", lambda),
                   lty = my_lty, col = my_col, lwd = 2)
  return(invisible())
}
