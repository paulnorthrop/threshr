# =================================== rbc =====================================

#' Simulation with (Inverse) Box-Cox Transformation
#'
#' Random generation for a random variable that has the following property:
#' a Box-Cox transformation with parameter \code{lambda} of the variable has a
#' specified distribution.
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

MN_median <- function(lambda = 1, npy, N, quantile_fn = stats::qexp, ...) {
  return(inv_bc(quantile_fn(-log(2) / (N * npy), log.p = TRUE, ...),
         lambda = lambda))
}

# ============================== bc_sim_study =================================

#' Box-Cox threshold selection simulation study
#'
#' Performs a simulation study to compare the extreme value predictive
#' performance over a range of training threshold levels after different
#' Box-Cox transformation.
#' @param sims An integer scalar.  The number of simulated datasets.
#' @param rbc_args A named list of arguments to be passed to \code{rbc},
#'   including any arguments to be passed to \code{sim_fn} via \code{...}
#'   in \code{rbc}.
#' @param MN_args A named list of arguments containing the number of
#'   observations per years \code{npy} and the time horizon \code{N}.
#' @param bcthresh_args A named list of arguments to be passed to
#'   \code{bcthresh}.  In particular: \code{probs, lambda, n, prior}.
#'   The defaults for \code{probs, lambda} and \code{n} are
#'   \code{probs = seq(0, 0.9, 0.1)}, \code{lambda = 1} and \code{n = 1000}.
#'   By default \code{n_v = 1}.  If \code{n_v} is supplied then this value is
#'   ignored, without warning.  Similarly, if \code{data} is supplied then
#'   this is also ignored, without warning.
#' @details The general idea is to simulated repeatedly data from a
#'   distribution for which a certain value of the Box-Cox transformation
#'   parameter will produce known behaviour.  This is achieved by simulating
#'   data using the function \code{\link{rbc}}.
#'
#'   In the example below data are simulated for a random variable that has the
#'   following property: a Box-Cox transformation with parameter
#'   \eqn{\lambda = 2} of the variable has a unit exponential distribution,
#'   which is in the Generalised Pareto family.
#'   Therefore, we should find that, when transforming prior to a
#'   threshold-based extreme value anaysis, using \eqn{\lambda = 2} produces
#'   better performance, for lower threshold levels, that other values of
#'   \eqn{\lambda}.
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
#' bcthresh_args <- c(list(lambda = c(1, 1.5, 2, 2.5, 3)), prior_args)
#' rbc_args <- list(n = 1000, lambda = 2)
#' MN_args <- list(npy = 20, N = 10 ^ seq(2, 4, len = 15))
#' res <- bc_sim_study(2, rbc_args, bcthresh_args, MN_args)
#' plot(res)
#' @export
bc_sim_study <- function(sims, rbc_args, bcthresh_args, MN_args) {
  # Record the call for later use
  Call <- match.call()
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
  # Set up an array in which to store the predctive performance results
  res_array <- array(dim = c(length(bcthresh_args$prob), sims, n_lambda))
  # ... and for predictive medians of N-year maxima, best threshold
  best_array <- array(dim = c(length(MN_args$N), sims, n_lambda))
#  # ... and for predictive medians of N-year maxima, averaged over thresholds
#  ave_array <- array(dim = c(length(MN_args$N), sims, n_lambda))
  # Create a sims by n_lambda matrix in which to store the best thresholds
  best_u <- matrix(NA, ncol = n_lambda, nrow = sims)
  for (i in 1:sims) {
    if (i %% 10 == 0) {
      print(i)
    }
    # Simulate a new dataset, overwrite bcthresh_args$data and call bcthresh()
    bcthresh_args$data <- do.call(rbc, rbc_args)
    res <- do.call(bcthresh, bcthresh_args)
    # Store lambda-specific results
    for (j in 1:n_lambda) {
      res_array[, i, j] <- res$pred_perf[, , j]
    }
    # Estimate predictive medians of N-year maxima for best threshold
    qbest <- predict(res, npy = MN_args$npy, n_years = MN_args$N,
                     which_u = "best", type = "q", x = 1 / 2)
    for (j in 1:n_lambda) {
      best_array[, i, j] <- qbest$y[, , j]
    }
    best_u[i, ] <- qbest$best_u
#    for (k in 1:length(MN_args$N)) {
#      qall <- predict(res, npy = MN_args$npy, n_years = MN_args$N[k],
#                      which_u = "all", type = "q", x = 1 / 2)
#      for (j in 1:n_lambda) {
#        all_array[, i, j] <- qall$y[, , j]
#      }
#    }
  }
  # Remove data from the returned list bcthresh_args
  bcthresh_args$data <- NULL
  # Calculate the true medians
  MN_args <- c(MN_args, lambda = rbc_args$lambda)
  true_medians <- do.call(MN_median, MN_args)
  temp <- list(pred_perf = res_array, rbc_args = rbc_args,
               bcthresh_args = bcthresh_args, MN_args = MN_args,
               best_array = best_array, best_u = best_u,
               true_medians = true_medians, call = Call)
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
#' @param type An integer scalar.  Determine the type of plot.
#' \itemize{
#'   \item{type = 1: }{plots averages measures of predictive performance
#'     against threshold level by the value of \eqn{\lambda}.}
#'   \item{type = 2: }{plots kernel density estimates of the predictive
#'     density of the median of \eqn{N}-year maxima by the value of
#'     \eqn{\lambda}, using the best threshold for each \eqn{\lambda}, with the
#'     true value of the median indicated using a vertical line.}
#' }
#' @param stat A numeric scalar.  Determines the statistic(s) used to
#'   summarise measures of performance across different simulated datasets,
#'   for each each combination of \eqn{\lambda} and threshold level.
#'   If \code{stat = -1} then the mean is used.  Otherwise, the
#'   100\code{stat}\% sample quantile is used, with \code{stat} being passed
#'   as the argument \code{probs} to \code{\link[stats]{quantile}}.
#'   For example, if \code{stat = 0.5} then the sample median is used.
#' @param which_N A numeric scalar  Specifies which value of \eqn{N}, that is,
#'   which component of \code{x$MN_args$N} is used as the time horizon when
#'   \code{type = 2}.
#' @param which_lambdas A numeric vector.  Specifies which values of
#'   \eqn{\lambda}, that is, the components of \code{x$bcthresh_args$lambda},
#'   to include in the plot.  The default is to use all these values.
#' @param normalise A logical scalar.  Controls what is plotted on the
#'   vertical axis of the plot.  If \code{normalise = FALSE} then the
#'   measure of performance defined in equation (7) of Northrop et al. (2017)
#'   is plotted. Otherwise, these measures are normalised in the manner of
#'   equation (14) Northrop et al. (2017).  The same normalisation constant
#'   is used across of values of \eqn{\lambda} that appear in the plot,
#'   with the property that the values sum to 1 for one value of \eqn{\lambda}
#'   and less than 1 for all other values of \eqn{\lambda}.
#' @param legend_pos The position of the legend that indicates which curves
#'   relates to which value of \eqn{\lambda}, to be passed as the argument
#'   \code{x} in \code{\link[graphics]{legend}}.  The default is
#'   \code{"bottom"} for \code{type = 1} and \code{"topleft"} for
#'   \code{type = 2}.
#' @param summary A logical scalar.  If \code{type = 2} should we add a legend,
#'   in the top right of the plot, giving the estimated bias, standard
#'   deviation and root mean squared error for each value of \eqn{lambda}?
#' @param digits An integer scalar.  The argument \code{digits} to be passed
#'   to \code{\link[base:Round]{signif}} when rounding the
#'   statistics produced when \code{summary = TRUE}.
#' @param ... Additional graphical parameters to be passed to
#'   \code{\link[graphics]{matplot}} and \code{\link[graphics]{legend}}.
#'   The default setting plots solid lines (\code{lty = 1}) of width 2
#'   (\code{lwd = 2}) using \code{col = 1:length(x$bcthresh_args$lambda)}.
#' @details If \code{stat = -1} (the plot contains sample means) or if
#'   \code{stat} is close to zero (plot contains estimates of low quantiles)
#'   then it is possible that a curve on the plot may be incomplete.  This
#'   indicates that, for a particular \eqn{\lambda} and threshold level, a
#'   measure of predictive performance is \code{-Inf} for at least one
#'   simulated dataset. This occurs when an observation in the data lies above
#'   the estimated upper end point of the predictive distribution produced when
#'   this observation is removed.
#' @return Nothing.
#' @section Examples:
#' See the examples in \code{\link{bc_sim_study}}.
#' @export
plot.bc_sim_study <- function(x, type = 1, stat = -1, which_N = 1,
                              which_lambdas = 1:length(x$bcthresh_args$lambda),
                              normalise = FALSE, legend_pos, summary = FALSE,
                              digits = 2, ...) {
  if (type != 1 & type != 2) {
    stop("type must be equal to 1 or 2")
  }
  if (missing(legend_pos)) {
    legend_pos <- ifelse(type == 1, "bottom", "topleft")
  }
  if (type == 1) {
    bcsim_type1_plot(x, stat, which_lambdas, normalise, legend_pos, ...)
  } else {
    bcsim_type2_plot(x, which_N, which_lambdas, legend_pos, summary, digits,
                     ...)
  }
  return(invisible())
}

# Plots of type 1

bcsim_type1_plot <- function(x, stat, which_lambdas, normalise, legend_pos,
                             ...) {
  # Choose the values of lambda
  x$pred_perf <- x$pred_perf[, , which_lambdas, drop = FALSE]
  lambda <- x$bcthresh_args$lambda[which_lambdas]
  n_lambda <- length(lambda)
  n_probs <- length(x$bcthresh_args$probs)
  if (length(stat) > 1 || !is.numeric(stat)) {
    stop("''stat'' must be a numeric scalar")
  }
  my_ylim <- NULL
  if (stat == -1) {
    my_ylab <- "mean CV performance"
    ymat <- apply(x$pred_perf, MARGIN = c(1, 3), mean)
    if (normalise){
      ymat_shoof <- ymat - mean(ymat * !is.infinite(ymat), na.rm = TRUE)
      ymat <- exp(ymat_shoof) / sum(exp(ymat_shoof), na.rm = TRUE)
      ymat <- ymat / max(colSums(ymat, na.rm = TRUE), na.rm = TRUE)
      my_ylab <- paste0("normalised mean CV performance")
      my_ylim <- c(0, max(ymat, na.rm = TRUE))
    }
  } else {
    my_ylab <- paste0(100 * stat, "% quantile of CV performance")
    ymat <- apply(x$pred_perf, MARGIN = c(1, 3), stats::quantile, probs = stat)
    if (normalise){
      ymat_shoof <- ymat - mean(ymat * !is.infinite(ymat), na.rm = TRUE)
      ymat <- exp(ymat_shoof) / sum(exp(ymat_shoof), na.rm = TRUE)
      ymat <- ymat / max(colSums(ymat, na.rm = TRUE), na.rm = TRUE)
      my_ylab <- paste0("normalised ",100*stat,"% quantile of CV performance")
      my_ylim <- c(0, max(ymat, na.rm = TRUE))
    }
  }
  my_lty <- 1
  my_col <- 1:n_lambda
  my_matplot <- function(x, y, ..., lty = my_lty, col = my_col, lwd = 2,
                         xlab = "quantile of training threshold / %",
                         ylab = my_ylab, type = "l", ylim = my_ylim) {
    graphics::matplot(x, y, ..., lty = lty, col = col, lwd = lwd, xlab = xlab,
                      ylab = ylab, type = type, ylim = ylim)
  }
  my_title <- expression(lambda)
  my_legend <- function(..., x = legend_pos, legend = lambda, lty = my_lty,
                        col = my_col, lwd = 2, title = my_title) {
    graphics::legend(x = x, legend = legend, lty = lty, col = col, lwd = lwd,
                     title = title, ...)
  }
  # Apply the averaging function over dimensions 1 and 3 of x$pred_perf
  # i.e. average over simulations (dimension 2) for for different threshold
  # levels (dimension 1) and values of lambda (dimension 3)
  my_matplot(x$bcthresh_args$probs, ymat, ...)
  my_legend(...)
  return(invisible())
}

# Plots of type 2

bcsim_type2_plot <- function(x, which_N, which_lambdas, legend_pos, summary,
                             digits, ...) {
  # Choose the values of lambda
  x$best_array <- as.matrix(x$best_array[which_N, , which_lambdas])
  lambda <- x$bcthresh_args$lambda[which_lambdas]
  n_lambda <- length(lambda)
  N <- x$MN_args$N[which_N]
  true_median <- x$true_medians[which_N]
  dens <- apply(x$best_array, 2, stats::density)
  my_xlim <- range(sapply(dens, "[", "x"))
  my_ylim <- range(sapply(dens, "[", "y"))
  my_lty <- 1
  my_col <- 1:n_lambda
  if (summary) {
    bias <- apply(x$best_array - true_median, 2, mean)
    stdev <- apply(x$best_array - true_median, 2, stats::sd)
    rmse_fun <- function(x) sqrt(mean(x ^ 2))
    RMSE <- apply(x$best_array - true_median, 2, rmse_fun)
    col1 <- c(expression(lambda), lambda)
    col2 <- c("bias", signif(bias, digits))
    col3 <- c("SD", signif(stdev, digits))
    col4 <- c("RMSE", signif(RMSE, digits))
    perf_stats <- c(col1, col2, col3, col4)
  } else {
  }
  my_title <- expression(lambda)
  my_leg <- lambda
  my_plot <- function(x, ..., lty = my_lty, col = my_col, lwd = 2,
                      xlab = paste0("median of ", N, "-year maximum"),
                      ylab = "density", type = "l", xlim = my_xlim,
                      ylim = my_ylim, title = NULL) {
    graphics::plot(NA, ..., xlab = xlab, ylab = ylab, type = type, xlim = xlim,
                   ylim = ylim)
    mapply(graphics::lines, x, lty = lty, col = col, lwd = lwd, ...)
  }
  my_legend <- function(..., x = legend_pos, legend = my_leg, lty = my_lty,
                        col = my_col, lwd = 2, title = my_title) {
    graphics::legend(x = x, legend = legend, lty = lty, col = col, lwd = lwd,
                     title = title, ...)
  }
  # Calculate RMSEs
  my_plot(dens, ...)
  my_legend(...)
  legend("topright", legend = perf_stats, ncol = 4)
  graphics::abline(v = true_median, lwd = 3, lty = 2)
  return(invisible())
}
