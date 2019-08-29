# ================================ print.bcthresh =============================

#' Print method for objects of class "bcthresh"
#'
#' \code{print} method for class "bcthresh".
#'
#' @param x an object inheriting from class "bcthresh", a result of a call to
#'   \code{\link{bcthresh}}.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base]{format}} and \code{\link[base:Round]{signif}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @details For each value in \code{x$lambda}, a matrix of the estimated
#'   threshold weights is printed.
#'   Each row gives the weights for each training threshold for a given
#'   validation threshold.  The row and column names are the approximate
#'   quantile levels of the thresholds.
#'
#'   \strong{Each row is normalised to sum to one.  Therefore, this output
#'   cannot be used to compare the peformances of different values of
#'   \eqn{\lambda}, only different training thresholds for a given
#'   \eqn{\lambda} (and validation threshold).}
#'
#'   Use \code{\link{plot.bcthresh}} to compare performance measures across
#'   different values of \eqn{\lambda}.
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#' @section Examples:
#' See the examples in \code{\link{bcthresh}}.
#' @export
print.bcthresh <- function(x, digits = 2, ...) {
  if (!inherits(x, "bcthresh")) {
    stop("use only with \"bcthresh\" objects")
  }
  cat("Threshold weights:\n")
  n_lambda <- length(x$lambda)
  for (lambda in x$lambda) {
    temp <- choose_lambda(x, lambda = lambda)
    class(temp) <- c("ithresh", "bcthresh")
    cat("lambda =", lambda, "\n")
    print(temp)
  }
  return(invisible(x))
}

# =============================== summary.bcthresh ============================

#' Summarizing measures of threshold predictive performance
#'
#' \code{summary} method for class "bcthresh"
#'
#' @param object an object inheriting from class "bcthresh", a result of a call
#'   to \code{\link{bcthresh}}.
#' @param digits An integer. Used for number formatting with
#'   \code{\link[base]{format}} and \code{\link[base:Round]{signif}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @return Returns a numeric matrix with 6 columns and \code{n_v} rows,
#'   where \code{n_v} is an argument to \code{\link{ithresh}} that
#'   determines how many of the largest training thresholds are used
#'   a validation thresholds.  The first column is the value of \eqn{\lambda}.
#'   The remaining 5 columns are those returned from
#'   \code{\link{summary.ithresh}}.
#' @section Examples:
#' See the examples in \code{\link{bcthresh}}.
#' @export
summary.bcthresh <- function(object, digits = 2, ...) {
  if (!inherits(object, "bcthresh")) {
    stop("use only with \"bcthresh\" objects")
  }
  res <- NULL
  n_lambda <- length(object$lambda)
  for (lambda in object$lambda) {
    temp <- choose_lambda(object, lambda = lambda)
    res <- rbind(res, summary(temp))
  }
  lambda_vec <- rep(object$lambda, each = length(object$v_vec))
  res <- cbind(lambda = lambda_vec, res)
  return(res)
}

# ================================ plot.bcthresh =============================

#' Plot method for objects of class "bcthresh"
#'
#' \code{plot} method for class "bcthresh".
#'
#' @param x an object inheriting from class "bcthresh", a result of a call to
#'   \code{\link{bcthresh}}.
#' @param which_v A numeric scalar.  Specifies the validation threshold, that
#'   is the components of \code{x$v_vec}, to use in the plot.
#' @param which_lambdas A numeric vector.  Specifies which values of
#'   \eqn{\lambda}, that is, the components of \code{x$lambda}, to include in
#'   the plot.  The default is to use all these values.
#' @param legend_pos The position of the legend specified using the argument
#'   \code{x} in \code{\link[graphics]{legend}}.
#' @param ... Additional graphical parameters to be passed to
#'   \code{\link[graphics]{matplot}} and \code{\link[graphics]{legend}}.
#'   The default setting plots solid lines (\code{lty = 1}) of width 2
#'   (\code{lwd = 2}) using \code{col = 1:length(x$lambda)}.
#' @details Plots the measure of predictive performance against quantile
#'   of training threshold for different values of \eqn{\lambda}: those
#'   stored in \code{x$lambda}.  It is possible that a curve on the plot
#'   may be incomplete.  This indicates that, for a particular \eqn{\lambda}
#'   and threshold level, a measure of predictive performance is \code{-Inf}.
#'   This occurs when an observation in the data lies above the estimated upper
#'   end point of the predictive distribution produced when this observation is
#'   removed.
#' @return Nothing.
#' @section Examples:
#' See the examples in \code{\link{bcthresh}}.
#' @export
plot.bcthresh <- function(x, which_v = 1,
                          which_lambdas = 1:length(x$lambda),
                          legend_pos = "bottom", ...) {
  # Choose the values of lambda
  x$pred_perf <- x$pred_perf[, , which_lambdas, drop = FALSE]
  the_lambdas <- x$lambda[which_lambdas]
  n_v <- length(x$v_vec)
  # Check that which_v is admissible
  if (!is.numeric(which_v) || !(which_v %in% 1:n_v)) {
    stop("'which_v' must be in 1:length(x$v_vec)")
  }
  my_lty <- 1
  my_col <- 1:length(x$lambda)
  my_lwd <- 2
  my_matplot <- function(x, y, ...,
                         xlab = "quantile of training threshold / %",
                         ylab = "CV performance", type = "l",
                         lty = my_lty, col = my_col, lwd = my_lwd) {
    graphics::matplot(x, y, ..., xlab = xlab, ylab = ylab, type = type,
                      lty = lty, col = col, lwd = lwd)
  }
  my_legend <- function(..., x = legend_pos,
                        legend = paste0("lambda = ", the_lambdas), lty = my_lty,
                        col = my_col, lwd = my_lwd) {
    graphics::legend(x = x, legend = legend, lty = lty, col = col, lwd = lwd)
  }
  # x[, i, , drop = FALSE] gives a matrix of predictive performances
  # Each column gives the values for a different value of lambda
  ymat <- x$pred_perf[, which_v, ]
  my_matplot(x$u_ps, ymat, ...)
  my_legend(...)
  return(invisible())
}
