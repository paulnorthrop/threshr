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
  # Numbers of training and validation thresholds
  n_u <- length(x$u_vec)
  n_v <- length(x$v_vec)
  # Calculate threshold weights
  # 1. Find the largest value of predictive performance, and the lambda
  which_max <- arrayInd(which.max(x$pred_perf), .dim = dim(x$pred_perf))
  max_val <- x$pred_perf[which_max]
  which_lambda <- which_max[3]
  # 2. Shoof using the means of the values for this value of lambda
  temp <- x$pred_perf[, , which_lambda, drop = FALSE]
  shoof <- matrix(colMeans(temp * !is.infinite(temp), na.rm = TRUE),
                  ncol = n_v, nrow = n_u, byrow = TRUE)
  cat("Threshold weights:\n")
  for (i in 1:length(x$lambda)) {
    cat("lambda =", x$lambda[i], "\n")
    tw <- apply(exp(x$pred_perf[, , i] - shoof), 2,
                function(x) x / sum(x, na.rm = TRUE))
    rownames(tw) <- x$u_ps
    colnames(tw) <- x$v_ps
    print.default(format(t(tw), digits = digits), print.gap = 1L,
                  quote = FALSE)
  }
  return(invisible(x))
}

# ================================ plot.bcthresh =============================

#' Plot method for objects of class "bcthresh"
#'
#' \code{plot} method for class "bcthresh".
#'
#' @param x an object inheriting from class "bcthresh", a result of a call to
#'   \code{\link{bcthresh}}.
#' @param y Not used.
#' @param which_v A numeric scalar.  Specifies the validation threshold, that
#'   is the components of \code{x$v_vec}, to use in the plot.
#' @param legend_pos The position of the legend specified using the argument
#'   \code{x} in \code{\link[graphics]{legend}}.
#' @param ... Additional arguments to be passed to
#'   \code{\link[graphics]{matplot}}.
#' @details Plots the measure of predictive performance against threshold
#' @return Nothing.
#' @section Examples:
#' See the examples in \code{\link{bcthresh}}.
#' @export
plot.bcthresh <- function(x, y, which_v = 1, legend_pos = "bottom", ...) {
  n_v <- length(x$v_vec)
  # Check that which_v is admissible
  if (!is.numeric(which_v) || !(which_v %in% 1:n_v)) {
    stop("'which_v' must be in 1:length(x$v_vec)")
  }
  my_matplot <- function(x, y, ...,
                         xlab = "quantile of training threshold / %",
                         ylab = "CV performance", type = "l") {
    graphics::matplot(x, y, ..., xlab = xlab, ylab = ylab, type = type)
  }
  # x[, i, , drop = FALSE] gives a matrix of predictive performances
  # Each column gives the values for a different value of lambda
  ymat <- x$pred_perf[, which_v, ]
  my_matplot(x$u_ps, ymat, ...)
  n_lambda <- length(x$lambda)
  graphics::legend(legend_pos, legend = paste0("lambda = ", x$lambda),
                   lty = 1:n_lambda, col = 1:n_lambda, lwd = 2)
  return(invisible())
}
