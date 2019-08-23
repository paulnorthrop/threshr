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
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
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
  twmax <- apply(exp(x$pred_perf[, , which_lambda] - shoof), 2, sum)
  print(twmax)
  cat("Threshold weights:\n")
  for (i in 1:length(x$lambda)) {
    cat("lambda =", x$lambda[i], "\n")
    tw <- apply(exp(x$pred_perf[, , i] - shoof), 2,
                function(x) x / sum(x, na.rm = TRUE))
#    tw <- apply(exp(x$pred_perf[, , i] - shoof), 2,
#                function(x) x / twmax)
    print(sum(tw))
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
#'   is the indices of the argument \code{v_vec} to \code{\link{bcthresh}},
#'   to use in the plot.
#' @param ... Additional arguments to be passed to
#'   \code{\link[graphics]{matplot}}.
#' @details Add details.
#' @return Nothing.
#' @export
plot.bcthresh <- function(x, y, which_v = 1, ...) {
  n_v <- length(x$v_vec)
  # Check that which_v is admissible
  if (!is.numeric(which_v) || !(which_v %in% 1:n_v)) {
    stop("'which_v' must be in 1:length(object$v_vec)")
  }

  my_matplot <- function(x, y, ..., xlab = "threshold", ylab = "performance",
                         type = "l") {
    graphics::matplot(x, y, ..., xlab = xlab, ylab = ylab, type = type)
  }
  # x[, i, , drop = FALSE] gives a matrix of predictive performances
  # Each column gives the values for a different value of lambda
  ymat <- x$pred_perf[, which_v, ]
  print(ymat)
  my_matplot(x$u_vec, ymat, ...)
  n_lambda <- length(x$lambda)
  graphics::legend("bottom", legend = paste0("lambda = ", x$lambda),
                   lty = 1:n_lambda, col = 1:n_lambda, lwd = 2)
  return(invisible())
}
