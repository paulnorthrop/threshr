# =========================== plot.bcthreshpred ===========================

#' Plot diagnostics an bcthreshpred object
#'
#' \code{plot} method for class "bcthreshpred".  Produces plots to summarise
#' the predictive inferences made by \code{\link{predict.bcthresh}}.
#'
#' @param x an object of class "bthreshpred", a result of a call to
#'   \code{\link{ithresh}}.
#' @param lambda A numeric vector.  Specifies which values of
#'   \eqn{\lambda} to include in the plot.  These values must be present
#'   in \code{x$lambda}.  The default is to use all these values.
#' @param n_years A numeric vector.  Specifies which values in
#'   \code{x$n_years} to use for the time horizon.  If
#'   \code{length(lambda) > 1} then \code{n_years} must have length 1.
#'   The defaults are that \code{x$n_years[1]} is used if
#'   \code{length(lambda) > 1} and \code{x$years} otherwise.
#' @param legend_pos The position of the legend specified using the argument
#'   \code{x} in \code{\link[graphics]{legend}}.  If \code{legend_pos} is
#'   missing then it is set internally, based on \code{x$type}.
#' @param ... Additional arguments passed to \code{\link[graphics]{matplot}}.
#' @export
plot.bcthreshpred <- function(x, lambda = x$lambda,
                              n_years = if (length(lambda) == 1) x$n_years else
                                x$n_years[1],
                              legend_pos, ...) {
  if (!inherits(x, "bcthreshpred")) {
    stop("use only with \"bcthreshpred\" objects")
  }
  if (missing(legend_pos)) {
    legend_pos <- switch(x$type,
                         d = "topright",
                         p = "bottomright",
                         q = "topleft")
  }
  n_lambda <- length(lambda)
  if (!all(lambda %in% x$lambda)) {
    stop("All lambda must be in x$lambda")
  }
  if (n_lambda > 1 && length(n_years) > 1) {
    stop("If length(lambda) > 1 then length(n_years) must be 1")
  }
  which_lambdas <- which(x$lambda %in% lambda)
  which_n_years <- which(x$n_years %in% n_years)
  my_col <- 1:n_lambda
  my_lty <- 1
  my_lwd <- 2
  my_matplot <- function(x, y, ..., type = "l", col = my_col, lty = my_lty,
                         lwd = my_lwd, xlab = my_xlab, ylab = my_ylab) {
    graphics::matplot(x, y, ..., type = type, col = col, lty = lty,
                      lwd = lwd, xlab = xlab, ylab = ylab)
  }
  my_legend <- function(..., x = legend_pos,
                        legend = paste0("lambda = ", lambda), lty = my_lty,
                        col = my_col, lwd = my_lwd) {
    graphics::legend(x = x, legend = legend, lty = lty, col = col, lwd = lwd)
  }
  if (x$type == "d") {
    my_xlab <- "quantile"
    my_ylab <- "density"
  }
  if (x$type == "p") {
    my_xlab <- "quantile"
    my_ylab <- "probability"
  }
  if (x$type == "q") {
    my_xlab <- "probability"
    my_ylab <- "quantile"
  }
  # If we averaged over thresholds then use only the threshold-averaged values
  if (x$type %in% c("d", "p", "q")) {
    if (x$which_u == "all") {
      x$x <- x$x[, 1, which_lambdas]
      x$y <- x$y[, length(x$u_vec), which_lambdas]
    } else {
      x$x <- x$x[, which_n_years, which_lambdas]
      x$y <- x$y[, which_n_years, which_lambdas]
    }
    my_matplot(x$x, x$y, ...)
    my_legend()
  }
  return(invisible())
}
