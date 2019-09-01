# =========================== plot.bcthreshpred ===========================

#' Plot diagnostics an bcthreshpred object
#'
#' \code{plot} method for class "bcthreshpred".  Produces plots to summarise
#' the predictive inferences made by \code{\link{predict.bcthresh}}.
#'
#' @param x an object of class "bthreshpred", a result of a call to
#'   \code{\link{ithresh}}.
#' @param which_lambdas A numeric vector.  Specifies which values of
#'   \eqn{\lambda}, that is, the components of \code{x$lambda},
#'   to include in the plot.  The default is to use all these values.
#' @param legend_pos The position of the legend specified using the argument
#'   \code{x} in \code{\link[graphics]{legend}}.  If \code{legend_pos} is
#'   missing then it is set internally, based on \code{x$type}.
#' @param ... Additional arguments passed to \code{\link[graphics]{matplot}}.
#' @export
plot.bcthreshpred <- function(x,
                              which_lambdas = 1:length(x$lambda),
                              legend_pos, ...) {
  if (!inherits(x, "bcthreshpred")) {
    stop("use only with \"bcthreshpred\" objects")
  }
  if (missing(legend_pos)) {
    legend_pos <- switch(x$type,
                         d = "topright",
                         p = "topleft",
                         q = "topleft")
  }
  lambda <- x$lambda[which_lambdas]
  my_col <- 1:length(lambda)
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
  if (x$type %in% c("d", "p", "q")) {
    x$x <- x$x[, which_lambdas, drop = FALSE]
    x$y <- x$y[, which_lambdas, drop = FALSE]
    my_matplot(x$x, x$y, ...)
    my_legend()
  }
  return(invisible())
}
