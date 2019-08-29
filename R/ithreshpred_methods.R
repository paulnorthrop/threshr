# =========================== plot.ithreshpred ===========================

#' Plot diagnostics an ithreshpred object
#'
#' \code{plot} method for class "ithreshpred".  Produces plots to summarise
#' the predictive inferences made by \code{\link{predict.ithresh}}.
#'
#' @param x an object of class "ithreshpred", a result of a call to
#'   \code{\link{ithresh}}.
#' @param y Not used.
#' @param ... Additional arguments passed on to
#'   \code{\link[revdbayes]{plot.evpred}}.
#' @param ave_only  A logical scalar.  Only relevant if
#'   \code{\link{predict.ithresh}} was called with \code{which_u = "all"}.
#'   If \code{TRUE} then plot only
#'   a curve for the weighted average over multiple training thresholds.
#'   If \code{FALSE} then also plot a curve for each training threshold.
#' @param add_best A logical scalar.  If \code{TRUE} then the best
#'   threshold, as judged using the validation threshold selected using the
#'   argument \code{which_v} supplied to \code{\link{predict.ithresh}}, is
#'   highlighted by plotting it with a different line style.
#' @details \emph{Single threshold case}, where
#'   \code{\link{predict.ithresh}} was called with numeric scalar
#'   \code{which_u} or \code{which_u = "best"}.
#'   \code{\link[revdbayes]{plot.evpred}} is called to produce the plot.
#'
#'   \emph{Multiple threshold} case, where
#'   \code{\link{predict.ithresh}} was called with \code{which_u = "all"}.
#'   Again, \code{\link[revdbayes]{plot.evpred}} is called but now the
#'   estimated predictive distribution function (\code{type = "p"} used
#'   in the call to \code{\link{predict.ithresh}}) or density function
#'   (\code{type = "d"}) is plotted for each of the training thresholds
#'   (grey lines) as is the result of the weighted average over the
#'   different training thresholds (black line).
#'   If graphical parameters, such as \code{lty}, \code{lwd} or \code{col}
#'   are passed via \code{...} then the first element relates to the
#'   weighted average and the remaining \code{length(x$u_vec)} elements to
#'   the respective training thresholds in \code{u_vec}.
#' @return A list containing the graphical parameters using in producing the
#'   plot including any arguments supplied via ... is returned (invisibly).
#' @examples
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
#' npy_gom <- length(gom)/105
#' best_p <- predict(gom_cv, n_years = c(100, 1000))
#' plot(best_p)
#'
#' # Predictive density function
#' best_d <- predict(gom_cv, type = "d", n_years = c(100, 1000))
#' plot(best_d)
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
#'
#' ### ... and highlight the best threshold
#'
#' plot(all_p, add_best = TRUE)
#' plot(all_d, add_best = TRUE)
#' @seealso \code{\link{ithresh}} for threshold selection in the i.i.d. case
#'   based on leave-one-out cross-validation.
#' @seealso \code{\link{predict.ithresh}} for predictive inference for the
#'   largest value observed in N years.
#' @seealso \code{\link{plot.ithresh}} for the S3 plot method for objects of
#'   class \code{ithresh}.
#' @seealso \code{\link{summary.ithresh}} Summarizing measures of threshold
#'   predictive performance.
#' @export
plot.ithreshpred <- function(x, y, ..., ave_only = FALSE, add_best = FALSE) {
  if (!inherits(x, "ithreshpred")) {
    stop("use only with \"ithreshpred\" objects")
  }
  # Single threshold
  if (x$which_u == "best" || is.numeric(x$which_u)) {
    temp <- x
    class(temp) <- "evpred"
    # plot is revdbayes:::plot.evpred
    revdbayes:::plot.evpred(temp, ...)
    for_plot <- list(...)
    return(invisible(for_plot))
  }
  best <- 1 + x$best_u
  # Multiple thresholds
  if (x$which_u == "all") {
    temp <- x
    class(temp) <- "evpred"
    n_col_y <- ncol(temp$y)
    # Create list of arguments for revdbayes' plot.evpred()
    # Set lty, lwd and col if they have not be supplied by the user
    # If they have been supplie dthen make them the correct length
    for_plot <- list(x = temp, ...)
    if (is.null(for_plot$lty)) {
      for_plot$lty <- rep(1, n_col_y)
      if (add_best) {
        for_plot$lty[best] <- 2
      }
    } else {
      for_plot$lty <- rep_len(for_plot$lty, n_col_y)
    }
    # If we only want one `averaged' line then make all the others disappear
    if (ave_only) {
      for_plot$lty <- c(1, rep(0, n_col_y - 1))
    }
    if (is.null(for_plot$lwd)) {
      for_plot$lwd <- rep(2, n_col_y)
    } else {
      for_plot$lwd <- rep_len(for_plot$lwd, n_col_y)
    }
    if (is.null(for_plot$col)) {
      for_plot$col <- c(1, rep("grey", n_col_y - 1))
      if (add_best) {
        for_plot$col[best] <- 1
      }
    } else {
      for_plot$col <- rep_len(for_plot$col, n_col_y)
    }
    if (is.null(for_plot$leg_text)) {
      if (ave_only) {
        for_plot$leg_text <- c("averaged")
      } else if (add_best) {
        flip <- c(1, best, (1:n_col_y)[-c(1, best)])
        for_plot$lty <- for_plot$lty[flip]
        for_plot$lwd <- for_plot$lwd[flip]
        for_plot$col <- for_plot$col[flip]
        for_plot$leg_text <- c("averaged", "best single", "single")
        for_plot$x$y <- for_plot$x$y[, flip]
      } else {
        for_plot$leg_text <- c("averaged", "single")
      }
    }
    # plot is revdbayes:::plot.evpred
    do.call(revdbayes:::plot.evpred, for_plot)
    # Replot the `average' line, so that it appears on top
    for_lines <- list(x = temp$x, y = temp$y[, 1], lty = for_plot$lty[1],
                      lwd = for_plot$lwd[1], col = for_plot$col[1])
    do.call(graphics::lines, for_lines)
    # Replot the `best' line, so that it appears on top
    if (add_best) {
      for_lines <- list(x = temp$x, y = temp$y[, best], lty = for_plot$lty[2],
                        lwd = for_plot$lwd[2], col = for_plot$col[2])
      do.call(graphics::lines, for_lines)
    }
  }
  return(invisible(for_plot))
}
