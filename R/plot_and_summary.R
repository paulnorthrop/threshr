# =========================== plot.stability ===========================

#' Plot diagnostics for a stability object
#'
#' \code{plot} method for objects of class "stability" returned from
#' \code{\link{stability}}
#'
#' @param x an object of class "stability", a result of a call to
#'   \code{\link{stability}}.
#' @param y Not used.
#' @param prob A logical scalar.  If \code{TRUE} then the levels of thresholds
#'   on the lower horizontal axis are represented by the proportion of
#'   observations that lie below a threshold.  If \code{prob = FALSE} then the
#'   values of the thresholds are used.
#' @param top_scale A character scalar.
#'   If \code{top_scale = "none"} then no axis labels appear on the upper
#'   horizontal axis.
#'   If \code{top_scale = "excesses"} then the number of threshold excesses
#'   at each threshold are indicated.
#'   If \code{top_scale = "opposite"} then the type of threshold level
#'   \emph{not} chosen using \code{prob} is indicated.
#' @param vertical A logical scalar.  Should the confidence intervals be
#'   depicted using a vertical line for each threshold (\code{TRUE}) or by
#'   joining up confidence limits across thresholds (\code{FALSE})?
#' @param ... Additional arguments passed on to
#'   \code{\link[graphics]{matplot}}, \code{\link[graphics]{axis}}
#'   and/or \code{\link[graphics]{segments}}.
#' @details Produces a simple threshold diagnostic plot based on the object
#'   returned from \code{\link{stability}}.
#'   The MLEs of the GP shape parameter $\eqn{\xi}$ and
#'   approximate \code{conf}\% confidence intervals
#'   for \eqn{\xi} are plotted against the threshold used to fit the GP model.
#'   This plot is used to choose a threshold above which the underlying GP
#'   shape parameter may be approximately constant. See Chapter 4 of
#'   Coles (2001).  See also the vignette "Introducing threshr".
#'   as described in .
#'   See also the vignette "Introducing threshr".
#' @return In addition to producing the plot a list of the arguments used
#'   by \code{\link[graphics]{matplot}}, \code{\link[graphics]{axis}} is
#'   returned (invisibly).
#' @seealso \code{\link{stability}}.
#' @examples
#' u_vec_gom <- quantile(gom, probs = seq(0, 0.9, by = 0.05))
#' gom_stab <- stability(data = gom, u_vec = u_vec_gom)
#' plot(gom_stab)
#' @export
plot.stability <- function(x, y, ..., prob = TRUE,
                           top_scale = c("none", "excesses", "opposite"),
                           vertical = TRUE) {
  if (!inherits(x, "stability")) {
    stop("use only with \"stability\" objects")
  }
  top_scale <- match.arg(top_scale)
  y_data <- cbind(x$lower, x$ests, x$upper)
  y_lab <- expression(xi)
  if (prob) {
    x_data <- x$u_ps
    t_data <- x$u_vec
    x_lab <- "quantile of training threshold / %"
  } else {
    x_data <- x$u_vec
    t_data <- x$u_ps
    x_lab <- "threshold"
  }
  if (top_scale == "excesses") {
    t_data <- x$nexc
  }
  xy_args <- list(x = x_data, y = y_data)
  # Look for user-supplied arguments to matplot.
  user_args <- list(...)
  m_cond <- names(user_args) %in% methods::formalArgs(graphics::matplot)
  a_cond <- names(user_args) %in% methods::formalArgs(graphics::axis)
  s_cond <- names(user_args) %in% methods::formalArgs(graphics::segments)
  matplot_args <- user_args[!a_cond | m_cond]
  axis_args <- user_args[!m_cond | a_cond]
  segments_args <- user_args[s_cond]
  axis_args$col <- 1
  if (is.null(matplot_args$xlab)) {
    matplot_args$xlab <- x_lab
  }
  if (is.null(matplot_args$ylab)) {
    matplot_args$ylab <- y_lab
  }
  if (is.null(matplot_args$type)) {
    matplot_args$type <- c("l", "b", "l")
  }
  if (is.null(matplot_args$pch)) {
    matplot_args$pch <- c(0, 16, 0)
  }
  if (is.null(matplot_args$col)) {
    matplot_args$col <- 1
  }
  if (is.null(matplot_args$lty)) {
    if (vertical) {
      matplot_args$lty <- c(0, 1, 0)
    } else {
      matplot_args$lty <- c(2, 1, 2)
    }
  }
  all_args <- c(xy_args, matplot_args)
  do.call(graphics::matplot, c(all_args, axes = FALSE))
  if (vertical) {
    segments_args$x0 <- x_data
    segments_args$x1 <- x_data
    segments_args$y0 <- x$lower
    segments_args$y1 <- x$upper
    do.call(graphics::segments, segments_args)
  }
  axis_args$side <- 2
  do.call(graphics::axis, axis_args)
  axis_args$side <- 1
  axis_args$at <- pretty(x_data)
  do.call(graphics::axis, axis_args)
  if (!is.null(axis_args$lwd)) {
    graphics::box(lwd = axis_args$lwd)
  } else {
    graphics::box()
  }
  # Add top scale?
  if (top_scale != "none") {
    axis_args$side <- 3
    if (top_scale == "excesses") {
      axis_args$labels <- t_data
      axis_args$at <- x_data
    } else {
      if (prob) {
        x_pos <- c(pretty(x_data), max(x_data))
        axis_args$at <- x_pos
        axis_args$labels <- signif(stats::quantile(x$data, probs = x_pos
                                                   / 100), 2)
      } else {
        top_vals <- pretty(t_data)
        top_vals <- unique(c(top_vals, max(x$u_ps)))
        axis_args$at <- stats::quantile(x$data, probs = top_vals / 100)
        axis_args$labels <- top_vals
      }
    }
    do.call(graphics::axis, axis_args)
  }
  return(invisible(list(matplot_args = matplot_args, axis_args = axis_args)))
}

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
