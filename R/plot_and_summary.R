# =========================== plot.ithresh ===========================

#' Plot diagnostics an ithresh object
#'
#' \code{plot} method for class "ithresh".  Produces an extreme value
#' threshold diagnostic plot based on an analysis performed by
#' \code{\link{ithresh}}.
#'
#' @param x an object of class "ithresh", a result of a call to
#'   \code{\link{ithresh}}.
#' @param y Not used.
#' @param which_val A numeric vector specifying the validation thresholds, that
#'   is the indices of the argument \code{v_vec} to \code{\link{ithresh}}.
#' @param prob A logical scalar.  If \code{TRUE} then the levels of thresholds
#'   are represented by the proportion of observations that lie below a
#'   threshold.  If \code{prob = FALSE} then the values of the thresholds are
#'   used.
#' @param top_scale A logical scalar indicating Whether or not to add a scale
#'   to the top horizontal axis.  If this is added it gives the threshold on
#'   the scale not chosen by \code{prob}.
#' @param add_legend A logical scalar indicating whether or not to add a
#'   legend to the plot.  If \code{method = "cv"} then the legend gives the
#'   levels of the validation thresholds.
#' @param legend_pos The position of the legend (if required) specified using
#'   the argument \code{x} in \code{\link[graphics]{legend}}.
#' @param ... Additional arguments passed on to \code{\link[graphics]{matplot}}
#'   and/or \code{\link[graphics]{legend}} and/or \code{\link[graphics]{axis}}.
#' @details Produces plots of the \emph{threshold weights}, defined in
#'  equation (14) of
#'   \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)},
#'   against training threshold.  A line is produced for each of the validation
#'   thresholds chosen in \code{which_v}.  The result is a plot like those in
#'   the top row of Figure 7 in
#'   \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}.
#' @examples
#' \dontrun{
#' u_vec <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
#' gom_cv <- ithresh(data = gom, u_vec = u_vec, n_v = 4)
#' plot(gom_cv, lwd = 2, add_legend = TRUE, legend_pos = "topleft")
#' mtext("significant wave height / m", side = 3, line = 2.5)
#' }
#' @seealso \code{\link{ithresh}} for threshold selection in the i.i.d. case
#'   based on leave-one-out cross-validaton.
#' @seealso \code{\link{summary.ithresh}} Summarizing measures of threshold
#'   predictive performance.
#' @export
plot.ithresh <- function(x, y, ..., which_val = NULL, prob = TRUE, top_scale = TRUE,
                        add_legend = FALSE, legend_pos = "topleft") {
  if (!inherits(x, "ithresh")) {
    stop("use only with \"ithresh\" objects")
  }
  # Use only the validation thresholds in columns which_val.
  if (!is.null(which_val)) {
    x$pred_perf <- x$pred_perf[, which_val, drop = FALSE]
    x$v_ps <- x$v_ps[which_val]
    x$v_vec <- x$v_vec[which_val]
  }
  # Calculate threshold weights.  Shift to avoid underflow.
  n_u <- length(x$u_vec)
  n_v <- length(x$v_vec)
  shoof <- matrix(colMeans(x$pred_perf * !is.infinite(x$pred_perf),
                           na.rm = TRUE), ncol = n_v, nrow = n_u,
                  byrow = TRUE)
  y_data <- apply(exp(x$pred_perf - shoof), 2,
                  function(x) x / sum(x, na.rm =TRUE))
  y_lab <- "threshold weight"
  if (prob) {
    v_data <- x$v_ps
  } else {
    v_data <- x$v_vec
  }
  # General aspects.
  if (prob) {
    x_data <- x$u_ps
    t_data <- x$u_vec
    x_lab <- "quantile of training threshold / %"
  } else {
    x_data <- x$u_vec
    t_data <- x$u_ps
    x_lab <- "threshold"
  }
  xy_args <- list(x = x_data, y = y_data)
  # Look for user-supplied arguments to matplot.
  user_args <- list(...)
  m_cond <- names(user_args) %in% methods::formalArgs(graphics::matplot)
  a_cond <- names(user_args) %in% methods::formalArgs(graphics::axis)
  l_cond <- names(user_args) %in% methods::formalArgs(graphics::legend)
  legend_args <- user_args[l_cond]
  matplot_args <- user_args[(!l_cond & !a_cond) | m_cond]
  axis_args <- user_args[(!l_cond & !m_cond) | a_cond]
  axis_args$col <- 1
  if (is.null(matplot_args$xlab)) {
    matplot_args$xlab <- x_lab
  }
  if (is.null(matplot_args$ylab)) {
    matplot_args$ylab <- y_lab
  }
  if (is.null(matplot_args$type)) {
    matplot_args$type <- "l"
  }
  if (is.null(matplot_args$col)) {
    matplot_args$col <- 1
    legend_args$col <- 1
  }
  if (is.null(matplot_args$lty)) {
    matplot_args$lty <- 1:ncol(y_data)
    legend_args$lty <- 1:ncol(y_data)
  }
  if (is.null(legend_args$type)) {
    legend_args$title <- "highest threshold"
  }
  legend_args$x <- legend_pos
  all_args <- c(xy_args, matplot_args)
  do.call(graphics::matplot, c(all_args, axes = FALSE))
  axis_args$side <- 2
  do.call(graphics::axis, axis_args)
  axis_args$side <- 1
  if (prob) {
    axis_args$at <- unique(c(pretty(x_data), v_data))
  } else {
    axis_args$at <- pretty(x_data)
  }
  do.call(graphics::axis, axis_args)
  if (!is.null(axis_args$lwd)) {
    graphics::box(lwd = axis_args$lwd)
  } else {
    graphics::box()
  }
  # Add top scale?
  if (top_scale) {
    top_vals <- pretty(t_data)
    if (prob) {
      top_vals <- top_vals[top_vals < max(x$u_vec)]
      axis_args$side <- 3
      axis_args$at <- 100 * stats::ecdf(x$data)(top_vals)
      axis_args$labels <- top_vals
      do.call(graphics::axis, axis_args)
    } else {
      top_vals <- unique(c(top_vals, x$v_ps))
      axis_args$side <- 3
      axis_args$at <- stats::quantile(x$data, probs = top_vals / 100)
      axis_args$labels <- top_vals
      do.call(graphics::axis, axis_args)
    }
  }
  # Add a legend?
  if (add_legend){
    if (prob) {
      legend_args$legend <- paste(v_data, "%")
      do.call(graphics::legend, legend_args)
    } else {
      legend_args$legend <- signif(v_data, 2)
      do.call(graphics::legend, legend_args)
    }
  }
}

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
#' @seealso \code{\link{stability}}.
#' @examples
#' u_vec <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
#' gom_stab <- stability(data = gom, u_vec = u_vec)
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
}


# ============================== summary.ithresh ==============================

#' Summarizing measures of threshold predictive performance
#'
#' \code{summary} method for class "ithresh"
#'
#' @param object an object of class "ithresh", a result of a call to
#'   \code{ithresh}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @return Returns a numeric matrix with 3 columns and \code{n_v} rows,
#'   where \code{n_v} is an argument to \code{\link{ithresh}} that
#'   determines how many of the largest training thresholds are used
#'   a validation thresholds.  The columns contain:
#' \itemize{
#'   \item {column 1:} {the validation threshold v}
#'   \item {column 2:} {the best training threshold u judged using the
#'     validation threshold v}
#'   \item {column 3:} {the index of the vector \code{u_vec} of training
#'     thresholds to which the threshold in column2 corresponds}
#' }
#' @examples
#' \dontrun{
#' u_vec <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
#' gom_cv <- ithresh(data = gom, u_vec = u_vec, n_v = 4)
#' summary(gom_cv)
#' }
#' @seealso \code{\link{ithresh}} for threshold selection in the i.i.d. case
#'   based on leave-one-out cross-validaton.
#' @seealso \code{\link{plot.ithresh}} for the S3 plot method for objects of
#'   class \code{ithresh}.
#' @export
summary.ithresh <- function(object, ...) {
  if (!inherits(object, "ithresh")) {
    stop("use only with \"ithresh\" objects")
  }
  # Find the best training threshold for each validation threshold
  which_best <- apply(object$pred_perf, 2, which.max)
  u_best <- object$u_vec[which_best]
  res <- cbind(object$v_vec, u_best, which_best)
  rownames(res) <- 1:4
  colnames(res) <- c("v", "best u", "index of u_vec")
  return(res)
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
#' @param ave_only.  A logical scalar.  Only relevant if
#'   \code{\link{predict.ithresh}} was called with \code{which_u = "all"}.
#'   If \code{TRUE} then plot only
#'   a curve for the weighted average over multiple training thresholds.
#'   If \code{FALSE} then also plot a curve for each training threshold.
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
#'   different training thresholds.
#'
#' @examples
#' \dontrun{
#' u_vec <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
#' npy_gom <- length(gom)/105
#' gom_cv <- ithresh(data = gom, u_vec = u_vec, n_v = 4)
#'
#' ### Best training threshold based on the lowest validation threshold
#'
#' # Predictive distribution function
#' best_p <- predict(gom_cv, npy = npy_gom, n_years = c(100, 1000))
#' plot(best_p)
#'
#' # Predictive density function
#' best_d <- predict(gom_cv, npy = npy_gom, type = "d", n_years = c(100, 1000))
#' plot(best_d)
#'
#' ### All thresholds plus weighted average of inferences over all thresholds
#'
#' all_p <- predict(gom_cv, npy = npy_gom, which_u = "all")
#' plot(all_p)
#'
#' # All thresholds plus weighted average of inferences over all thresholds
#' all_d <- predict(gom_cv, npy = npy_gom, which_u = "all", type = "d")
#' plot(all_d)
#' }
#' @seealso \code{\link{ithresh}} for threshold selection in the i.i.d. case
#'   based on leave-one-out cross-validaton.
#' @seealso \code{\link{predict.ithresh}} for predictive inference for the
#'   largest value observed in N years.
#' @seealso \code{\link{summary.ithresh}} Summarizing measures of threshold
#'   predictive performance.
#' @export
plot.ithreshpred <- function(x, y, ..., ave_only = FALSE) {
  if (!inherits(x, "ithreshpred")) {
    stop("use only with \"ithreshpred\" objects")
  }
  # Single threshold
  if (x$which_u == "best" || is.numeric(x$which_u)) {
    temp <- x
    class(temp) <- "evpred"
    revdbayes:::plot.evpred(temp, ...)
  }
  # Multiple thresholds
  if (x$which_u == "all") {
    temp <- x
    class(temp) <- "evpred"
    n_col_y <- ncol(temp$y)
    # Create list of arguments for revdbayes' plot.evpred()
    # Set lty, lwd and col only if they have not be supplied by the user
    for_plot <- list(x = temp, ...)
    if (is.null(for_plot$lty)) {
      for_plot$lty <- 1
    }
    # If we only want one `averaged' line then make all the others disappear
    if (ave_only) {
      for_plot$lty <- c(1, rep(0, n_col_y - 1))
    }
    if (is.null(for_plot$lwd)) {
      for_plot$lwd <- 2
    }
    if (is.null(for_plot$col)) {
      my_col <- c(1, rep(grey(0.75), n_col_y - 1))
      for_plot$col <- my_col
    }
    if (is.null(for_plot$leg_text)) {
      if (ave_only) {
        for_plot$leg_text <- c("averaged")
      } else {
        for_plot$leg_text <- c("averaged", "single")
      }
    }
    do.call(revdbayes:::plot.evpred, for_plot)
    # Replot the `average' line, so that it appears on top
    for_lines <- list(x = temp$x, y = temp$y[, 1], lty = for_plot$lty[1],
                      lwd = for_plot$lwd[1], col = for_plot$col[1])
    do.call(graphics::lines, for_lines)
  }
}
