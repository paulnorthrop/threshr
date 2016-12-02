# =========================== plot.thresh ===========================

#' Plot diagnostics for a thresh
#'
#' \code{plot} method for class "thresh".
#'
#' @param x an object of class "thresh", a result of a call to
#'   \code{\link{ithresh}}.
#' @param y Not used.
#' @param prob A logical scalar.  If \code{TRUE} then the levels of thresholds
#'   are represented by the proportion of observations that lie below a
#'   threshold.  If \code{prob = FALSE} then the values of the thresholds are
#'   used.
#' @param top_scale A logical scalar indicating Whether or not to add a scale
#'   to the top horizontal axis.  If this is added if gives the threshold on
#'   the scale not chosen by \code{prob}.
#' @param add_legend A logical scalar indicating whether or not to add a
#'   legend to the plot.  If \code{method = "cv"} then the legend gives the
#'   levels of the validation thresholds.
#' @param legend_pos The position of the legend (if required) specified using
#'   the argument \code{x} in \code{\link[graphics]{legend}}.
#' @param ... Additional arguments passed on to \code{\link[graphics]{matplot}
#'   and/or \code{\link[graphics]{legend} and/or \code{\link[graphics]{axis}.
#' @details Add some details.
#' @examples
#' # cv 2016
#' library(revdbayes)
#' data(gom)
#' u_vec <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
#' v_vec <- quantile(gom, probs = c(0.8, 0.85, 0.9, 0.95))
#' cv_control <- list(prior_args = list(max_xi = 1))
#' gom_cv <- ithresh(data = gom, method = "cv", u_vec = u_vec, v_vec = v_vec)
#' plot(gom_cv)
#' @export
plot.thresh <- function(x, y, prob = TRUE, top_scale = TRUE, add_legend = FALSE,
                        legend_pos = "topleft", ...) {
  if (!inherits(x, "thresh")) {
    stop("use only with \"thresh\" objects")
  }
  # Aspects that are specific to the method.
  if (x$method == "cv") {
    y_data <- x$tweights
    y_lab <- "threshold weight"
    if (prob) {
      v_data <- x$v_ps
      x_lab <- "quantile of training threshold / %"
    } else {
      v_data <- x$v_vec
      x_lab <- "threshold"
    }
  }
  # General aspects.
  if (prob) {
    x_data <- x$u_ps
    t_data <- x$u_vec
  } else {
    x_data <- x$u_vec
    t_data <- x$u_ps
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
  if (x$method == "cv" & prob == TRUE) {
    do.call(graphics::matplot, c(all_args, axes = FALSE))
    axis_args$side <- 2
    do.call(graphics::axis, axis_args)
    axis_args$side <- 1
    axis_args$at <- unique(c(pretty(x$u_ps), x$v_ps))
    do.call(graphics::axis, axis_args)
    box()
  } else {
    do.call(graphics::matplot, all_args)
  }
  # Add top scale?
  if (top_scale) {
    top_vals <- pretty(t_data)
    if (prob) {
      top_vals <- top_vals[top_vals < max(x$u_vec)]
      axis_args$side <- 3
      axis_args$at <- 100 * ecdf(x$data)(top_vals)
      axis_args$labels <- top_vals
      do.call(graphics::axis, axis_args)
    } else {
      if (x$method == "cv") {
        top_vals <- unique(c(top_vals, x$v_ps))
      }
      axis_args$side <- 3
      axis_args$at <- quantile(x$data, probs = top_vals / 100)
      axis_args$labels <- top_vals
      do.call(graphics::axis, axis_args)
    }
  }
  # Add a legend?
  if (x$method == "cv" & add_legend){
    if (prob) {
      legend_args$legend <- paste(v_data, "%")
      do.call(graphics::legend, legend_args)
    } else {
      legend_args$legend <- signif(v_data, 2)
      do.call(graphics::legend, legend_args)
    }
  }
}
