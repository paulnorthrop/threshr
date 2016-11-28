# =========================== plot.thresh ===========================

#' Plot diagnostics for a thresh
#'
#' \code{plot} method for class "thresh".
#'
#' @param x an object of class "thresh", a result of a call to
#'   \code{\link{ithresh}}.
#' @param y Not used.
#' @param ... Additional arguments passed on to \code{matplot}.
#' @details Add some details.
#' @examples
#' # NAJ 2016
#' library(revdbayes)
#' data(gom)
#' u_vec <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
#' v_vec <- quantile(gom, probs = c(0.8, 0.85, 0.9, 0.95))
#' naj_control <- list(prior_args = list(max_xi = 1))
#' gom_naj <- ithresh(data = gom, method = "naj", u_vec = u_vec, v_vec = v_vec)
#' plot(gom_naj)
#' @export
plot.thresh <- function(x, y, ...) {
  if (!inherits(x, "thresh")) {
    stop("use only with \"thresh\" objects")
  }
  if (x$method == "naj") {
    print("2")
    x_data <- x$u_vec
    y_data <- x$tweights
  }
  graphics::matplot(x_data, y_data, type = "l", ...)
}
