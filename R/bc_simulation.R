#' Simulation with (Inverse) Box-Cox Transformation
#'
#' Random generation for a randm variable that has the following property:
#' a Box-Cox transformation of the variable has a specified target
#' distribution.
#' @param n An integer scalar.  The sample size required.
#' @param lambda A numeric scalar.  The Box-Cox transformation parameter.
#' @param sim_fn An R function that simulates from the target probability
#'   distribution.
#' @param ... Further arguments to be passed to \code{sim_fn}.
#' @details A sample of size \code{n} is simulated using the function supplied
#'   in \code{sim_fn}.  Then the inverse of the Box-Cox transformation is
#'   applied to each value in this sample.
#' @return A numeric vector of length \code{n}.
#' @examples
#' x <- rboxcox(n = 100, lambda = 2)
#' x <- rboxcox(n = 100, lambda = 2, sim_fn = revdbayes::rgp)
#' x <- rboxcox(n = 100, lambda = 2, sim_fn = revdbayes::rgp, shape = 0.1)
#' @export
rboxcox <- function(n = 1, lambda = 1, sim_fn = stats::rexp, ...) {
  return(inv_bc(sim_fn(n, ...), lambda = lambda))
}
