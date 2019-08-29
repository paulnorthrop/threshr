#' Extract results for a chosen value of \eqn{\lambda}
#'
#' Extract from an object of class \code{"bcthresh"} returned from
#' \code{\link{bcthresh}} measures of predictive performance and
#' posterior samples obtained after a Box-Cox transformation of the data
#' using a particular value of \eqn{\lambda}.
#'
#' @param x An object of class \code{"bcthresh"} returned from
#'   \code{\link{bcthresh}}.
#' @param lambda A numeric scalar.  Must be contained in \code{x$lambda}.
#' @details The function simply extracts the relevant objects from \code{x}
#'   and recreates the transformed data, training and validation thresholds
#'   using the values of \code{lambda}.  The returned object has class
#'   \code{"ithresh"} and the same structure as an object returned from
#'   \code{\link{ithresh}}, for which \code{plot}, \code{summary} and
#'   \code{predict} methods are available.
#' @return An object of class \code{"ithresh"}.  See \code{\link{ithresh}}.
#'   In addition, the input \code{lambda} is added to this returned list.
#' @seealso \code{\link{plot.ithresh}} for the S3 plot method for objects of
#'   class \code{ithresh}.
#' @seealso \code{\link{summary.ithresh}} Summarizing measures of threshold
#'   predictive performance.
#' @seealso \code{\link{predict.ithresh}} for predictive inference for the
#'   largest value observed in N years.
#' @examples
#' probs <- seq(0, 0.9, 0.1)
#' lambda <- seq(1, 3, 0.5)
#'
#' prior_args <- list(prior = "flatflat", bin_prior = "haldane",
#'                    h_prior = list(min_xi = -Inf))
#' gom_args <- list(data = gom, probs = probs, lambda = lambda)
#' gom_lambda <- do.call(bcthresh, c(gom_args, prior_args))
#'
#' # lambda = 1 (the default)
#' res1 <- choose_lambda(gom_lambda)
#' plot(res1)
#' # lambda = 2
#' res2 <- choose_lambda(gom_lambda, lambda = 2)
#' plot(res2)
#' @export
choose_lambda <- function(x, lambda = 1) {
  if (!inherits(x, "bcthresh")) {
    stop("use only with \"bcthresh\" objects")
  }
  if (!(lambda %in% x$lambda)) {
    stop(paste("lambda must be one of ", x$lambda))
  }
  # Find the index of lambda in x$lambda
  which_lambda <- which(lambda == x$lambda)
  # Extract the relevant predictive performances and the posterior sample
  x$pred_perf <- as.matrix(x$pred_perf[, , which_lambda])
  x$sim_vals <- as.matrix(x$sim_vals[, , which_lambda])
  # Create the correct data, training and validation thresholds
  x$data <- bc(x$data, lambda = lambda)
  x$u_vec <- bc(x$u_vec, lambda = lambda)
  x$v_vec <- bc(x$v_vec, lambda = lambda)
  x$lambda <- lambda
  class(x) <- "ithresh"
  return(x)
}
