# which_u = numeric vector (of indices) or "best" or all"
# which_v = number

# If which_u = numeric vector then use u_vec(which_u) [check which_u suitable]
# If which_u = "best" then use which_v to find single u
# If which_u = "all" use all and use which_v to find weighted average.

# Check which_u and which_v for admissibility.

# ============================== predict.ithresh ==============================

#' Predictive distribution function of the largest value observed in N years.
#'
#' \code{predict} method for class "ithresh".
#'
#' @param object An object of class "ithresh", a result of a call to
#'   \code{\link{ithresh}}.
#' @param npy A numeric scalar. The mean number of observations per year
#'   of data, after excluding any missing values, i.e. the number of
#'   non-missing observations divided by total number of years of non-missing
#'   data.
#' @param n_years A numeric vector. Values of N.
#' @param which_u Either a character scalar or a numeric vector.
#'
#'   If \code{which_u} is a character scalar it must be either "best" or "all".
#'   If \code{which_u = "best"} then the threshold achieving the largest
#'   measure of predictive performance in \code{object$pred_perf}, based
#'   on the validation threshold selected using \code{n_v}, is used to
#'   perform prediction.
#'   If \code{which_u = "all"} then \emph{all} the thresholds are used to
#'   perform prediction.  The inferences from each threshold are weighted
#'   according to equation (15) in
#'   \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#'   based on the prior probabilities of thresholds in \code{u_prior}
#'   and column \code{n_v} of the measures of predictive performance in
#'   \code{object$pred_perf}.
#'
#'   Otherwise, \code{which_u} is a numeric vector that indicates which
#'   elements of \code{u_vec} the user wishes to select as thresholds on
#'   which to base prediction, that is, it must contain integers in
#'   {1, ..., \code{length(object$u_vec)}}.
#' @param which_v A numeric scalar. Indicates which element of
#'   \code{object$v_vec} is used in selecting a single threshold
#'   (if \code{which_u = "best"}) or weighting the inferences from
#'   all thresholds (if \code{which_u = "all"}).
#' @param u_prior  A numeric vector.  Prior probabilities for the thresholds
#'   in \code{u_prior}.  By default this is set to a vector of length
#'   \code{length(u_vec)} with each element equal to 1/\code{length(u_vec)}.
#' @param type A character vector.  Indicates which type of inference is
#'   required:
#' \itemize{
#'   \item "p" for the predictive distribution function,
#'   \item "d" for the predictive density function,
#'   \item "q" for the predictive quantile function,
#'   \item "i" for predictive intervals,
#'   \item "r" for random generation from the predictive distribution.
#' }
#' @param ... Additional arguments passed on to
#'   \code{\link[revdbayes]{predict.evpost}}, in particular \code{x}.
#' @details Add details
#' @return An object of class "pred_ithresh"
#' @examples
#' # Gulf of Mexico significant wave heights, default priors.
#' u_vec <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
#' npy_gom <- length(gom)/105
#' gom_cv <- ithresh(data = gom, u_vec = u_vec, n_v = 4)
#' pjn <- predict(gom_cv, npy = npy_gom)
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{http://dx.doi.org/10.1111/rssc.12159}
#' @export
predict.ithresh <- function(object, npy = NULL, n_years = 100,
                            which_u = "best", which_v = 1, u_prior = NULL,
                            type = c("p", "d", "q", "i", "r"), ...) {
  if (!inherits(object, "ithresh")) {
    stop("object must be of class ''ithresh'', produced by ithresh()")
  }
  type <- match.arg(type)
  if (!is.null(object$npy) & !is.null(npy)) {
    warning(paste("Two values of npy supplied.  The value npy = ", npy,
                  " from the current call has been used.", sep=""))
  }
  if (!is.null(object$npy) & is.null(npy)) {
    npy <- object$npy
  }
  if (is.null(object$npy) & is.null(npy)) {
    stop("npy must be given, here or in call to ithresh.")
  }
  # Select the threshold
  if (which_u == "best") {
    which_u <- which.max(gom_cv$pred_perf[, which_v])
  }
  # Create a list object of class "evpost" for revdbayes::predict.evpost().
  evpost_obj <- list()
  class(evpost_obj) <- "evpost"
  # Extract posterior sample for this threshold
  n <- object$n
  which_rows <- (1 + (which_u - 1) * n):(which_u * n)
  evpost_obj$bin_sim_vals <- object$sim_vals[which_rows, 1]
  evpost_obj$sim_vals <- object$sim_vals[which_rows, 2:3]
  evpost_obj$model <- "bingp"
  evpost_obj$thresh <- object$u_vec[which_u]
  for_predict_evpost <- list(object = evpost_obj, n_years = n_years, npy = npy,
                             type = type)
  ret_obj <- do.call(predict, for_predict_evpost)
  class(ret_obj) <- "pred_ithresh"
  return(ret_obj)
}

