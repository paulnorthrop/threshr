# plot.ithreshpred
# summary.ithresh

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
#' @param n_years A numeric vector. Value(s) of N.  If \code{which_u = "all"}
#'   then \code{n_years} must have length one.
#' @param which_u Either a character scalar or a numeric scalar.
#'
#'   If \code{which_u} is a character scalar it must be either "best" or "all".
#'   If \code{which_u = "best"} then the threshold achieving the largest
#'   measure of predictive performance in \code{object$pred_perf}, based
#'   on the validation threshold selected using \code{which_v}, is used to
#'   perform prediction.
#'   If \code{which_u = "all"} then \emph{all} the thresholds are used to
#'   perform prediction.  The inferences from each threshold are weighted
#'   according to equation (15) in
#'   \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#'   based on the prior probabilities of thresholds in \code{u_prior}
#'   and column \code{which_v} of the measures of predictive performance in
#'   \code{object$pred_perf}.
#'
#'   Otherwise, \code{which_u} is a numeric scalar that indicates which
#'   element of \code{u_vec} the user wishes to select as a threshold on
#'   which to base prediction, that is, \code{which_u} must be an integer
#'   in {1, ..., \code{length(object$u_vec)}}.
#' @param which_v A numeric scalar. Indicates which element of
#'   \code{object$v_vec} is used in selecting a single threshold
#'   (if \code{which_u = "best"}) or weighting the inferences from
#'   all thresholds (if \code{which_u = "all"}).
#' @param u_prior  A numeric vector.  Prior probabilities for the thresholds
#'   in \code{u_prior}.  By default this is set to a vector of length
#'   \code{length(u_vec)} with each element equal to 1/\code{length(u_vec)}.
#'   \code{u_prior} is normalized to have sum equal to 1 inside
#'   \code{predict.ithresh}.
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
#' @return An object of class "ithreshpred" with a similar structure to
#'   an object of class "evpred" returned from
#'   \code{\link[revdbayes]{predict.evpost}}.
#'   In addition, the object contains \code{which_u} and \code{which_v}
#'   and \code{u_vec = object$u_vec} and \code{v_vec = object$v_vec}.
#' @examples
#' \dontrun{
#' # Gulf of Mexico significant wave heights, default priors.
#' u_vec <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
#' npy_gom <- length(gom)/105
#' gom_cv <- ithresh(data = gom, u_vec = u_vec, n_v = 4)
#' pjn <- predict(gom_cv, npy = npy_gom)
#' plot(pjn)
#' }
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{http://dx.doi.org/10.1111/rssc.12159}
#' @export
predict.ithresh <- function(object, npy = NULL, n_years = 100,
                            which_u = c("best", "all"), which_v = 1L,
                            u_prior = NULL, type = c("p", "d", "q", "i", "r"),
                            ...) {
  if (!inherits(object, "ithresh")) {
    stop("object must be of class ''ithresh'', produced by ithresh()")
  }
  # Numbers of training and validation thresholds
  n_u <- length(object$u_vec)
  n_v <- length(object$v_vec)
  # Check inputs
  type <- match.arg(type)
  if (is.character(which_u)) {
    which_u <- match.arg(which_u)
  } else if (is.numeric(which_u)) {
    if (!(which_u %in% 1:n_u)) {
      stop("'which_u' must be in 1:length(object$u_vec)")
    }
  } else {
    stop("'which_u' must be ''best'' or ''all'' or a numeric scalar")
  }
  if (!is.numeric(n_years)) {
    stop("'n_years' must be numeric")
  }
  if (!is.numeric(which_v) || !(which_v %in% 1:n_v)) {
    stop("'which_v' must be in 1:length(object$v_vec)")
  }
  # u_prior must be a numeric vector of length n_u with no negative entries
  p_bad <- !is.numeric(u_prior) || length(u_prior) != n_u || any(u_prior < 0)
  if (!is.null(u_prior)) {
    if (p_bad) {
      stop("'u_prior' must be a non-negative length(object$u_vec)-vector")
    } else {
      u_prior <- u_prior / sum(u_prior)
    }
  }
  # Check that npy has been supplied
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
  if (!is.numeric(npy)) {
    stop("'npy' must be numeric")
  }
  # Select the user's option based on which_u -----------
  if (which_u == "best" || is.numeric(which_u)) {
    # Best threshold
    if (which_u == "best") {
      which_u <- which.max(object$pred_perf[, which_v])
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
    for_predict_evpost <- list(object = evpost_obj, n_years = n_years,
                               npy = npy, type = type)
    ret_obj <- do.call(revdbayes:::predict.evpost, for_predict_evpost)
  }
  # Need to do this: type = "p" only, do for all thresholds
  # Calculate weights (u_prior and normalized pred_perf)
  # Calculate threshold-averaged
  # Create $x and $y as matrices: averaged column first
  if (which_u == "all") {
    # All thresholds
    which_u <- which.max(object$pred_perf[, which_v])
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
    ret_obj <- do.call(revdbayes:::predict.evpost, for_predict_evpost)
  }
  ret_obj$which_u <- which_u
  ret_obj$u_vec <- object$u_vec
  ret_obj$which_v <- which_v
  ret_obj$v_vec <- object$v_vec
  class(ret_obj) <- "ithreshpred"
  return(ret_obj)
}

