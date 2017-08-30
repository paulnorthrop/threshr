# =========================== stability ===========================
#
#' Generalised Pareto parameter estimate stability
#'
#' Uses maximum likelihood estimation to fit a Generalised Pareto (GP)
#' model to threshold excesses over a range of thresholds.
#' The threshold excesses are treated as independent and identically
#' distributed (i.i.d.) observations.
#' The resulting estimates and confidence intervals can be plotted,
#' using \code{\link{plot.stability}},
#' to produde a crude graphical diagnostic for threshold choice.
#'
#' @param data  A numeric vector of observations.
#' @param u_vec A numeric vector of thresholds to be applied to the data.
#'   Any duplicated values will be removed.  These could be set at sample
#'   quantiles of \code{data} using \code{\link[stats]{quantile}}.
#' @param prof A logical scalar.  Whether to calculate confidence intervals
#'   for the GP shape parameter \eqn{\xi} based on the profile-likelihood
#'   for \eqn{\xi} or using the MLE plus or minus a multiple of the estimated
#'   standard error (SE) of the MLE.  The intervals produced by the former
#'   may be better but they take longer to calculate.
#'   Default: \code{FALSE}.
#' @param conf A numeric scalar in (0, 100).  Confidence level for the
#'   confidence intervals.  Default: 95\%.
#' @param mult A numeric vector of length 2.  The range of values over
#'   which the profile log-likelihood for \eqn{\xi} is calculated is
#'   (MLE - \code{mult[1]} c SE, MLE + \code{mult[2]} c SE),
#'   where MLE and SE are the MLE and estimated standard error of \eqn{\xi}
#'   and c is the constant for which this interval gives
#'   an approximate 100\code{conf}\% level confidence interval for \eqn{\xi}
#'   when \code{mult = c(1, 1)}.  The default, \code{mult = c(1, 2)}, works
#'   well in most cases.  If the routine fails because the range of \eqn{\xi}
#'   is not sufficiently wide then the relevant components of \code{mult}
#'   should be increased.
#' @param plot_prof A logical scalar.  Only relevant if \code{prof = TRUE}.
#'   If \code{plot_prof = TRUE} then the profile log-likelihood is plotted
#'   for each threshold.  If \code{FALSE} then nothing is plotted.
#' @param ... Further (optional) arguments to be passed to the
#'   \code{\link[stats]{optim}} function for the optimizations on which
#'   the profile-likelihood for \eqn{xi} is based.
#' @details
#'   For each threshold in \code{u_vec} a GP model is fitted by maximum
#'   likelihood estimation to the threshold excesses, i.e. the amounts
#'   by which the data exceed that threshold.  The MLEs of the GP shape
#'   parameter $\eqn{\xi}$ and approximate \code{conf}\% confidence intervals
#'   for \eqn{\xi} are stored for plotting (by \code{\link{plot.stability}})
#'   to produce a simple graphical diagnostic to inform threshold selection.
#'   This plot is used to choose a threshold above which the underlying GP
#'   shape parameter may be approximately constant. See Chapter 4 of
#'   Coles (2001).  See also the vignette "Introducing threshr".
#' @return
#'   An object (list) of class "stability" with components:
#'     \item{ests}{MLEs of the GP shape parameter \eqn{\xi}.}
#'     \item{ses}{Estimated SEs of the MLEs of \eqn{\xi}.}
#'     \item{lower}{Lower limit of 100\code{conf}\% confidence intervals
#'      for \eqn{\xi}.}
#'     \item{upper}{Upper limit of 100\code{conf}\% confidence intervals
#'      for \eqn{\xi}.}
#'     \item{nexc}{The number of threshold excesses.}
#'     \item{u_vec}{The thresholds supplied by the user.}
#'     \item{u_ps}{The approximate sample quantiles to which the thresholds
#'      in \code{u_vec} correspond.}
#'     \item{data}{The input \code{data}.}
#'     \item{conf}{The input \code{conf}.}
#'   Each of these components is a numeric vector of length
#' \code{length(u_vec)}.
#' @seealso \code{\link{ithresh}} for threshold selection in the i.i.d. case
#'   based on leave-one-out cross-validaton.
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   \url{http://dx.doi.org/10.1007/978-1-4471-3675-0_3}
#' @seealso \code{\link{plot.stability}} for the S3 \code{plot} method for
#'   objects of class \code{stability}.
#' @seealso \code{\link[stats]{quantile}}.
#' @examples
#' # Set a vector of thresholds
#' u_vec_gom <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
#'
#' # Symmetric confidence intervals
#' gom_stab <- stability(data = gom, u_vec = u_vec_gom)
#' plot(gom_stab)
#'
#' \dontrun{
#' # Profile-likelihood-based confidence intervals
#' gom_stab <- stability(data = gom, u_vec = u_vec_gom, prof = TRUE)
#' plot(gom_stab)
#' }
#' @export
stability <- function(data, u_vec, prof = FALSE, conf = 95, mult = 1:2,
                      plot_prof = FALSE, ...){
  # Put thresholds in ascending order and remove any repeated values.
  u_vec <- unique(sort(u_vec))
  n_u <- length(u_vec)
  temp <- list()
  temp$ests <- temp$ses <- temp$upper <- temp$lower <- numeric(n_u)
  # Numbers of excesses of each threshold
  temp$nexc <- unlist(lapply(u_vec, function(x) sum(data > x)))
  temp$u_vec <- u_vec
  if (prof) {
    cat(paste("Fitting at threshold number ..."), fill = TRUE)
  }
  conf <- conf / 100
  for (i in 1:n_u) {
    if (prof) {
      cat(paste(i, ""))
    }
    z_u <- data[data > u_vec[i]] - u_vec[i]
    z <- gp_mle(gp_data = z_u)
    z$data <- data[data > u_vec[i]]
    z$threshold <- u_vec[i]
    temp$ests[i] <- z$mle[2]
    d <- matrix(c(1, -u_vec[i]), ncol = 1)
    v <- t(d) %*% z$cov %*% d
    temp$ses[i] <- z$se[2]
    z_level <- -stats::qnorm((1 - conf) / 2)
    temp$upper[i] <- temp$ests[i] + z_level * temp$ses[i]
    temp$lower[i] <- temp$ests[i] - z_level * temp$ses[i]
    if (prof){
      xlow <- temp$ests[i] - z_level * mult[1] * temp$ses[i]
      xup <- temp$ests[i] + z_level * mult[2] * temp$ses[i]
      prof_res <- gp_profxi(z, xlow = xlow, xup = xup, conf = conf,
                            plot_prof = plot_prof, thresh_number = i, ...)
      temp$upper[i] <- prof_res[1]
      temp$lower[i] <- prof_res[2]
    }
  }
  if (is.null(names(u_vec))) {
    temp$u_ps <- round(100 * sapply(u_vec, function(x) mean(data < x)))
  } else {
    temp$u_ps <- as.numeric(substr(names(u_vec), 1, nchar(names(u_vec),
                                                          type = "c") - 1))
  }
  temp$data <- data
  temp$conf <- conf
  class(temp) <- "stability"
  return(temp)
}

# =========================== gp_profxi ===========================

gp_profxi <- function (z, xlow, xup, conf = 0.95, nint = 100,
                       plot_prof = FALSE, thresh_number = NULL, ...){
  xdat <- z$data
  u <- z$threshold
  v1 <- numeric(nint)
  v2 <- numeric(nint)
  #
  gpd_plikxi <- function(a) {
    if (abs(xi) < 10 ^ (-4)) {
      if (a <= 0) {
        l <- 10 ^ 6
      } else {
        l <- length(xdat) * log(a) + sum(xdat - u) / a
      }
    }
    else {
      y <- (xdat - u) / a
      y <- 1 + xi * y
      if (any(y <= 0) || a <= 0) {
        l <- 10 ^ 6
      } else {
        l <- length(xdat) * log(a) + sum(log(y)) * (1 / xi + 1)
      }
    }
    l
  }
  # Check whether or not the use has supplied method for optim().
  user_dots <- list(...)
  ### Upper tail ...
  #
  x2 <- seq(z$mle[2], xup, length = nint)
  sol <- z$mle[1]
  for (i in 1:nint) {
    xi <- x2[i]
    if (is.null(user_dots$method)) {
      opt <- stats::optim(sol, gpd_plikxi, method = "BFGS", ...)
    } else {
      opt <- stats::optim(sol, gpd_plikxi, ...)
    }
    sol <- opt$par
    v2[i] <- opt$value
  }
  #
  ### Lower tail ...
  #
  x1 <- seq(z$mle[2], xlow, length = nint)
  sol <- z$mle[1]
  for (i in 1:nint) {
    xi <- x1[i]
    if (is.null(user_dots$method)) {
      opt <- stats::optim(sol, gpd_plikxi, method = "BFGS", ...)
    } else {
      opt <- stats::optim(sol, gpd_plikxi, ...)
    }
    sol <- opt$par
    v1[i] <- opt$value
  }
  x <- c(rev(x1), x2)
  v <- c(rev(v1), v2)
  ma <- -z$nllh
  if (plot_prof) {
    graphics::plot(x, -v, type = "l", xlab = expression(xi),
         ylab = "profile log-likelihood")
    if (!is.null(thresh_number)) {
      graphics::title(main = paste("threshold number", thresh_number))
    }
    graphics::abline(h = ma, col = 4)
    graphics::abline(h = ma - 0.5 * stats::qchisq(conf, 1), col = 4)
    u <- graphics::par("usr")								### extract plotting coords
  }
  yaxis <- -v
  xaxis <- x
  conf_line <- ma - 0.5 * stats::qchisq(conf, 1)
  temp <- diff(yaxis-conf_line > 0)			### to find where curve crosses CI line
  loc <- which(temp == -1)					### upper limit of CI
  if (length(loc) == 0) {
    stop("Please try again using an increased value of mult[2]")
  }
  x1 <- xaxis[loc]
  x2 <- xaxis[loc + 1]
  y1 <- yaxis[loc]
  y2 <- yaxis[loc + 1]
  up_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  loc <- which(temp == 1)
  if (length(loc) == 0) {
    stop("Please try again using an increased value of mult[1]")
  }
  x1 <- xaxis[loc]
  x2 <- xaxis[loc + 1]
  y1 <- yaxis[loc]
  y2 <- yaxis[loc + 1]
  low_lim <- x1 + (conf_line - y1) * (x2 - x1) / (y2 - y1)
  xi <- z$mle[3]
  if (plot_prof){
    graphics::abline(v = up_lim, lty = 2)
    graphics::text(up_lim, u[3] - 0.02 * (u[4] - u[3]), round(up_lim, 2),
                   xpd = TRUE, cex = 0.75)
    graphics::abline(v = low_lim, lty = 2)
    graphics::text(low_lim, u[3] - 0.02 * (u[4] - u[3]), round(low_lim, 2),
                   xpd = TRUE, cex = 0.75)
    graphics::abline(v = xi, lty = 2)
    graphics::text(xi, u[3] - 0.02 * (u[4] - u[3]), round(xi, 2), xpd = TRUE,
         cex = 0.75)
  }
  return(c(low_lim, up_lim))
}
