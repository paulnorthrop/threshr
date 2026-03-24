#' Internal threshr functions
#'
#' Internal threshr functions
#' @details
#' These functions are not intended to be called by the user.
#' @name threshr-internal
#' @keywords internal
NULL

# ================================ Box-Cox ====================================

#' @keywords internal
#' @rdname threshr-internal
bc <- function(x, lambda = 1, lambda_tol = 1 / 50, m = 4) {
  #
  # Computes the Box-Cox transformation of a vector x.  If lambda is very close
  # to zero then a first order Taylor series approximation is used.
  #
  # This function is vectorised with respect to both x and lambda.
  #
  # Args:
  #   x          : A numeric vector. (Non-negative) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric vector.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  #   m          : order of TS expansion (1 for linear in lambda, 2 for
  #                quadratic etc)
  # Returns:
  #   A numeric vector.  The transformed value (x^lambda - 1) / lambda
  #
  # Example
  # x <- rexp(100)
  # y <- bc(x, lambda = 2)
  # hist(y, prob = TRUE)
  if (any(x < 0)) {
    stop("Invalid x: x must be non-negative")
  }
  max_len <- max(length(x), length(lambda))
  x <- rep_len(x, max_len)
  la <- rep_len(lambda, max_len)
  lnx <- log(x)
  i <- 0:m
  retval <- ifelse(abs(la) > lambda_tol, (x ^ la - 1) / la,
              ifelse(la == 0, lnx,
                ifelse(is.infinite(x),
                  ifelse(la < 0, -1 / la, Inf),
                    ifelse(x == 0, ifelse(la > 0, -1 / la, -Inf),
                      sum(lnx ^ (i + 1) * la ^ i / factorial(i + 1))))))
  return(retval)
}

# ============================= Inverse Box-Cox ===============================

#' @keywords internal
#' @rdname threshr-internal
inv_bc <- function(x, lambda = 1, lambda_tol = 1 / 50, m = 4) {
  #
  # Computes the inverse Box-Cox transformation of a vector.  If lambda is very
  # close to zero then a first order Taylor series approximation is used.
  #
  # Args:
  #   x          : A numeric vector. (Non-negative) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  #   m          : order of TS expansion for the log of the inverse function
  #                (1 for linear in lambda, 2 for quadratic etc)
  # Returns:
  #   A numeric vector.  The inverse transformed value
  #     (1 + lambda x) ^ (1 / lambda)
  #
  # Example
  # x <- rexp(100)
  # y <- bc(x, lambda = 2)
  # x2 <- inv_bc(y, lambda = 2)
  # summary(x - x2)
  # hist(y, prob = TRUE)
  if (lambda == 0) {
    return(exp(x))
  }
  if (lambda > 0 && any(x <= -1 / lambda)) {
    stop("At least one component of x is too small")
  }
  if (lambda < 0 && any(x >= -1 / lambda)) {
    stop("At least one component of x is too large")
  }
  if (abs(lambda) > lambda_tol) {
    retval <- (1 + lambda * x) ^ (1 / lambda)
  } else {
    j <- 0:m
    fun <- function(x) {
      return(exp(x * sum((-1) ^ j * (lambda * x) ^ j / (j + 1))))
    }
    retval <- vapply(x, fun, 0.0)
  }
  return(retval)
}

# ================================ Yeo-Johnson ================================

#' @keywords internal
#' @rdname threshr-internal
yj <- function(x, lambda = 1, lambda_tol = 1 / 50, m = 4) {
  #
  # Computes the Yeo-Johnson transformation of a vector.  If lambda is very
  # close to zero then a first order Taylor series approximation is used.
  #
  # Args:
  #   x          : A numeric vector. Values to be Yeo-Johnson transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  #   m          : order of TS expansion (1 for linear in lambda, 2 for
  #                quadratic etc)
  # Returns:
  #   A numeric vector.  The transformed value
  #     bc(1 + x, lambda),     if x >= 0,
  #     bc(1 - x, 2 - lambda), if x < 0
  #   where bc() is the Box-Cox transformation
  #
  # Example
  # x <- rnorm(100)
  # y <- yj(x, lambda = 2)
  # hist(y, prob = TRUE)
  retval <- numeric(length(x))
  xlt0 <- x < 0
  xgt0 <- !xlt0
  retval[xlt0] <- -bc(1 - x[xlt0], 2 - lambda, lambda_tol, m)
  retval[xgt0] <- bc(1 + x[xgt0], lambda, lambda_tol, m)
  return(retval)
}

# ============================ Inverse Yeo-Johnson ============================

#' @keywords internal
#' @rdname threshr-internal
inv_yj <- function(x, lambda = 1, lambda_tol = 1 / 50, m = 4) {
  #
  # Computes the Yeo-Johnson transformation of a vector.  If lambda is very
  # close to zero then a first order Taylor series approximation is used.
  #
  # Args:
  #   x          : A numeric vector. Values to be Yeo-Johnson transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  #   m          : order of TS expansion (1 for linear in lambda, 2 for
  #                quadratic etc)
  # Returns:
  #   A numeric vector.  The transformed value
  #     inv_bc(1 + x, lambda),     if x >= 0,
  #     inv_bc(1 - x, 2 - lambda), if x < 0
  #   where bc() is the Box-Cox transformation
  #
  # Example
  # x <- rnorm(100)
  # y <- yj(x, lambda = 2)
  # x2 <- inv_yj(y, lambda = 2)
  # summary(x - x2)
  # hist(y, prob = TRUE)
  retval <- numeric(length(x))
  xlt0 <- x < 0
  xgt0 <- !xlt0
  retval[xlt0] <- 1 - inv_bc(-x[xlt0], 2 - lambda, lambda_tol, m)
  retval[xgt0] <- inv_bc(x[xgt0], lambda, lambda_tol, m) - 1
  return(retval)
}

# ============================= Used in gp_mle() ============================ #

#' @keywords internal
#' @rdname threshr-internal
diag_pos <- function(x) {
  # Diagonal elements of a matrix, returning NA for any negative values
  y <- diag(x)
  y[y < 0] <- NA
  return(y)
}

# =========================== gp_mle ===========================

gp_mle <- function(gp_data) {
  # Maximum likelihood estimation for the generalized Pareto distribution
  #
  # Performs maximum likelihood estimation for the generalized Pareto
  # distribution.  Uses the function \code{gpdmle} associated with
  # Grimshaw (1993), which returns MLEs of sigma and k = - \xi.
  #
  # Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates
  #   for the Generalized Pareto Distribution.  Technometrics, 35(2), 185-191.
  #   and Computing (1991) 1, 129-133. https://doi.org/10.1007/BF01889987.
  #
  # Args:
  #   gp_data : A numeric vector containing positive values, assumed to be a
  #             random sample from a generalized Pareto distribution.
  #
  # Returns:
  #   A list with components
  #     mle  : A numeric vector.  MLEs of GP parameters sigma and xi.
  #     nllh : A numeric scalar.  The negated log-likelihood at the MLE.
  #
  # Remove missing values
  gp_data <- as.numeric(stats::na.omit(gp_data))
  # Call Grimshaw (1993) function, note: k is -xi, a is sigma
  pjn <- revdbayes::grimshaw_gp_mle(gp_data)
  temp <- list()
  temp$mle <- c(pjn$a, -pjn$k)  # mle for (sigma,xi)
  sc <- rep(temp$mle[1], length(gp_data))
  xi <- temp$mle[2]
  temp$nllh <- sum(log(sc)) + sum(log(1 + xi * gp_data / sc) * (1 / xi + 1))
  gp_info <- gp_obs_info(temp$mle, gp_data)
  temp$cov <- solve(gp_info)
  temp$se <- sqrt(diag_pos(temp$cov))
  return(temp)
}

# =========================== gp_obs_info ===========================

gp_obs_info <- function(gp_pars, y, eps = 1e-5, m = 3) {
  # Observed information for the generalized Pareto distribution
  #
  # Calculates the observed information matrix for a random sample \code{y}
  # from the generalized Pareto distribution, i.e. the negated Hessian matrix
  # of the generalized Pareto log-likelihood, evaluated at \code{gp_pars}.
  #
  # Args:
  #   gp_pars : A numeric vector. Parameters sigma and xi of the
  #             generalized Pareto distribution.
  #   y       : A numeric vector. A sample of positive data values.
  #   eps     : A (small, positive) numeric scalar.  If abs(xi) is smaller than
  #             eps then we approximate the [2, 2] element of the information
  #             matrix using a Taylor series approximation.
  #   m       : A non-negative integer.  The order of the polynomial used to
  #             make this approximation.
  #
  # Returns:
  #   A 2 by 2 numeric matrix.  The observed information matrix.
  #
  if (eps <= 0) {
    stop("'eps' must be positive")
  }
  if (m < 0) {
    stop("'m' must be non-negative")
  }
  # sigma
  s <- gp_pars[1]
  # xi
  x <- gp_pars[2]
  i <- matrix(NA, 2, 2)
  i[1, 1] <- -sum((1 - (1 + x) * y * (2 * s + x * y) / (s + x * y) ^ 2) / s ^ 2)
  i[1, 2] <- i[2, 1] <- -sum(y * (1 - y / s) / (1 + x * y / s) ^ 2 / s ^ 2)
  # Direct calculation of i22 is unreliable for x close to zero.
  # If abs(x) < eps then we expand the problematic terms (all but t4 below)
  # in powers of xi up to xi ^ m. The terms in 1/z and 1/z^2 cancel.
  z <- x / s
  t0 <- 1 + z * y
  t4 <- y ^ 2 / t0 ^ 2
  if (any(t0 <= 0)) {
    stop("The log-likelihood is 0 for this combination of data and parameters")
  }
  if (abs(x) < eps) {
    j <- 0:m
    zy <- z * y
    sum_fn <- function(zy) {
      return(sum((-1) ^ j * (j ^ 2 + 3 * j + 2) * zy ^ j / (j + 3)))
    }
    tsum <- vapply(zy, sum_fn, 0.0)
    i[2, 2] <- sum(y ^ 3 * tsum / s ^ 3 - t4 / s ^ 2)
  } else {
    t1 <- 2 * log(t0) / z ^ 3
    t2 <- 2 * y / (z ^ 2 * t0)
    t3 <- y ^ 2 / (z * t0 ^ 2)
    i[2, 2] <- sum((t1 - t2 - t3) / s ^ 3 - t4 / s ^ 2)
  }
  dimnames(i) <- list(c("sigma[u]", "xi"), c("sigma[u]", "xi"))
  return(i)
}
