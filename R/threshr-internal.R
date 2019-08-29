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
