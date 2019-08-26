#' Internal threshr functions
#'
#' Internal threshr functions
#' @details
#' These functions are not intended to be called by the user.
#' @name threshr-internal
#' @keywords internal
NULL

# =========================== box_cox_gm ===========================

#' @keywords internal
#' @rdname threshr-internal
bc_gm <- function(x, lambda = 1, lngm = 0, lambda_tol = 1 / 50, m = 4) {
  #
  # Computes the Box-Cox transformation of a vector.  If lambda is very close
  # to zero then a first order Taylor series approximation is used.
  #
  # Args:
  #   x          : A numeric vector. (Non-negative) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric vector.  Transformation parameter.
  #   lngm       : standardisation constant, often the mean of the logs of the
  #                original data vector.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  #   m          : order of TS expansion (1 for linear, 1 for quadratic etc)
  # Returns:
  #   A numeric vector.  The transformed value
  #     (x^lambda - 1) / (lambda * gm ^ (lambda - 1),
  #   where gm = exp(lngm).
  #
  # Example
  # x <- rexp(100)
  # lngm <- mean(log(x))
  # y <- bc_gm(x, lambda = 2, lngm = lngm)
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
  retval <- retval * exp((1 - lambda) * lngm)
  return(retval)
}
