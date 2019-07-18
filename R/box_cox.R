# =========================== box_cox_gm ===========================

bc_gm <- function(x, lambda = 1, lngm = 0,
                  lambda_tol = .Machine$double.eps) {
  #
  # Computes the Box-Cox transformation of a vector.  If lambda is very close
  # to zero then a first order Taylor series approximation is used.
  #
  # Args:
  #   x          : A numeric vector. (Non-negative) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lngm       : standardisation constant, often the mean of the logs of the
  #                original data vector.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda.tol use
  #                a Taylor series expansion.
  # Returns:
  #   A numeric vector.  The transformed value
  #     (x^lambda - 1) / (lambda * gm ^ (lambda - 1),
  #   where gm = exp(lngm).
  #
  # Example
  # x <- rexp(100)
  # lngm <- mean(log(x))
  # y <- box_cox_gm(x, lambda = 2, lngm = lngm)
  # hist(y, prob = TRUE)
  if (any(x < 0)) {
    stop("Invalid x: x must be non-negative")
  }
  max_len <- max(length(x), length(lambda))
  x <- rep_len(x, max_len)
  lambda <- rep_len(lambda, max_len)
  retval <- ifelse(abs(lambda) > lambda_tol, (x ^ lambda - 1) / lambda,
                   ifelse(lambda == 0, log(x),
                          ifelse(is.infinite(x),
                                 ifelse(lambda < 0, -1 / lambda, Inf),
                                 ifelse(x == 0, ifelse(lambda > 0, -1 / lambda,
                                                       -Inf),
                                        log(x) * (1 + lambda / 2)))))
  retval <- retval * exp((1 - lambda) * lngm)
  return(retval)
}
