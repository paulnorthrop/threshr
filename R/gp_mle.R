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

gp_obs_info <- function(gp_pars, y, eps = 1e-5) {
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
  #
  # Returns:
  #   A 2 by 2 numeric matrix.  The observed information matrix.
  #
  if (eps <= 0) {
    stop("'eps' must be positive")
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
  # in powers of z up to z ^ 2. The terms in 1/z and 1/z^2 cancel leaving
  # only a quadratic in z.
  z <- x / s
  t0 <- 1 + z * y
  if (any(t0 <= 0)) {
    stop("The log-likelihood is 0 for this combination of data and parameters")
  }
  if (abs(x) < eps) {
    s1 <- 12 * z ^ 2 * y ^ 2 / 5
    s2 <- 3 * z * y / 2
    s3 <- 2 / 3
    i[2, 2] <- sum(y ^ 3 * (s1 - s2 + s3) / s ^ 3 - (y / s) ^ 2 / t0 ^ 2)
  } else {
    t1 <- 2 * log(t0) / z ^ 3
    t2 <- 2 * y / (z ^ 2 * t0)
    t3 <- y ^ 2 / (z * t0 ^ 2)
    t4 <- y ^ 2 / t0 ^ 2
    i[2, 2] <- sum((t1 - t2 - t3) / s ^ 3 - t4 / s ^ 2)
  }
  dimnames(i) <- list(c("sigma[u]", "xi"), c("sigma[u]", "xi"))
  return(i)
}
