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
  pjn <- revdbayes:::grimshaw_gp_mle(gp_data)
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

gp_obs_info <- function(gp_pars, y) {
  # Observed information for the generalized Pareto distribution
  #
  # Calculates the observed information matrix for a random sample \code{y}
  # from the generalized Pareto distribution, i.e. the negated Hessian matrix
  # of the generalized Pareto log-likelihood, evaluated at \code{gp_pars}.
  #
  # Args:
  #   gp_pars : A numeric vector. Parameters sigma and xi of the
  #   generalized Pareto distribution.
  #   y       : A numeric vector. A sample of positive data values.
  #
  # Returns:
  #   A 2 by 2 numeric matrix.  The observed information matrix.
  #
  s <- gp_pars[1]
  x <- gp_pars[2]
  i <- matrix(NA, 2, 2)
  i[1,1] <- -sum((1 - (1 + x) * y * (2 * s + x * y) / (s + x * y) ^ 2) / s ^ 2)
  i[1,2] <- i[2,1] <- -sum(y * (1 - y / s) / (1 + x * y / s) ^ 2 / s ^ 2)
  i[2,2] <- sum(2 * log(1 + x * y / s) / x ^ 3 - 2 * y / (s + x * y) / x ^ 2 -
                  (1 + 1 / x) * y ^ 2 / (s + x * y) ^ 2)
  return(i)
}
