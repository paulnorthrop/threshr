# Look for names of u_vec and v_vec
# If there store them.
# Otherwise, calculate the non-exceedance probs and store them (0, 100)
# Add to returned object.
# Use in plot.thresh()
# Also xlab, ylab etc.
# In plot.thresh, options to have quantile on lower axis, value on upper axis

# =========================== thresh ===========================
#
#' Threshold selection in the i.i.d. case (peaks over threshold).
#'
#' Produces diagnostic plots to assist in the selection of an extreme value
#' threshold in the case where the data can be treated as independent and
#' identically distributed (i.i.d.) observations.  For example, it could be
#' that these observations are the cluster maxima resulting from the
#' declustering of time series data.  The diagnostic plots are based on
#' making inferences from a Generalised Pareto (GP) model for threshold
#' excesses.
#' @param data  A numeric vector of observations.
#' @param method A character scalar.  The method to be used.
#'   \itemize{
#'     \item {"stability"} {A parameter estimate stability plot for the GP
#'       shape parameter \eqn{\xi}}
#'     \item {"naj"} {Northrop, Attalides and Jonathan (2016)}
#'     \item {"nc"} {Northrop and Coleman (2014)}
#'     \item {"wadsworth"} {Wadsworth (2015)}
#'   }
#' @param u_vec A numeric vector. A vector of \emph{training} thresholds
#'   at which inferences are made from the GP model.
#' @param v_vec A numeric vector. Relevant only if \code{method = "naj"}.
#'   A vector of \emph{validation} thresholds used to quantify the predictive
#'   performance of the GP models fitted at the thresholds in \code{u_vec}.
#' @param naj_control A list of (optional) arguments for use if
#'   \code{method = "naj"}.  In particular:
#' \itemize{
#'   \item {\code{prior}} {A prior for the GP parameters, set using
#'     \code{\link[revdbayes]{set_prior}}.  Default: \code{prior = "flat"}.
#'     with \code{min_xi = -1}.}
#'   \item {\code{h_prior}} {A list of further arguments (hyperparameters)
#'     for the GP prior specified in \code{prior}.}
#'   \item {\code{bin_prior}} {A prior for the threshold exceedance
#'     probability \eqn{p}, set using \code{\link[revdbayes]{set_bin_prior}}.
#'     Default: \code{prior = "jeffreys"}.}
#'   \item {\code{h_bin_prior}} {A list of further arguments (hyperparameters)
#'     for the binomial prior specified in \code{bin_prior}.}
#'   \item {\code{n}} {The size of the posterior sample used to perform
#'     predictive inference.  Default: \code{n = 1000}.}
#'   \item {...} {Further arguments to \code{prior} or \code{bin_prior}.}
#' }
#' @param ... Further argments to be passed to
#' @details Add some details.
#' \emph{parameter estimate stability plot}:
#' \emph{Northrop, Attalides and Jonathan (2016)}:
#' \emph{Northrop and Coleman (2014)}:
#' \emph{Wadsworth (2015)}:
#' See the threshr vignette for further details and examples.
#' @return An object (list) of class \code{"thresh"}.
#' @seealso \code{\link[revdbayes]{rpost}} in the
#'   \code{\link[revdbayes]{revdbayes}} package for details of the arguments
#'   that can be passed to \code{rpost}.
#' @seealso \code{\link[revdbayes]{set_prior}} in the
#'   \code{\link[revdbayes]{revdbayes}} package for details of how to set a
#'   prior distribution for GP parameters.
#' @seealso \code{\link[revdbayes]{set_bin_prior}} in the
#'   \code{\link[revdbayes]{revdbayes}} package for details of how to set a
#'   prior distribution for the exceedance probability \eqn{p}.
#' @examples
#' # NAJ 2016
#' library(revdbayes)
#' data(gom)
#' u_vec <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
#' v_vec <- quantile(gom, probs = c(0.8, 0.85, 0.9, 0.95))
#' naj_control <- list(prior_args = list(max_xi = 1))
#' gom_naj <- ithresh(data = gom, method = "naj", u_vec = u_vec, v_vec = v_vec)
#' @export
ithresh <- function(data, method = c("stability", "naj", "nc", "wadsworth"),
                   u_vec, v_vec = max(u_vec), naj_control = list(), ...) {
  method <- match.arg(method)
  temp <- naj_fn(data = data, u_vec = u_vec, v_vec = v_vec,
                 naj_control = naj_control, ...)
  class(temp) <- "thresh"
  return(temp)
}

# =========================== naj_fn ===========================

naj_fn <- function(data, u_vec, v_vec = max(u_vec), naj_control = list(),
                   ...) {
  # Extract arguments for passing to revdbayes function rpost -----------------
  # GP prior.
  if (is.null(naj_control$prior)) {
    naj_control$prior <- "flat"
  }
  if (is.null(naj_control$h_prior$min_xi)) {
    naj_control$h_prior$min_xi <- -1
  }
  for_set_prior <- c(list(prior = naj_control$prior, model = "gp"),
                     naj_control$h_prior)
  gp_prior <- do.call(set_prior, for_set_prior)
  # Binomial prior.
  if (is.null(naj_control$bin_prior)) {
    naj_control$bin_prior <- "jeffreys"
  }
  for_set_bin_prior <- c(list(prior = naj_control$bin_prior),
                         naj_control$h_bin_prior)
  bin_prior <- do.call(set_bin_prior, for_set_bin_prior)
  # Size of samples from posterior distributions.
  if (is.null(naj_control$n)) {
    n <- 1000
  } else {
    n <- naj_control$n
  }
  # Set up quantities for passing t
  # Put thresholds in ascending order ...
  u_vec <- sort(u_vec)
  v_vec <- sort(v_vec)
  # Validation tresholds must be above the lowest training threshold ...
  if (min(v_vec) < min(u_vec)) {
    warning("validation thresholds below min(u_vec) have been ignored",
            immediate. = TRUE)
    v_vec <- v_vec[v_vec >= min(u_vec)]
  }
  # Set the numbers of training and validation thresholds.
  n_u <- length(u_vec)
  n_v <- length(v_vec)
  # Locate the largest observation(s) in the data.
  j_max <- which(data == max(data))
  data_max <- data[j_max]
  # Remove (one of the) the largest observation(s)
  data_rm <- data[-j_max[1]]
  #
  # A matrix to store the measures of predictive performance.
  pred <- matrix(NA, nrow = n_u, ncol = n_v)
    # Loop over the training thresholds -------------------------------------
  for (i in 1:n_u){
    # Simulate from (full) bin-GP posterior.
    u <- u_vec[i]
    temp <- rpost(n = n, model = "bingp", data = data, prior = gp_prior,
                  thresh = u, bin_prior = bin_prior)
    theta <- cbind(temp$bin_sim_vals, temp$sim_vals)
    # Simulate from the bin-GP posterior after removal of the maximum value.
    temp_rm <- rpost(n = n, model = "bingp", data = data, prior = gp_prior,
                     thresh = u, bin_prior = bin_prior)
    theta_rm <- cbind(temp_rm$bin_sim_vals, temp_rm$sim_vals)
    # Carry out leave-one-out Bayesian cross-validation -------------------
    # Which validation thresholds are >= u ?
    which_v <- which(v_vec >= u - 1e-10)
    v_vals <- v_vec[which_v]           # select the validation thresholds
    pred[i, which_v] <- bloocv(z = data, theta = theta, theta.rm = theta_rm,
                                  u1 = u, u2.vec = v_vals, z.max = data_max,
                                  z.rm = data_rm)
  }
  # Calculate threshold weights.  Shift to avoid underflow.
  shoof <- matrix(colMeans(pred*!is.infinite(pred), na.rm = TRUE),
                  ncol = n_v, nrow = n_u, byrow = TRUE)
  tweights <- apply(exp(pred - shoof), 2, function(x) x / sum(x, na.rm =TRUE))
  temp <- list(tweights = tweights, u_vec = u_vec, v_vec = v_vec,
               method = "naj")
  return(temp)

}

# =========================== bloocv ===========================

bloocv <- function(z, theta, theta.rm, u1, u2.vec, z.max, z.rm){
  #
  N1 <- sum(z.rm<=u1)             # number of values that do not exceed u1
  z.gt.u1 <- z.rm[z.rm>u1]        # data greater than u1
  p1 <- theta[,1]                 # p.u from posterior sample
  temp <- my.stuff(u1=u1,theta=theta,u2=u2.vec,z.gt.u1=z.gt.u1)
  t.col <- ncol(temp)             # ... so that we can extract the correct columns
  n2 <- length(u2.vec)            # number of validation thresholds considered
  p2 <- as.matrix(temp[,(t.col-n2+1):t.col])
  #
  #----------------------- Function pred.u2 ---------------------------#
  pred.u2 <- function(k){
    # Contribution for each z that is less than or equal to u1
    t1 <- 1
    if (N1>0) t1 <- mean((1-p2[,k])/(1-p1))/mean(1/(1-p1))
    # Contributions from zs that are in (u1,u2]
    t2 <- 1
    if (min(z.gt.u1)<=u2.vec[k]){
      W2 <- which(z.gt.u1<=u2.vec[k])
      fr2 <- as.matrix(temp[,W2])
      t2 <- colMeans((1-p2[,k])/fr2/p1)/colMeans(1/fr2/p1)
    }
    # Contributions from zs that are greater than u2
    t3 <- 1
    if (max(z.gt.u1)>u2.vec[k]){
      W3 <- which(z.gt.u1>u2.vec[k])
      fr3 <- as.matrix(temp[,W3])
      t3 <- 1/colMeans(1/fr3/p1)
    }
    # Contribution from z.max (NB. it is possible that z.max < u2)
    p1 <- theta.rm[,1]
    temp.rm <- my.stuff(u1=u1,theta=theta.rm,u2=u2.vec,z.max=z.max[1],for.max=T)
    p2 <- as.matrix(temp.rm[,-1])
    fn <- temp.rm[,1] # bin-GP density at max(z)
    t4 <- ifelse(z.max[1]>u2.vec[k],mean(p1*fn),mean(1-p2[,k]))
    n.max <- length(z.max)          # cardinality of max(z)
    #
    N1*log(t1)+sum(log(t2))+sum(log(t3))+n.max*log(t4)
  }
  #--------------------- End of function pred.u2 ----------------------#
  #
  sapply(1:n2,pred.u2)
  #
}

#################################################################################

my.stuff <- function(u1,theta,u2=u1,z.max=NULL,z.gt.u1=NULL,for.max=F){
  #
  # Quantities required to estimate f(y_r | y_(r))
  # u1    : training threshold
  # u2    : validation threshold (u1 <= u2) (Default: u2=u1)
  # theta : (m by 3) matrix: (p.u.1, sigma.1, xi)
  #
  # Posterior samples of (p,sigma,xi) at training threshold u1 ...
  p.1 <- theta[,1]; sigma.1 <- theta[,2]; xi <- theta[,3]
  # Probability of exceeding validation threshold u2, for each u2
  # (Note: this could be zero)
  p.2.vec <- sapply(u2,function(xx) p.1*apply(cbind(1+xi*(xx-u1)/sigma.1,0),1,max)^(-1/xi))
  # bin-GP density at max(z)
  fGP.n <- ifelse(1+xi*(z.max-u1)/sigma.1>0,(1+xi*(z.max-u1)/sigma.1)^(-(1+1/xi))/sigma.1,0)
  if (for.max) return(cbind(fGP.n,p.2.vec)) # n.sim by (1+length(u2.vec)) matrix
  # bin-GP density at each observation greater than u1
  fGP.r <- sapply(z.gt.u1,function(xx) (1+xi*(xx-u1)/sigma.1)^(-(1+1/xi))/sigma.1)
  #
  cbind(fGP.r,p.2.vec) # n.sim by (length(z.gt.u)+length(u2.vec)) matrix
}# .............. end of function my.stuff()

#################################################################################
