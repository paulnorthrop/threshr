#context("yjthresh vs bcthresh")

# Set a tolerance for the comparison of the simulated values

my_tol <- 1e-5
seed <- 29082019

# We check that, for positive data, yjthresh() and bcthresh() give equivalent
# results, when we call the former using data and the latter with data + 1

# Set a prior: flat for GP parameters, Haldane for P(exceedance)
prior_args <- list(prior = "flatflat", bin_prior = "haldane",
                   h_prior = list(min_xi = -Inf))

## Gulf of Mexico significant wave heights ------------------

gprobs <- c(0.1, 0.5, 0.9)
glambda <- seq(1, 3, 1)

set.seed(seed)
gom_yj_args <- list(data = gom, probs = gprobs, lambda = glambda, n_v = 2,
                    n = 1000)
gom_yj <- do.call(yjthresh, c(gom_yj_args, prior_args))
set.seed(seed)
gom_bc_args <- list(data = gom + 1L, probs = gprobs, lambda = glambda, n_v = 2,
                    n = 1000)
gom_bc <- do.call(bcthresh, c(gom_bc_args, prior_args))

test_that("yjthresh vs bcthresh: pred_perf", {
  testthat::expect_equal(gom_yj$pred_perf, gom_bc$pred_perf,
                         tolerance = my_tol)
})

test_that("yjthresh vs bcthresh: sim_vals", {
  testthat::expect_equal(gom_yj$sim_vals, gom_bc$sim_vals,
                         tolerance = my_tol)
})
