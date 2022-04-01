#context("choose_lambda")

seed <- 29082019
my_tol <- 1e-5

# Set a prior: flat for GP parameters, Haldane for P(exceedance)
prior_args <- list(prior = "flatflat", bin_prior = "haldane",
                   h_prior = list(min_xi = -Inf))

gprobs <- seq(0.1, 0.9, 0.4)
glambda <- c(1, 1.5)
gom_args <- list(data = gom, probs = gprobs, lambda = glambda)

set.seed(seed)
gom_lambda <- do.call(bcthresh, c(gom_args, prior_args))
from_bcthresh <- choose_lambda(gom_lambda, lambda = 1)

set.seed(seed)
ithresh_gom_args <- c(gom_args, list(u_vec = quantile(gom, probs = gprobs)))
ithresh_gom_args$probs <- NULL
from_ithresh <- do.call(ithresh, c(ithresh_gom_args, prior_args))

test_that("bcthresh (after choose_lambda) vs ithresh: pred_perf", {
  testthat::expect_equal(from_bcthresh$pred_perf, from_ithresh$pred_perf,
                         tolerance = my_tol)
})

test_that("bcthresh (after choose_lambda) vs ithresh: sim_vals", {
  testthat::expect_equal(from_bcthresh$sim_vals, from_ithresh$sim_vals,
                         tolerance = my_tol)
})

test_that("bcthresh (after choose_lambda) vs ithresh: u_ps", {
  testthat::expect_equal(from_bcthresh$u_ps, from_ithresh$u_ps,
                         tolerance = my_tol)
})

test_that("bcthresh (after choose_lambda) vs ithresh: v_ps", {
  testthat::expect_equal(from_bcthresh$v_ps, from_ithresh$v_ps,
                         tolerance = my_tol)
})

# Note: BC(1) subtracts 1 from the data

test_that("bcthresh (after choose_lambda) vs ithresh: data", {
  testthat::expect_equal(from_bcthresh$data + 1, from_ithresh$data,
                         tolerance = my_tol)
})

test_that("bcthresh (after choose_lambda) vs ithresh: u_vec", {
  testthat::expect_equal(from_bcthresh$u_vec + 1, from_ithresh$u_vec,
                         tolerance = my_tol)
})

test_that("bcthresh (after choose_lambda) vs ithresh: v_vec", {
  testthat::expect_equal(from_bcthresh$v_vec + 1, from_ithresh$v_vec,
                         tolerance = my_tol)
})

