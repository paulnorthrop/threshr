context("predict: bcthresh vs ithresh")

# We check that predict.bcthresh aress predict.ithresh, for different choices of the argument type.

set.seed(1092019)

# Set a tolerance for the comparison of the simulated values
my_tol <- 1e-5

# Set a prior: flat for GP parameters, Haldane for P(exceedance)
prior_args <- list(prior = "flatflat", bin_prior = "haldane",
                   h_prior = list(min_xi = -Inf))

## Gulf of Mexico significant wave heights ------------------

gprobs <- seq(0.1, 0.9, 0.1)
glambda <- seq(1, 3, 0.5)
gom_args <- list(data = gom, probs = gprobs, lambda = glambda)
gom_lambda <- do.call(bcthresh, c(gom_args, prior_args))

b1d <- predict(gom_lambda, lambda = 1, type = "d")
b2d <- predict(gom_lambda, lambda = 2, type = "d")
b3d <- predict(gom_lambda, lambda = 3, type = "d")

res1 <- choose_lambda(gom_lambda, lambda = 1)
res2 <- choose_lambda(gom_lambda, lambda = 2)
res3 <- choose_lambda(gom_lambda, lambda = 3)

c1d <- predict(res1, type = "d")
# Fiddle with the xs to make them the same as those used in predict.bcthresh()
c2d <- predict(res2, type = "d", x = threshr:::bc(b2d$x, res2$lambda))
c3d <- predict(res3, type = "d", x = threshr:::bc(b3d$x, res3$lambda))

# Equal apart from the class, so used equivalent
test_that("bcthresh vs ithresh: pred_perf", {
  testthat::expect_equivalent(c1d, b1d, tolerance = my_tol)
})
# Why are these different? ... because the xs are equally-spaced on different scales
test_that("bcthresh vs ithresh: pred_perf", {
  testthat::expect_equivalent(c2d, b2d, tolerance = my_tol)
})
test_that("bcthresh vs ithresh: pred_perf", {
  testthat::expect_equivalent(c3d, b3d, tolerance = my_tol)
})
