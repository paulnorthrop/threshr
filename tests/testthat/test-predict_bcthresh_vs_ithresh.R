#context("predict: bcthresh vs ithresh")

# We check that predict.bcthresh aress predict.ithresh, for different choices
# of the argument type.

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

## Density

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

# Make $x and $y returned from predict.ithresh arrays, like predict.bcthresh
c1d$x <- array(c1d$x, dim = c(nrow(c1d$x), 1, 1))
c1d$y <- array(c1d$y, dim = c(nrow(c1d$y), 1, 1))
c2d$x <- array(c2d$x, dim = c(nrow(c2d$x), 1, 1))
c2d$y <- array(c2d$y, dim = c(nrow(c2d$y), 1, 1))
c3d$x <- array(c3d$x, dim = c(nrow(c3d$x), 1, 1))
c3d$y <- array(c3d$y, dim = c(nrow(c3d$y), 1, 1))

# Equal apart from the class, so used equivalent
test_that("bcthresh vs ithresh: pred_perf", {
  testthat::expect_equal(c1d, b1d, tolerance = my_tol, ignore_attr = TRUE)
})
test_that("bcthresh vs ithresh: pred_perf", {
  testthat::expect_equal(c2d, b2d, tolerance = my_tol, ignore_attr = TRUE)
})
test_that("bcthresh vs ithresh: pred_perf", {
  testthat::expect_equal(c3d, b3d, tolerance = my_tol, ignore_attr = TRUE)
})

## Distribution function

b1p <- predict(gom_lambda, lambda = 1, type = "p")
b2p <- predict(gom_lambda, lambda = 2, type = "p")
b3p <- predict(gom_lambda, lambda = 3, type = "p")

res1 <- choose_lambda(gom_lambda, lambda = 1)
res2 <- choose_lambda(gom_lambda, lambda = 2)
res3 <- choose_lambda(gom_lambda, lambda = 3)

c1p <- predict(res1, type = "p")
# Fiddle with the xs to make them the same as those used in predict.bcthresh()
c2p <- predict(res2, type = "p", x = threshr:::bc(b2p$x, res2$lambda))
c3p <- predict(res3, type = "p", x = threshr:::bc(b3p$x, res3$lambda))

# Make $x and $y returned from predict.ithresh arrays, like predict.bcthresh
c1p$x <- array(c1p$x, dim = c(nrow(c1p$x), 1, 1))
c1p$y <- array(c1p$y, dim = c(nrow(c1p$y), 1, 1))
c2p$x <- array(c2p$x, dim = c(nrow(c2p$x), 1, 1))
c2p$y <- array(c2p$y, dim = c(nrow(c2p$y), 1, 1))
c3p$x <- array(c3p$x, dim = c(nrow(c3p$x), 1, 1))
c3p$y <- array(c3p$y, dim = c(nrow(c3p$y), 1, 1))

# Equal apart from the class, so used equivalent
test_that("bcthresh vs ithresh: pred_perf", {
  testthat::expect_equal(c1p, b1p, tolerance = my_tol, ignore_attr = TRUE)
})
# Why are these different? ... because the xs are equally-spaced on different scales
test_that("bcthresh vs ithresh: pred_perf", {
  testthat::expect_equal(c2p, b2p, tolerance = my_tol, ignore_attr = TRUE)
})
test_that("bcthresh vs ithresh: pred_perf", {
  testthat::expect_equal(c3p, b3p, tolerance = my_tol, ignore_attr = TRUE)
})

