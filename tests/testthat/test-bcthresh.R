#context("bcthresh")

# Set a tolerance for the comparison of the simulated values

my_tol <- 1e-5
seed <- 18072019

# 1. We check that the results from calling bcthresh() with lambda = 1 agree
#    with those from an equivalent call to ithresh().

probs <- c(0.05, 0.5, 0.9)
u_vec_gom <- quantile(gom, probs = probs)

set.seed(seed)
res1 <- ithresh(data = gom, u_vec = u_vec_gom, n_v = 2)
set.seed(seed)
res2 <- bcthresh(data = gom, probs = probs, lambda = 1, n_v = 2)

test_that("bcthresh vs ithresh: pred_perf", {
  testthat::expect_equal(res1$pred_perf, res2$pred_perf[, , 1],
                         tolerance = my_tol)
})

# 2. Repeat for trans = "BC"

set.seed(seed)
res1 <- ithresh(data = gom, u_vec = u_vec_gom, n_v = 2, trans = "BC")
set.seed(seed)
res2 <- bcthresh(data = gom, probs = probs, lambda = 1, n_v = 2, trans = "BC")

test_that("bcthresh vs ithresh: pred_perf, trans = ''BC''", {
  testthat::expect_equal(res1$pred_perf, res2$pred_perf[, , 1],
                         tolerance = my_tol)
})

