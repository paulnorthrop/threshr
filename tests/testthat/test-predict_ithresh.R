context("predict.ithresh function")

# We check predict.ithresh, for different choices of the argument type.

# Set a tolerance for the comparison of the simulated values
my_tol <- 1e-5

# Check that calling with type = "q" (quantiles) with the default
# probabilities p = c(0.025, 0.25, 0.5, 0.75, 0.975) and then calling
# with the results and with type = "p" gets us back to the initial
# probabilities.

u_vec_gom <- quantile(gom, probs = seq(0, 0.9, by = 0.1))
gom_cv <- ithresh(data = gom, u_vec = u_vec_gom)

# (Default: which_u = "best")
qs <- predict(gom_cv, type = "q", n_years = c(100, 1000))$y
ps <- predict(gom_cv, type = "p", x = qs, n_years = c(100, 1000))$y
check_ps <- matrix(c(0.025, 0.25, 0.5, 0.75, 0.975), 5, 2)

testthat::expect_equal(ps, check_ps, tolerance = my_tol)
