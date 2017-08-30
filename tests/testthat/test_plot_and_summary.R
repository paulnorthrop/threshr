context("plot and summary")

# Set a tolerance for the comparison of the simulated values

my_tol <- 1e-5
seed <- 27082017

# 1. Check that the output from summary.ithresh is a double matrix and has the
#    correct dimensions.

set.seed(seed)
n_v <- 4
u_vec_gom <- quantile(gom, probs = seq(0, 0.9, len = n_v))
gom_cv <- ithresh(data = gom, u_vec = u_vec_gom, n_v = 4)
check_summary <- summary(gom_cv)

testthat::expect_type(check_summary, "double")
testthat::expect_is(check_summary, "matrix")
testthat::expect_equal(nrow(check_summary), n_v)
testthat::expect_equal(ncol(check_summary), 5)

# 2. Check that the columns of the object returned by plot.ithresh (when
#    which_u is not supplied) sum to 1 if na.rm = TRUE and c(NA, ..., NA, 1)
#    if na.rm = FALSE.

testthat::expect_equal(colSums(pjn$y, na.rm = TRUE), rep(1, n_v),
                       tolerance = my_tol)
testthat::expect_equal(colSums(pjn$y), c(rep(NA, n_v - 1), 1),
                       tolerance = my_tol)
