context("Box-Cox")

# Check that bc() is correct for lambda = 0

x <- 1:10
test_that("bc for lambda = 0, x > 0", {
  testthat::expect_equal(bc(x = x, lambda = 0), log(x))
})

# Check that bc() throws an error when x < 0

x <- -1
temp <- try(bc(x = x, lambda = 0), silent = TRUE)
test_that(paste("bc for lambda = 0, x = ", x = x), {
  testthat::expect_equal(attr(temp, "class"), "try-error")
})

# Check that bc() is correct for lambda very slightly smaller in
# magnitude than lambda_tol = 1 / 50 and m (Taylor series polynomial order)
# is large

eps <- 1e-10
m <- 10
lambda_tol <- 1 / 50
x <- 1:10

# lambda very slightly less than lambda_tol

lambda <- lambda_tol - eps
check_val <- (x ^ lambda - 1) / lambda
for (i in 1:length(x)) {
  test_that(paste("bc, 0 < lambda < lambda_tol, x = ", x = x[i]), {
    testthat::expect_equal(bc(x = x[i], lambda = lambda,
                                 lambda_tol = lambda_tol, m = m),
                           check_val[i])
  })
}

# lambda very slightly greater than -lambda_tol

lambda <- -lambda_tol + eps
check_val <- (x ^ lambda - 1) / lambda
for (i in 1:length(x)) {
  test_that(paste("bc, -lambda_tol < lambda < 0, x = ", x = x[i]), {
    testthat::expect_equal(bc(x = x[i], lambda = lambda,
                                 lambda_tol = lambda_tol, m = m),
                           check_val[i])
  })
}
