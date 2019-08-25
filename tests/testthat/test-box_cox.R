context("Box-Cox")

# Check that bc_gm() is correct for lambda = 0

x <- 1:10
for (i in length(x)) {
  test_that(paste("bc_gm for lambda = 0, x = ", x = x), {
    testthat::expect_equal(bc_gm(x = x, lambda = 0), log(x))
  })
}

# Check that bc_gm() throws an error when x < 0

x <- -1
for (i in length(x)) {
  temp <- try(bc_gm(x = x, lambda = 0))
  test_that(paste("bc_gm for lambda = 0, x = ", x = x), {
    testthat::expect_equal(attr(temp, "class"), "try-error")
  })
}


# Check that box_cox_deriv is correct for lambda very slightly smaller in
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
  test_that(paste("box_cox, 0 < lambda < lambda_tol, x = ", x = x[i]), {
    testthat::expect_equal(bc_gm(x = x[i], lambda = lambda,
                                 lambda_tol = lambda_tol, m = m),
                           check_val[i])
  })
}

# lambda very slightly greater than -lambda_tol

lambda <- -lambda_tol + eps
check_val <- (x ^ lambda - 1) / lambda
for (i in 1:length(x)) {
  test_that(paste("box_cox_deriv, -lambda_tol < lambda < 0, x = ", x = x[i]), {
    testthat::expect_equal(bc_gm(x = x[i], lambda = lambda,
                                 lambda_tol = lambda_tol, m = m),
                           check_val[i])
  })
}
