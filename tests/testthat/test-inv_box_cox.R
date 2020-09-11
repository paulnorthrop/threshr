context("Inverse Box-Cox")

# Check that inv_bc() is correct for lambda = 0

x <- 1:10
test_that("inv_bc for lambda = 0, x > 0", {
  testthat::expect_equal(inv_bc(x = x, lambda = 0), exp(x))
})

# Check that inv_bc() throws an error when x is out of range

lambda <- 2
x <- -1 / lambda
temp <- try(inv_bc(x = x, lambda = lambda), silent = TRUE)
test_that(paste("bc for lambda = 2, x is out of range, x = ", x = x), {
  testthat::expect_equal(attr(temp, "class"), "try-error")
})

lambda <- 2
x <- -1 / lambda - 0.0001
temp <- try(inv_bc(x = x, lambda = lambda), silent = TRUE)
test_that(paste("bc for lambda = 2, x is out of range, x = ", x = x), {
  testthat::expect_equal(attr(temp, "class"), "try-error")
})

lambda <- -2
x <- -1 / lambda
temp <- try(inv_bc(x = x, lambda = lambda), silent = TRUE)
test_that(paste("bc for lambda = -2, x is out of range, x = ", x = x), {
  testthat::expect_equal(attr(temp, "class"), "try-error")
})

lambda <- -2
x <- -1 / lambda + 0.0001
temp <- try(inv_bc(x = x, lambda = lambda), silent = TRUE)
test_that(paste("bc for lambda = -2, x is out of range, x = ", x = x), {
  testthat::expect_equal(attr(temp, "class"), "try-error")
})

# Check that inv_box_cox() is correct for lambda very slightly smaller in
# magnitude than lambda_tol = 1 / 50 and m (Taylor series polynomial order)
# is large

eps <- 1e-10
m <- 100
lambda_tol <- 1 / 50
x <- 1:10

# lambda very slightly less than lambda_tol

lambda <- lambda_tol - eps
check_val <- (1 + lambda * x) ^ (1 / lambda)
for (i in 1:length(x)) {
  test_that(paste("inv_box_cox, 0 < lambda < lambda_tol, x = ", x = x[i]), {
    testthat::expect_equal(inv_bc(x = x[i], lambda = lambda,
                                  lambda_tol = lambda_tol, m = m),
                           check_val[i])
  })
}

# lambda very slightly greater than -lambda_tol

lambda <- -lambda_tol + eps
check_val <- (1 + lambda * x) ^ (1 / lambda)
for (i in 1:length(x)) {
  test_that(paste("inv_box_cox, -lambda_tol < lambda < 0, x = ", x = x[i]), {
    testthat::expect_equal(inv_bc(x = x[i], lambda = lambda,
                                  lambda_tol = lambda_tol, m = m),
                           check_val[i])
  })
}

# Check that inv_box_cox() is correct for lambda larger in magnitude than
# lambda_tol = 1 / 50

# lambda very slightly larger than lambda_tol

lambda <- lambda_tol + eps
check_val <- (1 + lambda * x) ^ (1 / lambda)
for (i in 1:length(x)) {
  test_that(paste("inv_box_cox, 0 < lambda < lambda_tol, x = ", x = x[i]), {
    testthat::expect_equal(inv_bc(x = x[i], lambda = lambda,
                                  lambda_tol = lambda_tol),
                           check_val[i])
  })
}

# lambda very slightly smaller than -lambda_tol

lambda <- -lambda_tol - eps
check_val <- (1 + lambda * x) ^ (1 / lambda)
for (i in 1:length(x)) {
  test_that(paste("inv_box_cox, -lambda_tol < lambda < 0, x = ", x = x[i]), {
    testthat::expect_equal(inv_bc(x = x[i], lambda = lambda,
                                  lambda_tol = lambda_tol, m = m),
                           check_val[i])
  })
}
