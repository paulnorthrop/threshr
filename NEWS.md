# threshr 1.0.6.9000

## Bug fixes and minor improvements

* Corrected a typo in the description of of the component `sim_vals` of the object returned from `ithresh`. `u_vec[i]` should be `u_vec[j]`. Thank you to Ye Liu for spotting this. [Issue 3](https://github.com/paulnorthrop/threshr/issues/3#issue-2832668692)

# threshr 1.0.6

## Bug fixes and minor improvements

* Fixed 3 \link{} targets with missing (revdbayes) package anchors in the Rd file for `ithresh()`.

# threshr 1.0.5

## Bug fixes

* Fixed issues with the incorrect use of \itemize in some Rd files.

# threshr 1.0.4

## Bug fixes and minor improvements

* Create the help file for the package correctly, with alias threshr-package.

* Activated 3rd edition of the `testthat` package

* README.md: Used app.codecov.io as base for codecov link.

# threshr 1.0.3

## Bug fixes and minor improvements

* Tests in `test-box_cox.R` and `test-inv_box_cox.R` have been modified to avoid errors in the upcoming new release of the `testthat` package.

# threshr 1.0.2

## Bug fixes and minor improvements

* In `ithresh()` a user-supplied (log-)prior R function can now be set for the binomial probability p of threshold exceedance.  The functionality requires at least version 1.3.4 of the revdbayes package.

* A print method for class `ithresh` has been added.

* In `plot.ithresh()` a more informative error message is given if an inappropriate value of the argument `which_v` is supplied.

* In `predict.ithresh()` further arguments can now be passed to `revdbayes::predict.evpost`.  In particular, the level(s) of predictive intervals can be set.  An example has been added to the documentation.

* pkgdown documentation at [https://paulnorthrop.github.io/threshr/](https://paulnorthrop.github.io/threshr/)

* revdbayes:: is used instead of revdbayes::: to avoid CRAN package check NOTEs.

# threshr 1.0.1

## Bug fixes and minor improvements

* Some examples and tests are modified slightly to avoid using unrealistically high or low thresholds.

* Dependence on R version changed to 3.3.0 to avoid CRAN NOTE.



