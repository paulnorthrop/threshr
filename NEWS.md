# threshr 1.0.2

## Bug fixes and minor improvements

* In `ithresh()` a user-supplied (log-)prior R function can now be set for the binomial probability p of threshold exceedance.  The functionality requires at least version 1.3.4 of the revdbayes package.

* A print method for class `ithresh` has been added.

* In `plot.ithresh()` a more informative error message is given if an inappropriate value of the argument `which_v` is supplied.

* In `predict.ithresh()` further arguments can now be passed to `revdbayes::predict.evpost`.  In particular, the level(s) of predictive intervals can be set.  An example has been added to the documentation.

# threshr 1.0.1

## Bug fixes and minor improvements

* Some examples and tests are modified slightly to avoid using unrealistically high or low thresholds.

* Dependence on R version changed to 3.3.0 to avoid CRAN NOTE.



