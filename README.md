
<!-- README.md is generated from README.Rmd. Please edit that file -->
threshr: Threshold Selection and Uncertainty for Extreme Value Analysis
-----------------------------------------------------------------------

### What does threshr do?

The `threshr` package provides functions for the selection of extreme value threshold. At the moment only the simplest case, where the data can be treated as independent identically distributed observations, is considered. Future releases will tackle more general situations.

### A simple example

The main function in the threshr package is `ithresh`, which uses leave-one-out cross-validation in a Bayesian seup to compare the extreme value predictive ability resulting from the use of each of a user-supplied set of thresholds. The following code produces a threshold diagnostic plot using a dataset `gom` containing 315 storm peak significant waveheights. We set a vector `u_vec` of thresholds; call `ithresh`, supplying the data and thresholds; and use then plot the results. In this minimal example (`ithresh` has further arguments) thresholds are judged in terms of the quality of prediction of whether the validation observation lies above the highest threshold in `u_vec` and, if it does, how much it exceeds this highest threshold.

``` r
u_vec <- quantile(gom, probs = seq(0, 0.95, by = 0.05))
gom_cv <- ithresh(data = gom, u_vec = u_vec)
plot(gom_cv)
```

### Installation

To get the current released version from CRAN:

``` r
install.packages("threshr")
```

### Vignette

See `vignette("threshr-vignette", package = "threshr")` for an overview of the package.