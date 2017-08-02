#' threshr: Threshold Selection and Uncertainty for Extreme Value Analysis
#'
#' Provides functions for the selection of extreme value threshold.
#' See the 'threshr' website for more information, documentation and examples.
#'
#' @details The main function in the threshr package is
#'
#'   See \code{vignette("threshr-vignette", package = "threshr")} for an
#'   overview of the package.
#'
#' @references Northrop, P. J. (2016). revdbayes: Ratio-of-Uniforms Sampling
#'   for Bayesian Extreme Value Analysis. R package version 1.1.0.
#'   \url{https://cran.r-project.org/package=revdbayes}.
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{http://dx.doi.org/10.1111/rssc.12159}
#'
#' @seealso \code{\link[revdbayes]{set_prior}} and
#'   \code{\link[revdbayes]{rpost}}: functions from the
#'   \code{\link{revdbayes}} package used to set an extreme value prior density
#'   and to sample from a extreme value posterior density respectively.
#' @docType package
#' @name threshr
#' @import methods
NULL

#' Storm peak significant wave heights from the Gulf of Mexico
#'
#' A numeric vector containing 315 hindcasts of storm peak significant wave
#' heights, metres, from 1900 to 2005 at an unnamed location in the Gulf
#' of Mexico.
#'
#'@format A vector containing 315 observations.
#'@source Oceanweather Inc. (2005) GOMOS -- Gulf of Mexico hindcast study.
#'@references Northrop, P. J., N. Attalides, and P. Jonathan. (2016).
#'  Cross-Validatory Extreme Value Threshold Selection and Uncertainty with
#'  Application to Ocean Storm Severity. \emph{Journal of the Royal
#'  Statistical Society: Series C (Applied Statistics)}.
#'  doi:\href{https://doi.org/10.1111/rssc.12159}{10.1111/rssc.12159}.
"gom"

#' Storm peak significant wave heights from the North Sea
#'
#' A numeric vector containing 628 hindcasts of storm peak significant wave
#' heights, metres, from 1964 to 1995 at an unnamed location in the North Sea.
#'
#'@format A vector containing 628 observations.
#'@source Oceanweather Inc. (1995) NEXT -- North Sea hindcast study.
#'@references Northrop, P. J., N. Attalides, and P. Jonathan. (2016).
#'  Cross-Validatory Extreme Value Threshold Selection and Uncertainty with
#'  Application to Ocean Storm Severity. \emph{Journal of the Royal
#'  Statistical Society: Series C (Applied Statistics)}.
#'  doi:\href{https://doi.org/10.1111/rssc.12159}{10.1111/rssc.12159}.
"ns"
