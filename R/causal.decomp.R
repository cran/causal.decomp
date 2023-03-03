#' causal.decomp: Causal Decomposition Analysis.
#' 
#' The causal.decomp package provides four important functions:
#' mmi, smi, pocr, and sensitivity.
#' 
#' @docType package
#' @name causal.decomp
#' 
#' @importFrom MASS gamma.shape polr
#' @importFrom nnet multinom
#' @importFrom PSweight SumStat
#' @importFrom parallel detectCores mclapply
#' @importFrom stats Gamma as.formula binomial coef
#'   formula gaussian glm inverse.gaussian lm
#'   model.frame pcauchy plogis pnorm poisson
#'   predict quantile rbinom rgamma rmultinom
#'   rnorm rpois terms update
#' @importFrom utils packageDescription
NULL
#> NULL