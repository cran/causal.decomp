#' Synthetic Data for Illustration
#'
#' A ranomly generated dataset containing 1000 rows and 9 columns with no
#' missing values.
#'
#' @usage sdata
#'
#' @format A data frame containing the following variables. The data are
#'   provided only for explanatory purposes. The mediators are assumed to be
#'   independent of each other.
#' \describe{
#'   \item{C.num:}{ A quantitative covariate.}
#'   \item{C.bin:}{ A binary covariates with a value of 0 or 1.}
#'   \item{R:}{ A group indicator with four levels.}
#'   \item{S:}{ A quantitative intermediate confounder between a mediator and the outcome.}
#'   \item{M.num:}{ A quantitative mediator.}
#'   \item{M.bin:}{ A binary mediator with a value of 0 or 1.}
#'   \item{M.cat:}{ A categorical mediator with three levels.}
#'   \item{Y.num:}{ A quantitative outcome.}
#'   \item{Y.bin:}{ A binary outcome with a value of 0 or 1.}
#' }
#'
#' @details Note that all the variables are randomly generated using the dataset
#'   used in Park et al. (2021+).
#'
#' @references Park, S., Kang, S., and Lee, C. (2021+). "Choosing an Optimal Method
#'   for Causal Decomposition Analysis: A Better Practice for Identifying
#'   Contributing Factors to Health Disparities".
#'
#' @keywords datasets
"sdata"
