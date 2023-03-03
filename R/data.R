#' Synthetic Data for Illustration
#'
#' A randomly generated dataset containing 1000 rows and 9 columns with no
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
#'   \item{X:}{ A quantitative intermediate confounder between a mediator and the outcome.}
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


#' Synthetic Data Generated Based on the Midlife Development in the U.S. (MIDUS) Study
#'
#' This is a synthetic dataset that includes variables from the Midlife Development
#' in the U.S. (MIDUS) study. It has been artificially generated based on the actual
#' MIDUS data, which is not publicly available due to confidentiality concerns.
#' The synthetic data set consists of 1948 rows and 9 columns, with no missing values.
#'
#' @usage sMIDUS
#'
#' @format A data frame containing the following variables.
#' \describe{
#'   \item{health:}{ cardiovascular health score.}
#'   \item{racesex:}{ race-gender groups with four levels (1: White men, 2: White women, 3: Black men, 4: Black women).}
#'   \item{lowchildSES:}{ socioeconomic status (SES) in the childhood.}
#'   \item{abuse:}{ adverse experience in the childhood.}
#'   \item{edu:}{ education level.}
#'   \item{age:}{ age.}
#'   \item{stroke:}{ genetic vulnerability with a value of 0 or 1.}
#'   \item{T2DM:}{ genetic vulnerability with a value of 0 or 1.}
#'   \item{heart:}{ genetic vulnerability with a value of 0 or 1.}
#' }
#'
#' @details Note that all the variables are fabricated using the actual MIDUS data
#'   used in Park et al. (2023).
#'
#' @references Park, S., Kang, S., Lee, C., & Ma, S. (2023). Sensitivity Analysis
#'   for Causal Decomposition Analysis: Assessing Robustness Toward Omitted Variable Bias,
#'   Journal of Causal Inference. Forthcoming.
#'
#' @keywords datasets
"sMIDUS"
