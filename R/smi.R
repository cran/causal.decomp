#' Single-Mediator-Imputation Estimation Method
#'
#' 'smi' is used to estimate the initial disparity, disparity reduction, and
#' disparity remaining for causal decomposition analysis, using the
#' single-mediator-imputation estimation method proposed by Park et. al. (2021+).
#'
#' @details This function returns the point estimates of the initial disparity,
#'   disparity reduction, and disparity remaining for a categorical
#'   treatment and a variety of types of outcome and mediator(s) in causal
#'   decomposition analysis. It also returns nonparametric
#'   percentile bootstrap confidence intervals for each estimate.
#'
#'   The definition of the initial disparity, disparity reduction, and
#'   disparity remaining can be found in help('mmi'). As opposed to the 'mmi'
#'   function, this function uses the single-mediator-imputation method
#'   suggested by Park et al. (2021+). See the reference for more details.
#'
#'   If one wants to make the inference conditional on baseline covariates,
#'   set 'conditional = T' and center the data before fitting the models.
#'
#'   As of version 0.0.1, the mediator model ('fit.m') can be of class 'lm', 'glm',
#'   'multinom', or 'polr', corresponding respectively to the linear regression
#'   models and generalized linear models, multinomial log-linear models, and
#'   ordered response models. The outcome model ('fit.y') can be of class 'lm' or 'glm'.
#'   Also, the treatment model ('fit.r') can be of class 'CBPS' or 'SumStat', both of
#'   which use the propensity score weighting. It is only necessary when 'conditional = F'.
#'
#' @param fit.r a fitted model object for treatment. Can be of class 'CBPS' or
#'   'SumStat'. Default is 'NULL'. Only necessary if 'conditional' is 'FALSE'.
#' @param fit.m a fitted model object for mediator. Can be of class 'lm', 'glm',
#'   'multinom', or 'polr'.
#' @param fit.y a fitted model object for outcome. Can be of class 'lm' or 'glm'.
#' @param sims number of Monte Carlo draws for nonparametric bootstrap.
#' @param conf.level level of the returned two-sided confidence intervals,
#'   which are estimated by the nonparametric percentile bootstrap method.
#'   Default is .95, which returns the 2.5 and 97.5 percentiles of the simulated
#'   quantities.
#' @param treat a character string indicating the name of the treatment variable
#'   used in the models. The treatment can be categorical with two or more
#'   categories (two- or multi-valued factor).
#' @param covariates a vector containing the name of the covariate variable(s)
#'   used in the models. Each covariate can be categorical with two or more
#'   categories (two- or multi-valued factor) or continuous (numeric).
#' @param conditional a logical value. If 'TRUE', the function will return the
#'   estimates conditional on those covariate values; and all covariates in
#'   mediator and outcome models need to be centered prior to fitting.
#'   Default is 'TRUE'. If 'FALSE', 'fit.r' needs to be specified.
#' @param cluster a vector of cluster indicators for the bootstrap. If provided,
#'   the cluster bootstrap is used. Default is 'NULL'.
#' @param long a logical value. If 'TRUE', the output will contain the entire
#'   sets of estimates for all bootsrap samples. Default is 'FALSE'.
#' @param mc.cores The number of cores to use. Must be exactly 1 on Windows.
#'
#' @return
#'
#'   \item{result}{a matrix containing the point estimates of the initial disparity,
#'   disparity remaining, and disparity reduction, and the percentile bootsrap
#'   confidence intervals for each estimate.}
#'   \item{all.result}{a matrix containing the point estimates of the initial disparity,
#'   disparity remaining, and disparity reduction for all bootsrap samples. Returned
#'   if 'long' is 'TRUE'.}
#'
#' @author
#'   Suyeon Kang, University of California, Riverside, \email{skang062@@ucr.edu};
#'   Soojin Park, University of California, Riverside, \email{soojinp@@ucr.edu}.
#'
#' @seealso \code{\link{mmi}}
#'
#' @references
#'   Park, S., Kang, S., and Lee, C. (2021+). "Choosing an Optimal Method for Causal
#'   Decomposition Analysis: A Better Practice for Identifying Contributing Factors to
#'   Health Disparities".
#'
#' @export
#' @examples
#' data(sdata)
#'
#' #-----------------------------------------------------------#
#' # Example 1-a: Continuous Outcome
#' #-----------------------------------------------------------#
#' require(PSweight)
#' fit.r1 <- SumStat(R ~ C.num + C.bin, data = sdata, weight = "IPW")
#' \donttest{require(CBPS)
#' fit.r2 <- CBPS(R ~ C.num + C.bin, data = sdata, method = "exact",
#'           standardize = "TRUE")}
#'
#' # Continuous mediator
#' fit.m1 <- lm(M.num ~ R + C.num + C.bin, data = sdata)
#' fit.y1 <- lm(Y.num ~ R + M.num + S + C.num + C.bin, data = sdata)
#' res.1a1 <- smi(fit.r = fit.r1, fit.m = fit.m1,
#'           fit.y = fit.y1, sims = 40, conditional = FALSE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.1a1
#'
#' \donttest{# Binary mediator
#' fit.m2 <- glm(M.bin ~ R + C.num + C.bin, data = sdata,
#'           family = binomial(link = "logit"))
#' fit.y2 <- lm(Y.num ~ R + M.bin + S + C.num + C.bin, data = sdata)
#' res.1a2 <- smi(fit.r = fit.r1, fit.m = fit.m2,
#'           fit.y = fit.y2, sims = 40, conditional = FALSE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.1a2
#'
#' # Categorical mediator
#' require(MASS)
#' fit.m3 <- polr(M.cat ~ R + C.num + C.bin, data = sdata)
#' fit.y3 <- lm(Y.num ~ R + M.cat + S + C.num + C.bin, data = sdata)
#' res.1a3 <- smi(fit.r = fit.r1, fit.m = fit.m3,
#'           fit.y = fit.y3, sims = 40, conditional = FALSE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.1a3
#'
#' require(nnet)
#' fit.m4 <- multinom(M.cat ~ R + C.num + C.bin, data = sdata)
#' res.1a4 <- smi(fit.r = fit.r1, fit.m = fit.m4,
#'           fit.y = fit.y3, sims = 40, conditional = FALSE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.1a4}
#' 
#' #-----------------------------------------------------------#
#' # Example 1-b: Binary Outcome
#' #-----------------------------------------------------------#
#' \donttest{# Continuous mediator
#' fit.y1 <- glm(Y.bin ~ R + M.num + S + C.num + C.bin,
#'           data = sdata, family = binomial(link = "logit"))
#' res.1b1 <- smi(fit.r = fit.r1, fit.m = fit.m1,
#'           fit.y = fit.y1, sims = 40, conditional = FALSE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.1b1
#'
#' # Binary mediator
#' fit.y2 <- glm(Y.bin ~ R + M.bin + S + C.num + C.bin,
#'           data = sdata, family = binomial(link = "logit"))
#' res.1b2 <- smi(fit.r = fit.r1, fit.m = fit.m2,
#'           fit.y = fit.y2, sims = 40, conditional = FALSE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.1b2
#'
#' # Categorical mediator
#' fit.y3 <- glm(Y.bin ~ R + M.cat + S + C.num + C.bin,
#'           data = sdata, family = binomial(link = "logit"))
#' res.1b3 <- smi(fit.r = fit.r1, fit.m = fit.m3,
#'           fit.y = fit.y3, sims = 40, conditional = FALSE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.1b3
#'
#' res.1b4 <- smi(fit.r = fit.r1, fit.m = fit.m4,
#'           fit.y = fit.y3, sims = 40, conditional = FALSE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.1b4}
#'
#' #--------------------------------------------------------------#
#' # Example 2-a: Continuous Outcome, Conditional on Covariates
#' #--------------------------------------------------------------#
#' \donttest{# For conditional = T, need to create data with centered covariates
#' # copy data
#' sdata.c <- sdata
#' # center quantitative covariate(s)
#' sdata.c$C.num <- scale(sdata.c$C.num, center = TRUE, scale = FALSE)
#' # center binary (or categorical) covariates(s)
#' # only neccessary if the desired baseline level is NOT the default baseline level.
#' sdata.c$C.bin <- relevel(sdata.c$C.bin, ref = "1")
#'
#' # Continuous mediator
#' fit.y1 <- lm(Y.num ~ R + M.num + S + C.num + C.bin, data = sdata.c)
#' fit.m1 <- lm(M.num ~ R + C.num + C.bin, data = sdata.c)
#' res.2a1 <- smi(fit.m = fit.m1,
#'           fit.y = fit.y1, sims = 40, conditional = TRUE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.2a1
#'
#' # Binary mediator
#' fit.y2 <- lm(Y.num ~ R + M.bin + S + C.num + C.bin, data = sdata.c)
#' fit.m2 <- glm(M.bin ~ R + C.num + C.bin, data = sdata.c,
#'           family = binomial(link = "logit"))
#' res.2a2 <- smi(fit.m = fit.m2,
#'           fit.y = fit.y2, sims = 40, conditional = TRUE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.2a2
#'
#' # Categorical mediator
#' fit.y3 <- lm(Y.num ~ R + M.cat + S + C.num + C.bin, data = sdata.c)
#' fit.m3 <- polr(M.cat ~ R + C.num + C.bin, data = sdata.c)
#' res.2a3 <- smi(fit.m = fit.m3,
#'           fit.y = fit.y3, sims = 40, conditional = TRUE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.2a3
#'
#' fit.m4 <- multinom(M.cat ~ R + C.num + C.bin, data = sdata.c)
#' res.2a4 <- smi(fit.m = fit.m4,
#'           fit.y = fit.y3, sims = 40, conditional = TRUE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.2a4}
#'
#' #-----------------------------------------------------------#
#' # Example 2-b: Binary Outcome, Conditional on Covariates
#' #-----------------------------------------------------------#
#' \donttest{# Continuous mediator
#' fit.y1 <- glm(Y.bin ~ R + M.num + S + C.num + C.bin,
#'           data = sdata.c, family = binomial(link = "logit"))
#' res.2b1 <- smi(fit.m = fit.m1,
#'           fit.y = fit.y1, sims = 40, conditional = TRUE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.2b1
#'
#' # Binary mediator
#' fit.y2 <- glm(Y.bin ~ R + M.bin + S + C.num + C.bin,
#'           data = sdata.c, family = binomial(link = "logit"))
#' res.2b2 <- smi(fit.m = fit.m2,
#'           fit.y = fit.y2, sims = 40, conditional = TRUE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.2b2
#'
#' # Categorical mediator
#' fit.y3 <- glm(Y.bin ~ R + M.cat + S + C.num + C.bin,
#'           data = sdata.c, family = binomial(link = "logit"))
#' res.2b3 <- smi(fit.m = fit.m3,
#'           fit.y = fit.y3, sims = 40, conditional = TRUE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.2b3
#'
#' res.2b4 <- smi(fit.m = fit.m4,
#'           fit.y = fit.y3, sims = 40, conditional = TRUE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.2b4}
smi <- function(fit.r = NULL, fit.m, fit.y,
                treat, covariates, sims = 100, conf.level = .95,
                conditional = TRUE, cluster = NULL, long = FALSE,
                mc.cores = 1L){

  # Warning for inappropriate settings
  if(!is.null(fit.r) && conditional){
    stop("treatment model must be NULL when conditional = TRUE")
  }

  if(is.null(fit.r) && !conditional){
    stop("treatment model must be specified when conditional = FALSE")
  }

  if(length(treat) >= 2){
    stop("only one variable must be selected for 'treat'")
  }

  if(conf.level <= 0 | conf.level >= 1){
    stop("'conf.level' must be between 0 and 1")
  }

  if(.Platform$OS.typ == "windows" && as.integer(mc.cores) > 1L){
    warning("'mc.cores' > 1 is not supported on Windows. mc.cores = 1 forced")
    mc.cores <- 1L
  } else if (as.integer(mc.cores) > detectCores()){
    warning(paste("The maximum number of cores is ", detectCores(),
                  ". mc.cores = ", detectCores(), " forced", sep = ""))
    mc.cores <- detectCores()
  } else if (as.integer(mc.cores) < 1L){
    stop("'mc.cores' must be >= 1")
  }

  if(length(class(fit.m)) == 1 && class(fit.m) == "list"){
    stop("There must be only one mediator model.
         If you have multiple mediators, use 'mmi' function in the package")
  }

  # No possibility of multiple mediators
  isMultiConfounders <- FALSE
  num.ms <- 1

  # Model type indicators
  isCBPS.r <- inherits(fit.r, "CBPS")
  isSumStat.r <- inherits(fit.r, "SumStat")

  isGlm.y <- inherits(fit.y, "glm")
  isLm.y <- inherits(fit.y, "lm")

  isGlm.m <- inherits(fit.m, "glm")
  isLm.m <- inherits(fit.m, "lm")
  isNominal.m <- inherits(fit.m, "multinom")
  isOrdinal.m <- inherits(fit.m, "polr")
  if (isGlm.m) {
    FamilyM <- fit.m$family$family
  }

  if(!is.null(fit.r) && !isCBPS.r && !isSumStat.r){
    stop("unsupported treatment model")
  }

  if(!isGlm.y && !isLm.y){
    stop("unsupported outcome model")
  }

  if(!isGlm.m && !isLm.m && !isNominal.m && !isOrdinal.m){
    stop("unsupported mediator model(s)")
  }

  # Numbers of observations and model frame
  if(!is.null(fit.r) && isCBPS.r){
    n.r <- nrow(model.frame(fit.r))
  } else if (!is.null(fit.r) && isSumStat.r){
    n.r <- nrow(fit.r$propensity)
  } else {
    n.r <- NULL
  }
  n.m <- nrow(model.frame(fit.m))
  y.data <- model.frame(fit.y)
  n.y <- nrow(y.data)
  if(is.null(n.r)){
    if(n.m != n.y){
      stop("number of observations do not match between mediator
           and outcome models")
    }
  } else if (!is.null(n.r)){
    if(n.m != n.y | n.r != n.y | n.m != n.r){
      stop("number of observations do not match between treatment,
           mediator and outcome models")
    }
  }

  # Extract treatment and outcome variable names from models
  if(!is.null(fit.r)){
    if(isCBPS.r){
      treat0 <- all.vars(formula(fit.r))[[1]]
    } else if(isSumStat.r){
      treat0 <- names(fit.r$ps.weight)[[1]]
    }
    if(treat != treat0){
      stop("response variable of 'fit.r' and specified 'treat' do not match")
    } else {
      treat <- treat0
    }
  }
  r <- length(levels(y.data[, treat]))
  outcome <- all.vars(formula(fit.y))[[1]]

  # Bootstrap QoI
  message("Running nonparametric bootstrap\n")

  # Make objects available to all functions
  environment(combine.b) <- environment(combine) <- environment(select.treat) <- environment()

  out <- matrix(NA, 3 * (r - 1), 3)
  all.out <- matrix(NA, 3 * (r - 1), sims)
  colnames(out) <- c("estimate", paste(conf.level * 100, "% CI Lower", sep = ""),
                     paste(conf.level * 100, "% CI Upper", sep = ""))
  rn <- NULL
  for(rr in 2:r){
    rn <- c(rn, c(paste("Initial Disparity   (", levels(y.data[, treat])[1],
                        " vs ", levels(y.data[, treat])[rr], ")", sep = ""),
                  paste("Disparity Remaining (", levels(y.data[, treat])[1],
                        " vs ", levels(y.data[, treat])[rr], ")", sep = ""),
                  paste("Disparity Reduction (", levels(y.data[, treat])[1],
                        " vs ", levels(y.data[, treat])[rr], ")", sep = "")))

  }
  rownames(out) <- rownames(all.out) <- rn

  summary0 <- mclapply(1:sims, combine.b, fit.r = fit.r, fit.s = NULL, fit.m = fit.m, fit.y = fit.y,
                       cluster = cluster, func = "smi", mc.cores = mc.cores)
  out.full <- simplify2array(summary0)
  out[, 1] <- combine(fit.r = fit.r, fit.s = NULL, fit.m = fit.m, fit.y = fit.y, func = "smi")
  out[, 2] <- apply(out.full, 1, quantile, prob = (1 - conf.level)/2, na.rm = TRUE)
  out[, 3] <- apply(out.full, 1, quantile, prob = 1/2 + conf.level/2, na.rm = TRUE)
  all.out <- out.full

  if(long){
    out <- list(result = out, all.result = all.out)
  } else {
    out <- list(result = out)
  }
  return(out)

}
