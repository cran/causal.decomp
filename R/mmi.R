#' Multiple-Mediator-Imputation Estimation Method
#'
#' 'mmi' is used to estimate the initial disparity, disparity reduction, and
#' disparity remaining for causal decomposition analysis, using the
#' multiple-mediator-imputation estimation method proposed by Park et. al. (2020).
#'
#' @details This function returns the point estimates of the initial disparity,
#'   disparity reduction, and disparity remaining for a categorical
#'   treatment and a variety of types of outcome and mediator(s) in causal
#'   decomposition analysis. It also returns nonparametric
#'   percentile bootstrap confidence intervals for each estimate.
#'
#'   The initial disparity represents the expected difference in an outcome between
#'   a comparison group \eqn{R=j} and a reference group \eqn{R=i} where \eqn{i \neq j}{i != j}.
#'   That is,
#'   \deqn{\tau(i,j) \ = \ E\{Y|R=j\} - E\{Y|R=i\},}{%
#'         \delta(i,j) = E[Y|R=j] - E[Y|R=i],}
#'   where \eqn{R} and \eqn{Y} are the group indicator and the outcome variable, respectively.
#'   The disparity reduction represents the expected change in an outcome for the group \eqn{R=j}
#'   after adjusting the level of mediator(s) to the level of the reference group.
#'   That is,
#'   \deqn{\delta(j) \ = \ E\{Y|R=j\} - E\{Y(G_M(i))|R=j\},}{%
#'         \delta(j) = E[Y|R=j] - E[Y(G_M(i))|R=j],}
#'   where \eqn{G_M(i)}{G_M(i)} is a random draw from the mediator distribution of
#'   the reference group.
#'   The disparity remaining represents the remaining disparity for the group \eqn{R=j}
#'   even after adjusting the level of mediators to the reference group. Formally,
#'   \deqn{\zeta(i) \ = \ E\{Y(G_M(i))|R=j\} - E\{Y|R=i\}.}{%
#'         \zeta(i) = E[Y(G_M(i))|R=j] - E[Y|R=i].}
#'   The disparity reduction and remaining can be estimated using the
#'   multiple-mediator-imputation method suggested by Park et al. (2020).
#'   See the references for more details.
#'
#'   If one wants to make the inference conditional on baseline covariates,
#'   set 'conditional = T' and center the data before fitting the models.
#'
#'   As of version 0.0.1, the intetmediate confounder model ('fit.x') can be of
#'   class 'lm', 'glm', 'multinom', or 'polr', corresponding respectively to the
#'   linear regression models and generalized linear models, multinomial
#'   log-linear models, and ordered response models.
#'   The outcome model ('fit.y') can be of class 'lm' or 'glm'.
#'   Also, the treatment model ('fit.r') can be of class 'CBPS' or 'SumStat', both of
#'   which use the propensity score weighting. It is only necessary when 'conditional = F'.
#'
#' @param fit.r a fitted model object for treatment. Can be of class 'CBPS' or
#'   'SumStat'. Default is 'NULL'. Only necessary if 'conditional' is 'FALSE'.
#' @param fit.x a fitted model object for intermediate confounder(s). Each intermediate
#'   model can be of class 'lm', 'glm', 'multinom', or 'polr'. When multiple confounders
#'   are considered, can be of class 'list' containing multiple models.
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
#'   the cluster bootdtrap is used. Default is 'NULL'.
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
#' @seealso \code{\link{smi}}
#'
#' @references
#'   Park, S., Lee, C., and Qin, X. (2020). "Estimation and sensitivity analysis
#'   for causal decomposition in heath disparity research", arXiv preprint arXiv:2008.12812.
#'
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
#' fit.m1 <- lm(M.num ~ R + C.num + C.bin, data = sdata)
#' fit.m2 <- glm(M.bin ~ R + C.num + C.bin, data = sdata,
#'           family = binomial(link = "logit"))
#' require(MASS)
#' fit.m3 <- polr(M.cat ~ R + C.num + C.bin, data = sdata)
#' fit.s1 <- lm(S ~ R + C.num + C.bin, data = sdata)
#' require(nnet)
#' fit.m4 <- multinom(M.cat ~ R + C.num + C.bin, data = sdata)
#' fit.y1 <- lm(Y.num ~ R + M.num + M.bin + M.cat + S + C.num + C.bin,
#'           data = sdata)
#'
#' require(PSweight)
#' fit.r1 <- SumStat(R ~ C.num + C.bin, data = sdata, weight = "IPW")
#' \donttest{require(CBPS)
#' fit.r2 <- CBPS(R ~ C.num + C.bin, data = sdata, method = "exact",
#'           standardize = "TRUE")}
#'
#' res.1a <- mmi(fit.r = fit.r1, fit.x = fit.s1,
#'           fit.y = fit.y1, sims = 40, conditional = FALSE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.1a
#'
#' #-----------------------------------------------------------#
#' # Example 1-b: Binary Outcome
#' #-----------------------------------------------------------#
#' \donttest{fit.y2 <- glm(Y.bin ~ R + M.num + M.bin + M.cat + S + C.num + C.bin,
#'           data = sdata, family = binomial(link = "logit"))
#'
#' res.1b <- mmi(fit.r = fit.r1, fit.x = fit.s1,
#'           fit.y = fit.y2, sims = 40, conditional = FALSE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.1b}
#' 
#' #-----------------------------------------------------------#
#' # Example 2-a: Continuous Outcome, Conditional on Covariates
#' #-----------------------------------------------------------#
#' \donttest{# For conditional = TRUE, need to create data with centered covariates
#' # copy data
#' sdata.c <- sdata
#' # center quantitative covariate(s)
#' sdata.c$C.num <- scale(sdata.c$C.num, center = TRUE, scale = FALSE)
#' # center binary (or categorical) covariates(s)
#' # only neccessary if the desired baseline level is NOT the default baseline level.
#' sdata.c$C.bin <- relevel(sdata.c$C.bin, ref = "1")
#'
#' # fit mediator and outcome models
#' fit.m1 <- lm(M.num ~ R + C.num + C.bin, data = sdata.c)
#' fit.m2 <- glm(M.bin ~ R + C.num + C.bin, data = sdata.c,
#'           family = binomial(link = "logit"))
#' fit.m3 <- polr(M.cat ~ R + C.num + C.bin, data = sdata.c)
#' fit.s2 <- lm(S ~ R + C.num + C.bin, data = sdata.c)
#' fit.y1 <- lm(Y.num ~ R + M.num + M.bin + M.cat + S + C.num + C.bin,
#'           data = sdata.c)
#'
#' res.2a <- mmi(fit.x = fit.s2,
#'           fit.y = fit.y1, sims = 40, conditional = TRUE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.2a}
#' 
#' #-----------------------------------------------------------#
#' # Example 2-b: Binary Outcome, Conditional on Covariates
#' #-----------------------------------------------------------#
#' \donttest{fit.y2 <- glm(Y.bin ~ R + M.num + M.bin + M.cat + S + C.num + C.bin,
#'           data = sdata.c, family = binomial(link = "logit"))
#'
#' res.2b <- mmi(fit.x = fit.s2,
#'           fit.y = fit.y2, sims = 40, conditional = TRUE,
#'           covariates = c("C.num", "C.bin"), treat = "R")
#' res.2b}
mmi <- function(fit.r = NULL, fit.x, fit.y,
                treat, covariates, sims = 100, conf.level = .95,
                conditional = TRUE, cluster = NULL, long = FALSE,
                mc.cores = 1L){

  fit.m <- fit.x

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

  # Possibility of multiple mediators
  if(length(class(fit.m)) == 1 && class(fit.m) == "list"){
    isMultiConfounders <- TRUE
    num.ms <- length(fit.m)
  } else {
    isMultiConfounders <- FALSE
    num.ms <- 1
  }

  # Model type indicators
  isCBPS.r <- inherits(fit.r, "CBPS")
  isSumStat.r <- inherits(fit.r, "SumStat")

  isGlm.y <- inherits(fit.y, "glm")
  isLm.y <- inherits(fit.y, "lm")

  isGlm.m <- isLm.m <- isNominal.m <- isOrdinal.m <- FamilyM <- rep(NA, num.ms)
  for (mm in 1:num.ms) {
    if(isMultiConfounders){
      fit.mm <- fit.m[[mm]]
    } else if (!isMultiConfounders){
      fit.mm <- fit.m
    }

    isGlm.m[mm] <- inherits(fit.mm, "glm")
    isLm.m[mm] <- inherits(fit.mm, "lm")
    isNominal.m[mm] <- inherits(fit.mm, "multinom")
    isOrdinal.m[mm] <- inherits(fit.mm, "polr")

    if(!isGlm.m[mm] && !isLm.m[mm] && !isNominal.m[mm] && !isOrdinal.m[mm]){
      stop("unsupported mediator model(s)")
    }

    if (isGlm.m[mm]) {
      FamilyM[mm] <- fit.mm$family$family
    }
  }

  if(!is.null(fit.r) && !isCBPS.r && !isSumStat.r){
    stop("unsupported treatment model")
  }

  if(!isGlm.y && !isLm.y){
    stop("unsupported outcome model")
  }

  # Numbers of observations and model frame
  if(!is.null(fit.r) && isCBPS.r){
    n.r <- nrow(model.frame(fit.r))
  } else if (!is.null(fit.r) && isSumStat.r){
    n.r <- nrow(fit.r$propensity)
  } else {
    n.r <- NULL
  }
  n.m <- rep(NA, num.ms)
  for(mm in 1:num.ms){
    if(isMultiConfounders){
      fit.mm <- fit.m[[mm]]
    } else if (!isMultiConfounders){
      fit.mm <- fit.m
    }
    n.m[mm] <- nrow(model.frame(fit.mm))
  }
  if(length(unique(n.m)) > 1){
    stop("number of observations do not match between mediator models")
  } else {
    n.m <- n.m[1]
  }
  y.data <- model.frame(fit.y)
  n.y <- nrow(y.data)
  if(is.null(n.r)){
    if(n.m != n.y){
      stop("number of observations do not match between mediator and outcome models")
    }
  } else if (!is.null(n.r)){
    if(n.m != n.y | n.r != n.y | n.m != n.r){
      stop("number of observations do not match between treatment, mediator and outcome models")
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
                       cluster = cluster, func = "mmi", mc.cores = mc.cores)
  out.full <- simplify2array(summary0)
  out[, 1] <- combine(fit.r = fit.r, fit.s = NULL, fit.m = fit.m, fit.y = fit.y, func = "mmi")
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
