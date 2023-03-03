#' Sensitivity Analysis Using R-Squared Values for Causal Decomposition Analysis
#'
#' The function 'sensitivity' is used to implement sensitivity analysis for the
#' causal decomposition analysis, using R-squared values as sensitivity parameters.
#' 
#' @usage 
#' sensitivity(boot.res, fit.y, fit.m = NULL, mediator = NULL, covariates, treat,
#'             sel.lev.treat, conf.level = 0.95, max.rsq = 0.3)
#'
#' @details This function is used to implement sensitivity analysis for the
#' causal decomposition analysis, using two sensitivity parameters: (i) the
#' R-squared value between the outcome and unobserved confounder given the treatment,
#' intermediate confounder, mediator, and covariates; and (ii) the R-squared value
#' between the mediator and unobserved confounder given the treatment, intermediate
#' confounder, and covariates (Park et al., 2023).
#'
#' As of version 0.1.0, 'boot.res' must be fitted by the 'smi' function with a single
#' mediator, which can be of class 'lm', 'glm', or 'polr'. 
#'
#' @param boot.res bootstrap results from an object fitted by the 'smi' function.
#' @param fit.y outcome model used in fitting the 'smi' function.
#'   Can be of class 'lm' or 'glm'.
#' @param fit.m mediator model used in fitting the 'smi' function.
#'   Can be of class 'lm', 'glm', or 'polr'.
#' @param mediator a vector containing the name of mediator used in the models.
#' @param covariates a vector containing the name of the covariate variable(s)
#'   used in the models. Each covariate can be categorical with two or more
#'   categories (two- or multi-valued factor) or continuous (numeric).
#' @param treat a character string indicating the name of the treatment variable
#'   used in the models. The treatment can be categorical with two or more
#'   categories (two- or multi-valued factor).
#' @param sel.lev.treat a level of categorical treatment variable which is to be
#'   compared with the reference level.
#' @param conf.level level of the returned two-sided confidence intervals,
#'   which are estimated by the nonparametric percentile bootstrap method.
#'   Default is .95, which returns the 2.5 and 97.5 percentiles of the simulated
#'   quantities.
#' @param max.rsq upper limit of the two sensitivity parameters (R-squared values).
#'   Once it is set, the R-squared values between 0 and this upper limit are
#'   explored to draw the sensitivity contour plots. Default is 0.3.
#'
#' @return
#'
#'   \item{call}{original function call.}
#'   \item{disparity.reduction}{a matrix containing the estimated disparity reduction
#'   value along with lower and upper limits of the confidence interval, for each
#'   combination of the two sensitivity parameters, assuming those two are equal.}
#'   \item{disparity.remaining}{a matrix containing the estimated disparity remaining
#'   value along with lower and upper limits of the confidence interval, for each
#'   combination of the two sensitivity parameters, assuming those two are equal.}
#'   \item{rm}{R-squared values between the mediator and unobserved confounder
#'   (first sensitivity parameter), which are explored for the sensitivity analysis.}
#'   \item{ry}{R-squared values between the outcome and unobserved confounder,
#'   (second sensitivity parameter), which are explored for the sensitivity analysis.}
#'   \item{red}{a matrix containing the estimated disparity reduction values given
#'   each combination of the two sensitivity parameters.}
#'   \item{lower_red}{a matrix containing the lower limit of disparity reduction given
#'   each combination of the two sensitivity parameters.}
#'   \item{rem}{a matrix containing the estimated disparity remaining values given
#'   each combination of the two sensitivity parameters.}
#'   \item{lower_rem}{a matrix containing the lower limit of disparity remaining given
#'   each combination of the two sensitivity parameters.}
#'   \item{RV_red}{robustness value for disparity reduction, which represents the
#'   strength of association that will explain away the estimated disparity reduction.}
#'   \item{RV_red_alpha}{robustness value for disparity reduction, which represents
#'   the strength of association that will change the significance of the estimated
#'   disparity reduction at the given significance level, assuming an equal association
#'   to the mediator and the outcome.}
#'   \item{RV_rem}{robustness value for disparity remaining, which represents the
#'   strength of association that will explain away the estimated disparity remaining.}
#'   \item{RV_rem_alpha}{robustness value for disparity remaining, which represents
#'   the strength of association that will change the significance of the estimated
#'   disparity remaining at the given significance level, assuming an equal association
#'   to the mediator and the outcome.}
#'
#' @author
#'   Suyeon Kang, University of California, Riverside, \email{skang062@@ucr.edu};
#'   Soojin Park, University of California, Riverside, \email{soojinp@@ucr.edu}.
#'
#' @seealso \code{\link{smi}}
#'
#' @references
#'   Park, S., Kang, S., Lee, C., & Ma, S. (2023). Sensitivity analysis for causal
#'   decomposition analysis: Assessing robustness toward omitted variable bias,
#'   Journal of Causal Inference. Forthcoming.
#'
#' @export
#' @examples
#' data(sMIDUS)
#' 
#' # Center covariates
#' sMIDUS$age <- scale(sMIDUS$age, center = TRUE, scale = FALSE)
#' sMIDUS$stroke <- scale(sMIDUS$stroke, center = TRUE, scale = FALSE)
#' sMIDUS$T2DM <- scale(sMIDUS$T2DM, center = TRUE, scale = FALSE)
#' sMIDUS$heart <- scale(sMIDUS$heart, center = TRUE, scale = FALSE)
#' 
#' fit.m <- lm(edu ~ racesex + age + stroke + T2DM + heart, data = sMIDUS)
#' fit.y <- lm(health ~ racesex + lowchildSES + abuse + edu + racesex:edu +
#'             age + stroke + T2DM + heart, data = sMIDUS)
#' smiRes <- smi(fit.m = fit.m, fit.y = fit.y, sims = 40, conf.level = .95,
#'           covariates = c("age", "stroke", "T2DM", "heart"), treat = "racesex", seed = 227)
#' sensRes <- sensitivity(boot.res = smiRes, fit.m = fit.m, fit.y = fit.y, mediator = "edu",
#'                        covariates = c("age", "stroke", "T2DM", "heart"), treat = "racesex",
#'                        sel.lev.treat = "4", max.rsq = 0.3)
#' 
#' sensRes
#' names(sensRes) 
#' # Create sensitivity contour plots
#' plot(sensRes)
sensitivity <- function(boot.res, fit.y, fit.m = NULL, mediator = NULL, covariates,
                        treat, sel.lev.treat, conf.level = 0.95, max.rsq = 0.3){
  
  # Match arguments
  call <- match.call()
  
  # warning for inappropriate setting (version 0.1.0)
  if(!inherits(boot.res, "smi")){
    stop("class of 'boot.res' must be 'smi'.")
  }
  
  if(!is.null(fit.m)){
    isNominal.m <- inherits(fit.m, "multinom")
    if(isNominal.m){
      stop("class of 'fit.m' must be 'lm', 'glm'm or 'polr'.")
    }
  }
  
  data <- model.frame(fit.y)
  if(!(sel.lev.treat %in% levels(data[, treat]))){
    stop("'sel.lev.treat' must be one of the levels of treatment")
  }
  if(sel.lev.treat == levels(data[, treat])[1]){
    stop("'sel.lev.treat' must not be the reference level of treatment")
  }
  
  outcome <- all.vars(formula(fit.y))[[1]]
  len.treat <- length(levels(data[, treat]))
  # new data with releveled treatment
  data.new <- data
  data.new[, treat] <- relevel(data.new[, treat], ref = sel.lev.treat)
  # new outcome model with releveled treatment
  fit.y.new <- update(fit.y, data = data.new)
  lev <- which(levels(data[, treat]) == sel.lev.treat)
  
  if(is.null(fit.m) & !is.null(mediator)){
    
    # Possibility of multiple mediators
    if(inherits(fit.m, "list")){
      isMultiMediators <- TRUE
      num.meds <- length(fit.m)
    } else {
      isMultiMediators <- FALSE
      num.meds <- 1
    }
    ## Fit mediator model(s)
    # make an empty list to save mediator model(s)
    fit.meds <- vector(mode = "list", length = num.meds)
    for (i in 1:num.meds) {
      # The formula for each mediator model: lm(M1 ~ R + C); lm(M2 ~ R + C);...
      med.formula <- as.formula(paste(mediator[i], paste(c(treat, covariates), collapse = " + "), sep = " ~ "))
      if(is.factor(data[, mediator[i]])){
        fit.meds[[i]] <- lm(med.formula, data = data)
      } else {
        fit.meds[[i]] <- glm(med.formula, data = data)
      }
    }
    #fit.meds.new1 <- update(fit.meds[[1]], data = data.new)
    #vcov(fit.meds.new1)["racesex41", "racesex41"]
    
  } else if (!is.null(fit.m) & is.null(mediator)){
    
    # Possibility of multiple mediators
    if(inherits(fit.m, "list")){
      isMultiMediators <- TRUE
      num.meds <- length(fit.m)
    } else {
      isMultiMediators <- FALSE
      num.meds <- 1
    }
    
    fit.meds <- vector(mode = "list", length = num.meds)
    mediator <- rep(NA, num.meds)
    for (i in 1:num.meds) {
      fit.meds[[i]] <- fit.m[[i]]
      mediators[i] <- all.vars(formula(fit.m[[i]]))[[1]]
    }
    
  } else if (!is.null(fit.m) & !is.null(mediator)) {
    
    # Possibility of multiple mediators
    if(inherits(fit.m, "list")){
      isMultiMediators <- TRUE
      num.meds <- length(fit.m)
    } else {
      isMultiMediators <- FALSE
      num.meds <- 1
    }
    
    fit.meds <- vector(mode = "list", length = num.meds)
    mediator0 <- rep(NA, num.meds)
    for (i in 1:num.meds) {
      if(!isMultiMediators){
        fit.meds[[i]] <- fit.m
        mediator0[i] <- all.vars(formula(fit.m))[[1]]
      } else {
        fit.meds[[i]] <- fit.m[[i]]
        mediator0[i] <- all.vars(formula(fit.m[[i]]))[[1]]
      }
      if(mediator0[i] != mediator[i]){
        stop("Response variable of 'fit.m' and 'mediator' must be match.")
      }
    }
    
  } else if (is.null(fit.m) & is.null(mediator)) {
    stop("Either 'fit.m' or 'mediator' must not be NULL.")
  }
  
  ##1. Calculate the SE of gamma_dm and df from the coefficient of M for a comparison group in fit.y.
  if (num.meds == 1) {
    meds <- all.vars(formula(fit.m))[[1]]
    if(is.factor(data[, mediator])){
      meds <- which(meds == substring(colnames(vcov(fit.y.new)), 1, nchar(meds)))[1]
    }
    var_gamma <- vcov(fit.y.new)[meds, meds]
    beta.resm.est <- coef(fit.y.new)[meds]
  } else if (num.meds == 2) {
    meds1 <- all.vars(formula(fit.m[[1]]))[[1]]
    meds2 <- all.vars(formula(fit.m[[2]]))[[1]]
    if(is.factor(data[, mediator[1]])){
      meds1 <- which(meds1 == substring(colnames(vcov(fit.y.new)), 1, nchar(meds1)))[1]
    }
    if(is.factor(data[, mediator[2]])){
      meds2 <- which(meds2 == substring(colnames(vcov(fit.y.new)), 1, nchar(meds2)))[1]
    }
    var_gamma <- vcov(fit.y.new)[meds1, meds1] + vcov(fit.y.new)[meds2, meds2] +
      2 * vcov(fit.y.new)[meds1, meds2]
    beta.resm.est <- sum(coef(fit.y.new)[c(meds1, meds2)]) ##??? when there are two Ms.
  } else if (num.meds > 2) {
    # will be added.
    stop("The case of three or more mediators is not supported yet.")
  }
  se_gamma <- sqrt(var_gamma)
  df <- fit.y.new$df.residual
  
  ##2. Calculate the effect of R on D and M (the coefficient for R in fit.m).
  treat.lev <- paste(treat, sel.lev.treat, sep = "")
  # sample covariance of initial disparity and disparity reduction
  cov.ini.red <- cov(boot.res$all.result[3 * (lev - 1) - 2, ],
                     boot.res$all.result[3 * (lev - 1), ])
  # sample covariance of initial disparity and alpha_r
  cov.ini.alphar <- cov(boot.res$all.result[3 * (lev - 1) - 2, ],
                        boot.res$alpha.r[, lev - 1])
  # sample covariance of initial disparity and se_gamma
  cov.ini.segamma <- cov(boot.res$all.result[3 * (lev - 1) - 2, ],
                         boot.res$se.gamma[, lev - 1])
  # sample variance of alpha_r
  var_alphahat.r <- var(boot.res$alpha.r[, lev - 1]) # 09/05/2022, 2/28/2023
  
  if (num.meds == 1) {
    alpha.r <- coef(fit.meds[[1]])[treat.lev] # coef of R[lev] on M = abs(alpha_r)
    ab <- abs(alpha.r)
    var_ab <- vcov(fit.meds[[1]])[treat.lev, treat.lev]
  } else if (num.meds == 2) {
    alpha.r <- coef(fit.meds[[1]])[treat.lev] + coef(fit.meds[[2]])[treat.lev]  # 09/05/2022
    ab <- abs(alpha.r)
    var_ab <- vcov(fit.meds[[1]])[treat.lev, treat.lev] +
      vcov(fit.meds[[2]])[treat.lev, treat.lev] - 2 * cov.ini.red
  } else if (num.meds > 2) {
    # will be added.
    stop("The case of three or more mediators is not supported yet.")
  }
  
  ## Calculate the variance of initial disparity # Use the one from the estimator function. 
  # Y ~ R + C
  ini.formula <- as.formula(paste(outcome, paste(c(treat, covariates), collapse = " + "), sep = " ~ "))
  ini <- lm(ini.formula, data = data)
  #var_tau <- vcov(ini)[treat.lev, treat.lev]
  var_tau <- var(boot.res$all.result[3 * (lev - 1) - 2, ])  # added 09/05/2022
  
  ##3. The rest of the R squared values (how much unobserved confounders explain the mediator and outcome)
  #   are sensitivity parameters.  
  # Set sensitivity parameters
  ry <- seq(0, max.rsq, by = 0.01) 
  rm <- seq(0, max.rsq, by = 0.01) 
  
  # Point Estimate of Disparity Reduction for ref vs lev
  res.red <- boot.res$result[3 * (lev - 1), 1]
  # Point Estimate of Disparity Remaining for ref vs lev
  res.rem <- boot.res$result[3 * (lev - 1) - 1, 1]
  red <- rem <- lower_red <- lower_rem <- upper_red <- upper_rem <- var_ngamma <- matrix(NA, length(ry), length(rm))
  
  for (i in 1:length(ry)){
    for (j in 1:length(rm)){
      
      # true disp reduction (= est - bias)
      red[i, j] <- abs(res.red) - ab * se_gamma * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j])
      beta.m <- beta.resm.est + ifelse(beta.resm.est > 0, - 1, 1) * se_gamma * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j])
      # new variance of gamma
      var_ngamma[i, j] <- var_gamma * (1 - ry[i])/(1 - rm[j]) * (df/(df - 1)) 
      # new CI for disp reduction
      qn <- qnorm(1/2 + conf.level/2)
      se_red <- sqrt(ab^2 * var_ngamma[i, j] + beta.m^2 * var_alphahat.r)  # Eq (15)
      lower_red[i, j] <- red[i, j] - qn * se_red
      upper_red[i, j] <- red[i, j] + qn * se_red
      
      # true disp remaining
      rem[i, j] <- abs(res.rem) - ab * se_gamma * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j])
      # new CI for disp remaining
      qn <- qnorm(1/2 + conf.level/2)
      k <- 1 * sqrt(ry[i] * rm[j] * df)/sqrt(1 - rm[j])
      se_rem <- sqrt(var_tau + se_red^2 - 2 * cov.ini.red +
                       2 * k * mean(boot.res$se.gamma[, lev - 1]) * cov.ini.alphar + # se_gamma
                       2 * k * mean(boot.res$alpha.r[, lev - 1]) * cov.ini.segamma)
      lower_rem[i, j] <- rem[i, j] - qn * se_rem
      upper_rem[i, j] <- rem[i, j] + qn * se_rem
      
    }
  }
  
  out.red <- cbind(ry, rm, diag(red), diag(lower_red), diag(upper_red))
  colnames(out.red) <- c("r2.y", "r2.m", "Reduction", "95% CI Lower", "95% CI Upper")
  out.rem <- cbind(ry, rm, diag(rem), diag(lower_rem), diag(upper_rem))
  colnames(out.rem) <- c("r2.y", "r2.m", "Remaining", "95% CI Lower", "95% CI Upper")
  
  ### ROBUSTNESS VALUE (RV)
  g_red <- abs(res.red)/(ab * se_gamma * sqrt(df))
  g_rem <- abs(res.rem)/(ab * se_gamma * sqrt(df))
  RV_red <- as.numeric((1/2) * (sqrt(g_red^4 + 4 * g_red^2) - g_red^2))
  RV_rem <- as.numeric((1/2) * (sqrt(g_rem^4 + 4 * g_rem^2) - g_rem^2))
  
  #ii <- which(lower_red < 0 & upper_red > 0, arr.ind = TRUE)
  ii <- which(abs(red - 0) < 0.01, arr.ind = TRUE)
  ry1 <- ry[ii[, 1]]
  rm1 <- rm[ii[, 2]]
  RV_red_alpha <- mean(sqrt(rm1 * ry1)) # this is the minimum R squared of unobserved confounder 
  
  #kk <- which(lower_rem < 0 & upper_rem > 0, arr.ind = TRUE)
  kk <- which(abs(rem - 0) < 0.01, arr.ind = TRUE)
  ry1 <- ry[kk[, 1]]
  rm1 <- rm[kk[, 2]]
  RV_rem_alpha <- mean(sqrt(rm1 * ry1)) # this is the minimum R squared of unobserved confounder 
  
  output <- list(call = call, disparity.reduction = out.red, disparity.remaining = out.rem,
                 rm = rm, ry = ry, red = red, lower_red = lower_red, rem = rem, lower_rem = lower_rem,
                 RV_red = RV_red, RV_red_alpha = RV_red_alpha,
                 RV_rem = RV_rem, RV_rem_alpha = RV_rem_alpha)
  
  class(output) <- "sensitivity"
  return(output)
  
}
