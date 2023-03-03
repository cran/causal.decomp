# A function used in 'mmi', 'smi', and 'pocr' for comparing all pairs of R
combine <- function(fit.r, fit.x, fit.m, fit.y, func, weights){
  
  # Model frame Y model
  y.data <- model.frame(fit.y)
  # Number of observations
  n.y <- nrow(y.data)
  
  if(func == "mmi" | func == "smi"){
    ## Calculate weighted outcome from R model
    if(!is.null(fit.r) && isCBPS.r){
      y.data$w <- fit.r$weights
      prop.treat <- table(y.data[, treat])/nrow(y.data)
    } else if(!is.null(fit.r) && isSumStat.r){
      y.data$w <- fit.r$ps.weights[, 3]
      prop.treat <- table(y.data[, treat])/nrow(y.data)
    } else {
      y.data$w <- rep(1, n.y)
      prop.treat <- rep(1, length(levels(y.data[, treat])))
    }
    
    wy <- rep(NA, length(levels(y.data[, treat])))
    if(isGlm.y){ ## fixed 8/27/2021 for conditional = T
      y.data$pred.prob <- predict(fit.y, newdata = y.data, type = "response")
      if(conditional){
        wy.f <- as.formula(paste("pred.prob", paste(covariates, collapse = "+"), sep = "~"))
      } else if (!conditional){
        wy.f <- as.formula(paste("pred.prob", "~ 1"))
      }
    } else if (!isGlm.y){
      if(conditional){
        wy.f <- as.formula(paste(outcome, paste(covariates, collapse = "+"), sep = "~"))
      } else if (!conditional){
        wy.f <- as.formula(paste(outcome, "~ 1"))
      }
    }
    for(i in 1:length(levels(y.data[, treat]))){
      wy[i] <- lm(wy.f, weights = w * prop.treat[i],
                  data = y.data[y.data[, treat] == levels(y.data[, treat])[i], ])$coef[1]
    }
  }
  
  if(func == "mmi" | func == "smi"){
    environment(select.treat) <- environment()
    # Do all combinations, ref.lev.treat vs. sel.lev.treat
    res <- NULL
    for(l in 2:length(levels(y.data[, treat]))){
      res <- c(res, select.treat(fit.m = fit.m, fit.y = fit.y, sel.lev.treat = l,
                                     ref.lev.treat = 1,
                                     prop.treat = prop.treat, wy = wy, func = func, weights = weights))
    }
  } else if (func == "pocr"){
    environment(select.treat.pocr) <- environment()
    # Do all combinations, ref.lev.treat vs. sel.lev.treat
    res <- NULL
    for(l in 2:length(levels(y.data[, treat]))){
      res <- c(res, select.treat.pocr(fit.x = fit.x, fit.m = fit.m, fit.y = fit.y, sel.lev.treat = l,
                                      ref.lev.treat = 1))
    }
  }
  
  return(res) 
}
