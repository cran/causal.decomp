# A function used in 'mmi' and 'smi' for comparing a pair of sel.lev.R vs ref.lev.treat
select.treat.pocr <- function(fit.s, fit.m, fit.y, sel.lev.treat, ref.lev.treat = 1, func){
  
  #:::::::::::::::::::::::::::#
  # Predict mediator::::::::::#
  #:::::::::::::::::::::::::::#
  y.data.new <- y.data
  y.data.new[, treat] <- levels(y.data[, treat])[sel.lev.treat]
  
  # Compute predicted values of mediator
  
  ### Case 1: GLM Mediator
  if(isGlm.m){
    
    if(FamilyM == "poisson" | FamilyM == "Gamma" |
       FamilyM == "gaussian" | FamilyM == "inverse.gaussian"){
      
      mod.m.f <- as.formula(paste(M, "~ ", paste(covariates, collapse = "+")))
      data.ref <- subset(y.data, y.data[, treat] == levels(y.data[, treat])[ref.lev.treat])
      mod.m.ref <- glm(mod.m.f, data = data.ref, family = M.fun)
      data.sel <- subset(y.data, y.data[, treat] == levels(y.data[, treat])[sel.lev.treat])
      mod.m.sel <- glm(mod.m.f, data = data.sel, family = M.fun)
      
      # alpha1
      if(LinkM == "identity") {
        alpha1 <- coef(mod.m.ref)[1] - coef(mod.m.sel)[1]
      } else if(LinkM == "log") {
        alpha1 <- exp(coef(mod.m.ref)[1]) - exp(coef(mod.m.sel)[1])
      } else if(LinkM == "sqrt") {
        alpha1 <- (coef(mod.m.ref)[1])^2 - (coef(mod.m.sel)[1])^2
      } else if(LinkM == "inverse") {
        alpha1 <- (coef(mod.m.ref)[1])^(-1) - (coef(mod.m.sel)[1])^(-1)
      } else if(LinkM == "1/mu^2") {
        alpha1 <- (coef(mod.m.ref)[1])^(-1/2) - (coef(mod.m.sel)[1])^(-1/2)
      } 
      
      # beta3, beta4
      name.beta3 <- paste(M, levels(y.data[, M])[2], sep = "")
      name.beta41 <- paste(paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = ""),
                           paste(M, levels(y.data[, M])[2], sep = ""), sep = ":")
      name.beta42 <- paste(paste(M, levels(y.data[, M])[2], sep = ""),
                           paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = ""), sep = ":")
      
      ind.beta3 <- which(names(coef(fit.y)) == name.beta3)
      ind.beta4 <- which(names(coef(fit.y)) == name.beta41 | names(coef(fit.y)) == name.beta42)
      
      beta3 <- coef(fit.y)[ind.beta3]
      beta4 <- coef(fit.y)[ind.beta4]
      
      # phi1
      mod.y.f <- as.formula(paste(outcome, "~ ", paste(c(treat, covariates), collapse = "+")))
      mod.y <- lm(mod.y.f, data = y.data)
      name.phi1 <- paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = "")
      ind.phi1 <- which(names(coef(mod.y)) == name.phi1)
      coef.phi1 <- coef(mod.y)[ind.phi1]
      
      # disparity reduction
      disp.red <- alpha1 * (beta3 + beta4)
      # disparity remaining for a continuous mediator
      disp.rem <- phi1 - alpha1 * (beta3 + beta4)
      
    } else if (FamilyM == "binomial"){
      
      mod.m.f <- as.formula(paste(M, "~ ", paste(covariates, collapse = "+")))
      data.ref <- subset(y.data, y.data[, treat] == levels(y.data[, treat])[ref.lev.treat])
      mod.m.ref <- glm(mod.m.f, data = data.ref, family = M.fun)
      data.sel <- subset(y.data, y.data[, treat] == levels(y.data[, treat])[sel.lev.treat])
      mod.m.sel <- glm(mod.m.f, data = data.sel, family = M.fun)
      
      # alpha1, alpha0
      if(LinkM == "logit"){
        alpha1 <- plogis(coef(mod.m.ref)[1]) - plogis(coef(mod.m.sel)[1])
        alpha0 <- plogis(coef(mod.m.sel)[1])
      } else if(LinkM == "probit") {
        alpha1 <- pnorm(coef(mod.m.ref)[1]) - pnorm(coef(mod.m.sel)[1])
        alpha0 <- pnorm(coef(mod.m.sel)[1])
      } else if(LinkM == "cloglog") {
        alpha1 <- exp(- exp(coef(mod.m.sel)[1])) - exp(- exp(coef(mod.m.ref)[1]))
        alpha0 <- 1 - exp(- exp(coef(mod.m.sel)[1]))
      }
      
      # beta1, beta2, beta3, beta4
      name.beta1 <- paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = "")
      name.beta2 <- S
      name.beta3 <- paste(M, levels(y.data[, M])[2], sep = "")
      name.beta41 <- paste(paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = ""),
                           paste(M, levels(y.data[, M])[2], sep = ""), sep = ":")
      name.beta42 <- paste(paste(M, levels(y.data[, M])[2], sep = ""),
                           paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = ""), sep = ":")
      
      ind.beta1 <- which(names(coef(fit.y)) == name.beta1)
      ind.beta2 <- which(names(coef(fit.y)) == name.beta2)
      ind.beta3 <- which(names(coef(fit.y)) == name.beta3)
      ind.beta4 <- which(names(coef(fit.y)) == name.beta41 | names(coef(fit.y)) == name.beta42)
      
      beta1 <- coef(fit.y)[ind.beta1]
      beta2 <- coef(fit.y)[ind.beta2]
      beta3 <- coef(fit.y)[ind.beta3]
      beta4 <- coef(fit.y)[ind.beta4]
      
      # gamma1
      name.gamma1 <- paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = "")
      ind.gamma1 <- which(names(coef(fit.s)) == name.gamma1)
      gamma1 <- coef(fit.s)[ind.gamma1]
      
      # disparity reduction
      disp.red <- alpha1 * (beta3 + beta4)
      # disparity remaining for a binary mediator
      disp.rem <- beta1 + beta2 * gamma1 + beta4 * alpha0
      
    } else {
      stop("unsupported glm family")
    }
    
    ### Case 2: LM Mediator
  } else if(isLm.m & !isGlm.m){
    
    # alpha1
    name.R.in.M <- paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = "")
    ind.R.in.M <- which(names(coef(fit.m)) == name.R.in.M)
    coef.R.in.M <- coef(fit.m)[ind.R.in.M] 
    
    # beta3
    ind.M.in.Y <- which(names(coef(fit.y)) == M)
    coef.M.in.Y <- coef(fit.y)[ind.M.in.Y]
    
    # beta4
    name.RM.in.Y1 <- paste(paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = ""), M, sep = ":")
    name.RM.in.Y2 <- paste(M, paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = ""), sep = ":")
    ind.RM.in.Y <- which(names(coef(fit.y)) == name.RM.in.Y1 | names(coef(fit.y)) == name.RM.in.Y2)
    coef.RM.in.Y <- coef(fit.y)[ind.RM.in.Y]
    
    # phi1
    mod.y.f <- as.formula(paste(outcome, "~ ", paste(c(treat, covariates), collapse = "+")))
    mod.y <- lm(mod.y.f, data = y.data)
    name.R.in.y <- paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = "")
    ind.R.in.y <- which(names(coef(mod.y)) == name.R.in.y)
    coef.R.in.y <- coef(mod.y)[ind.R.in.y]
    
    # disparity reduction = alpha1 * (beta3 + beta4)
    disp.red <- coef.R.in.M * (coef.M.in.Y + coef.RM.in.Y)
    # disparity remaining = phi1 - alpha1 * (beta3 + beta4)
    disp.rem <- coef.R.in.y - coef.R.in.M * (coef.M.in.Y + coef.RM.in.Y)
    
    ### Case 3: Nominal or Ordinal Mediator
  } else if (isNominal.m | isOrdinal.m) {
    
    mod.m.f <- as.formula(paste(M, "~ ", paste(covariates, collapse = "+")))
    m <- length(levels(y.data[, M]))
    prob.m.ref <- prob.m.sel <- rep(NA, m)
    
    if(isNominal.m){
      
      data.ref <- subset(y.data, y.data[, treat] == levels(y.data[, treat])[ref.lev.treat])
      mod.m.ref <- multinom(mod.m.f, data = data.ref, trace = FALSE)
      data.sel <- subset(y.data, y.data[, treat] == levels(y.data[, treat])[sel.lev.treat])
      mod.m.sel <- multinom(mod.m.f, data = data.sel, trace = FALSE)
      
      prob.m.ref[1] <- 1
      prob.m.sel[1] <- 1
      for(mm in 2:m){
        prob.m.ref[mm] <- exp(coef(mod.m.ref)[mm - 1, 1])
        prob.m.sel[mm] <- exp(coef(mod.m.sel)[mm - 1, 1])
      }
      prob.m.ref <- prob.m.ref/(1 + sum(prob.m.ref[2:m]))
      prob.m.sel <- prob.m.sel/(1 + sum(prob.m.sel[2:m]))
      
    } else if (isOrdinal.m){
      
      data.ref <- subset(y.data, y.data[, treat] == levels(y.data[, treat])[ref.lev.treat])
      mod.m.ref <- polr(mod.m.f, data = data.ref, method = fit.m$method)
      data.sel <- subset(y.data, y.data[, treat] == levels(y.data[, treat])[sel.lev.treat])
      mod.m.sel <- polr(mod.m.f, data = data.sel, method = fit.m$method)
      
      prob.m.ref[1] <- pfun(mod.m.ref$zeta[1])
      prob.m.sel[1] <- pfun(mod.m.sel$zeta[1])
      for(mm in 2:(m - 1)){
        prob.m.ref[mm] <- pfun(mod.m.ref$zeta[mm]) - pfun(mod.m.ref$zeta[mm - 1])
        prob.m.sel[mm] <- pfun(mod.m.sel$zeta[mm]) - pfun(mod.m.sel$zeta[mm - 1])
      }
      prob.m.ref[m] <- 1 - pfun(mod.m.ref$zeta[m - 1])
      prob.m.sel[m] <- 1 - pfun(mod.m.sel$zeta[m - 1])
      
    }
    #prob.m.ref <- colMeans(predict(mod.m.ref, type = "probs"))
    #prob.m.sel <- colMeans(predict(mod.m.sel, type = "probs"))
    
    # gamma1
    name.R <- paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = "")
    ind.R <- which(names(coef(fit.y)) == name.R)
    ind.S <- which(names(coef(fit.y)) == S)
    ind.R.in.S <- which(names(coef(fit.s)) == name.R)
    
    disp.red <- disp.rem <- rep(NA, m)
    for(mm in 2:m){
      name.M <- paste(M, levels(y.data[, M])[mm], sep = "")
      name.RM1 <- paste(name.R, name.M, sep = ":")
      name.RM2 <- paste(name.M, name.R, sep = ":")
      ind.M <- which(names(coef(fit.y)) == name.M)
      ind.RM <- which(names(coef(fit.y)) == name.RM1 | names(coef(fit.y)) == name.RM2)
      
      # beta1, beta2, beta4
      name.beta1 <- paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = "")
      name.beta2 <- S
      name.beta41 <- paste(paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = ""),
                           paste(M, levels(y.data[, M])[2], sep = ""), sep = ":")
      name.beta42 <- paste(paste(M, levels(y.data[, M])[2], sep = ""),
                           paste(treat, levels(y.data[, treat])[sel.lev.treat], sep = ""), sep = ":")
      
      ind.beta1 <- which(names(coef(fit.y)) == name.beta1)
      ind.beta2 <- which(names(coef(fit.y)) == name.beta2)
      ind.beta4 <- which(names(coef(fit.y)) == name.beta41 | names(coef(fit.y)) == name.beta42)
      
      beta1 <- coef(fit.y)[ind.beta1]
      beta2 <- coef(fit.y)[ind.beta2]
      beta4 <- coef(fit.y)[ind.beta4]
      
      disp.red[mm] <- (coef(fit.y)[ind.M] + coef(fit.y)[ind.RM]) * (prob.m.sel[mm] - prob.m.ref[mm])
      disp.rem[mm] <- coef(fit.y)[ind.RM] * prob.m.ref[mm]
    }
    
    disp.red <- sum(disp.red, na.rm = TRUE)
    disp.rem <- coef(fit.y)[ind.R] + coef(fit.y)[ind.S] * coef(fit.s)[ind.R.in.S] + sum(disp.rem, na.rm = TRUE)
    
  } else {
    stop("mediator model is not yet implemented")
  }
  
  #:::::::::::::::#
  #:Results:::::::#
  #:::::::::::::::#
  # Initial disparity, Disparity remaining, and Disparity reduction in order
  out <- rep(NA, 3)
  out[2] <- disp.rem
  out[3] <- disp.red
  out[1] <- out[2] + out[3]
  
  return(out)
}