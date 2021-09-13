# A function used in 'mmi' and 'smi' for comparing a pair of sel.lev.R vs ref.lev.treat
select.treat <- function(fit.m, fit.y, sel.lev.treat, ref.lev.treat = 1, prop.treat, wy, func){
  
  #:::::::::::::::::::::::::::#
  # Predict mediator::::::::::#
  #:::::::::::::::::::::::::::#
  if(func == "mmi"){
    y.data.new <- y.data
    y.data.new[, treat] <- levels(y.data[, treat])[sel.lev.treat]
    y.data.new2 <- y.data.new
  } else if(func == "smi"){
    y.data.new <- y.data.new2 <- y.data 
    y.data.new[, treat] <- levels(y.data[, treat])[1]
  }
  
  # Compute predicted values of mediator
  M <- rep(NA, num.ms)
  PredictM <- data.frame(matrix(NA, n.y, num.ms))
  
  for(i in 1:num.ms){
    
    if(isMultiConfounders){
      fit.mm <- fit.m[[i]]
    } else if (!isMultiConfounders){
      fit.mm <- fit.m
    }
    M[i] <- all.vars(formula(fit.mm))[[1]]
    
    ### Case 1: GLM Mediator
    if(isGlm.m[i]){
      
      muM <- predict(fit.mm, newdata = y.data.new, type = "response")
      
      if(FamilyM[i] == "poisson"){
        PredictM[, i] <- rpois(n.y, lambda = muM)
      } else if (FamilyM[i] == "Gamma") {
        shape <- gamma.shape(fit.mm)$alpha
        PredictM[, i] <- rgamma(n.y, shape = shape, scale = muM/shape)
      } else if (FamilyM[i] == "binomial"){
        ind.lev <- rbinom(n.y, size = 1, prob = muM) + 1
        PredictM[, i] <- sapply(ind.lev, FUN = function(ind){levels(y.data.new[, M[i]])[ind]})
        PredictM[, i] <- as.factor(PredictM[, i])
      } else if (FamilyM[i] == "gaussian"){
        sigma <- sqrt(summary(fit.mm)$dispersion)
        error <- rnorm(n.y, mean = 0, sd = sigma)
        PredictM[, i] <- muM + error
      } else if (FamilyM[i] == "inverse.gaussian"){
        disp <- summary(fit.mm)$dispersion
        PredictM[, i] <- SuppDists::rinvGauss(n.y, nu = muM, lambda = 1/disp)
      } else {
        stop("unsupported glm family")
      }
      
      ### Case 2: LM Mediator
    } else if(isLm.m[i] & !isGlm.m[i]){
      
      sigma <- summary(fit.mm)$sigma
      error <- rnorm(n.y, mean = 0, sd = sigma)
      PredictM[, i] <- predict(fit.mm, type = "response", newdata = y.data.new) + error
      
      ### Case 3: Nominal or Ordinal Mediator
    } else if (isNominal.m[i] | isOrdinal.m[i]) {
      
      probs <- predict(fit.mm, newdata = y.data.new, type = "probs")
      m <- length(unique(y.data[, M[i]]))
      draws <- matrix(NA, n.y, m)
      for(ii in 1:n.y){
        draws[ii, ] <- t(rmultinom(1, 1, prob = probs[ii, ]))
      }
      ind.lev <- apply(draws, 1, which.max)
      PredictM[, i] <- sapply(ind.lev, FUN = function(ind){levels(y.data.new[, M[i]])[ind]})
      PredictM[, i] <- as.factor(PredictM[, i])
      
    } else {
      stop("mediator model(s) is(are) not yet implemented")
    }
    
  }
  
  #::::::::::::::::::::::::#
  # Predict outcomes:::::::#
  #::::::::::::::::::::::::#
  for(i in 1:num.ms){
    y.data.new2[, M[i]] <- PredictM[, i]
  }
  if(isGlm.y){
    y.data$muldm  <- predict(fit.y, newdata = y.data.new2, type = "response")
  } else if (!isGlm.y){
    y.data$muldm  <- predict(fit.y, newdata = y.data.new2)
  }
  if(func == "mmi"){
    sel.grp <- 1
  } else if(func == "smi"){
    sel.grp <- sel.lev.treat
  }
  if(conditional){
    a.f <- as.formula(paste("muldm ~ ", paste(covariates, collapse = "+")))
    a.w <- NULL
  } else if (!conditional){
    a.f <- as.formula(paste("muldm ~ 1"))
    a.w <- y.data[y.data[, treat] == levels(y.data[, treat])[sel.grp], ]$w * prop.treat[sel.grp]
  }
  
  # Compute outcomes after incorporating predicted values of mediator(s) 
  a <- lm(a.f, weights = a.w, data = y.data[y.data[, treat] == levels(y.data[, treat])[sel.grp], ])
  wmuldm <- a$coef[1]
  w.ref.sel <- mean(as.numeric(wmuldm))
  
  #:::::::::::::::#
  #:Results:::::::#
  #:::::::::::::::#
  # Initial disparity, Disparity remaining, and Disparity reduction in order
  out <- rep(NA, 3)
  out[1] <- wy[sel.lev.treat] - wy[ref.lev.treat]
  out[2] <- w.ref.sel - wy[ref.lev.treat]
  out[3] <- wy[sel.lev.treat] - w.ref.sel
  
  return(out)
}