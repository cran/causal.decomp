# Bootstrappting with 'combine' function
combine.b <- function(XX, fit.r, fit.x, fit.m, fit.y, cluster = NULL, func){

  if(is.null(cluster)){
    # sample the observations with replacement
    sb <- sample(1:nrow(y.data), n.y, replace = TRUE)
  } else {
    clusters <- unique(cluster)
    # sample the clusters with replacement
    units <- sample(clusters, size = length(clusters), replace = TRUE)
    # create bootstrap sample with sapply
    df.bs <- sapply(units, function(x) which(cluster == x))
    sb <- unlist(df.bs)
  }
  # create bootstrap sample
  y.data <- model.frame(fit.y)
  y.data.b <- y.data[sb, ]
  if(!"(weights)" %in% colnames(y.data.b)){  # version 0.1.0
    y.data.b <- cbind(y.data.b, weights[sb])
    colnames(y.data.b)[ncol(y.data.b)] <- "(weights)"
  }

  #assign("y.data.b", y.data.b, envir = .GlobalEnv)
  #assign("weights.b", weights.b, envir = .GlobalEnv)
  fit.y.b <- update(fit.y, data = y.data.b, weights = y.data.b[, "(weights)"])
  
  if(!is.null(fit.r) && isCBPS.r){
    fit.r.b <- update(fit.r, data = y.data.b, sample.weights = y.data.b[, "(weights)"])##
  } else if(!is.null(fit.r) && isSumStat.r){
    ps.xvar <- rownames(fit.r$unweighted.sumstat)
    for(name in 1:length(ps.xvar)){
      if(!ps.xvar[name] %in% colnames(y.data.b)){
        ps.xvar[name] <- substr(ps.xvar[name], 1, nchar(ps.xvar[name]) - 1)
      }
    }
    ps.formula <- as.formula(paste(treat, paste(ps.xvar, collapse = " + "), sep = " ~ "))
    fit.r.b <- SumStat(ps.formula = ps.formula, data = y.data.b,
                       trtgrp = fit.r$trtgrp, weight = "IPW")
  } else {
    fit.r.b <- NULL
  }
  
  if(!is.null(fit.x)){
    fit.x.b <- update(fit.x, data = y.data.b, weights = y.data.b[, "(weights)"])
  #  fit.x1.b <- update(fit.x1, data = y.data.b)
  #  fit.x2.b <- update(fit.x2, data = y.data.b)
  } else {
    fit.x.b <- NULL
  }
  
  if(!isMultiConfounders){
    if(isNominal.m){
      fit.m.b <- update(fit.m, data = y.data.b, trace = FALSE, weights = y.data.b[, "(weights)"])
    } else {
      fit.m.b <- update(fit.m, data = y.data.b)
    }
  } else if(isMultiConfounders){
    fit.m.b <- list()
    for(i in 1:length(fit.m)){
      if(isNominal.m[i]){
        fit.m.b[[i]] <- update(fit.m[[i]], data = y.data.b, trace = FALSE, weights = y.data.b[, "(weights)"])
      } else {
        fit.m.b[[i]] <- update(fit.m[[i]], data = y.data.b)#, weights = NULL)
      }
    }
  }
  
  environment(combine) <- environment()
  out <- combine(fit.r = fit.r.b, fit.x = fit.x.b, fit.m = fit.m.b, fit.y = fit.y.b, func = func,
                 weights = y.data.b[, "(weights)"])
  return(out)
}
