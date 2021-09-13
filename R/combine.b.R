# Bootstrappting with 'combine' function
combine.b <- function(X, fit.r, fit.s, fit.m, fit.y, cluster = NULL, func){
  
  if(is.null(cluster)){
    # sample the observations with replacement
    sb <- sample(1:nrow(y.data), n.y, replace = TRUE)
    # create bootstap sample
    y.data.b <- y.data[sb, ]
  } else {
    clusters <- unique(cluster)
    # sample the clusters with replacement
    units <- sample(clusters, size = length(clusters), replace = TRUE)
    # create bootstap sample with sapply
    df.bs <- sapply(units, function(x) which(cluster == x))
    y.data.b <- y.data[unlist(df.bs), ]
  }
  
  fit.y.b <- update(fit.y, data = y.data.b)
  
  if(!is.null(fit.r) && isCBPS.r){
    fit.r.b <- update(fit.r, data = y.data.b)
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
  
  if(!is.null(fit.s)){
    fit.s.b <- update(fit.s, data = y.data.b)
  } else {
    fit.s.b <- NULL
  }
  
  if(!isMultiConfounders){
    if(isNominal.m){
      fit.m.b <- update(fit.m, data = y.data.b, trace = FALSE)
    } else {
      fit.m.b <- update(fit.m, data = y.data.b)
    }
  } else if(isMultiConfounders){
    fit.m.b <- list()
    for(i in 1:length(fit.m)){
      if(isNominal.m[i]){
        fit.m.b[[i]] <- update(fit.m[[i]], data = y.data.b, trace = FALSE)
      } else {
        fit.m.b[[i]] <- update(fit.m[[i]], data = y.data.b)
      }
    }
  }
  
  out <- combine(fit.r = fit.r.b, fit.s = fit.s.b, fit.m = fit.m.b, fit.y = fit.y.b, func = func)
  return(out)
}