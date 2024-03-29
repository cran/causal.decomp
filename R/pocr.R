#' Product-of-Coefficients-Regression Estimation Method
#'
#' 'pocr' is used to estimate the initial disparity, disparity reduction, and
#' disparity remaining for causal decomposition analysis, using the
#' product-of-coefficients-regression estimation method proposed by Park et al. (2021+).
#' 
#' @usage 
#' pocr(fit.x = NULL, fit.m, fit.y, treat, covariates, sims = 100, conf.level = .95,
#'      cluster = NULL, long = TRUE, mc.cores = 1L, seed = NULL)
#'
#' @details This function returns the point estimates of the initial disparity,
#'   disparity reduction, and disparity remaining for a categorical
#'   treatment and a variety of types of outcome and mediator(s) in causal
#'   decomposition analysis. It also returns nonparametric
#'   percentile bootstrap confidence intervals for each estimate.
#'
#'   The definition of the initial disparity, disparity reduction, and
#'   disparity remaining can be found in help('mmi'). As opposed to the 'mmi' and
#'   'smi' function, this function uses the product-of-coefficients-regression method
#'   suggested by Park et al. (2021+). It always make the inference
#'   conditional on baseline covariates. Therefore, users need to center the data
#'   before fitting the models. See the reference for more details.
#'
#'   As of version 0.1.0, the mediator model ('fit.m') can be of class 'lm', 'glm',
#'   'multinom', or 'polr', corresponding respectively to the linear regression
#'   models and generalized linear models, multinomial log-linear models, and
#'   ordered response models. The outcome model ('fit.y') can be of class 'lm'.
#'   The intermediate confounder model ('fit.x') can also be of class 'lm' and only
#'   necessary when the mediator is categorical.
#'
#' @param fit.x a fitted model object for intermediate confounder. Can be of class 'lm'.
#'   Only necessary if the mediator is categorical. Default is 'NULL'.
#' @param fit.m a fitted model object for mediator. Can be of class 'lm', 'glm',
#'   'multinom', or 'polr'.
#' @param fit.y a fitted model object for outcome. Can be of class 'lm'.
#' @param sims number of Monte Carlo draws for nonparametric bootstrap.
#' @param conf.level level of the returned two-sided confidence intervals,
#'   which are estimated by the nonparametric percentile bootstrap method.
#'   Default is to return the 2.5 and 97.5 percentiles of the simulated quantities.
#' @param treat a character string indicating the name of the treatment variable
#'   used in the models. The treatment can be categorical with two or more
#'   categories (two- or multi-valued factor).
#' @param covariates a vector containing the name of the covariate variable(s)
#'   used in the models. Each covariate can be categorical with two or more
#'   categories (two- or multi-valued factor) or continuous (numeric).
#' @param cluster a vector of cluster indicators for the bootstrap. If provided,
#'   the cluster bootstrap is used. Default is 'NULL'.
#' @param long a logical value. If 'TRUE', the output will contain the entire
#'   sets of estimates for all bootstrap samples. Default is 'TRUE'.
#' @param mc.cores The number of cores to use. Must be exactly 1 on Windows.
#' @param seed seed number for the reproducibility of results. Default is `NULL'.
#'
#' @return
#'
#'   \item{result}{a matrix containing the point estimates of the initial disparity,
#'   disparity remaining, and disparity reduction, and the percentile bootstrap
#'   confidence intervals for each estimate.}
#'   \item{all.result}{a matrix containing the point estimates of the initial disparity,
#'   disparity remaining, and disparity reduction for all bootstrap samples. Returned
#'   if 'long' is 'TRUE'.}
#'
#' @author
#'   Suyeon Kang, University of California, Riverside, \email{skang062@@ucr.edu};
#'   Soojin Park, University of California, Riverside, \email{soojinp@@ucr.edu}.
#'
#' @seealso \code{\link{mmi}}, \code{\link{smi}}
#'
#' @references
#'   Park, S., Kang, S., and Lee, C. (2021+). "Choosing an optimal method for causal
#'   decomposition analysis: A better practice for identifying contributing factors to
#'   health disparities". arXiv preprint arXiv:2109.06940.
#'
#' @export
#' @examples
#' data(sdata)
#'
#' # To be conditional on covariates, first create data with centered covariates
#' # copy data
#' sdata.c <- sdata
#' # center quantitative covariate(s)
#' sdata.c$C.num <- scale(sdata.c$C.num, center = TRUE, scale = FALSE)
#' # center binary (or categorical) covariates(s)
#' # only neccessary if the desired baseline level is NOT the default baseline level.
#' sdata.c$C.bin <- relevel(sdata.c$C.bin, ref = "1")
#'
#' #---------------------------------------------------------------------------------#
#' # Example 1: Continuous Mediator
#' #---------------------------------------------------------------------------------#
#' fit.m1 <- lm(M.num ~ R + C.num + C.bin, data = sdata.c)
#' fit.y1 <- lm(Y.num ~ R + M.num + M.num:R + X +
#'           C.num + C.bin, data = sdata.c)
#' res1 <- pocr(fit.m = fit.m1, fit.y = fit.y1, sims = 40,
#'         covariates = c("C.num", "C.bin"), treat = "R", seed = 111)
#' res1
#'
#' #---------------------------------------------------------------------------------#
#' # Example 2: Binary Mediator
#' #---------------------------------------------------------------------------------#
#' \donttest{fit.x1 <- lm(X ~ R + C.num + C.bin, data = sdata.c)
#' fit.m2 <- glm(M.bin ~ R + C.num + C.bin, data = sdata.c,
#'           family = binomial(link = "logit"))
#' fit.y2 <- lm(Y.num ~ R + M.bin + M.bin:R + X +
#'           C.num + C.bin, data = sdata.c)
#' res2 <- pocr(fit.x = fit.x1, fit.m = fit.m2, fit.y = fit.y2,
#'         sims = 40, covariates = c("C.num", "C.bin"), treat = "R", seed = 111)
#' res2}
#'
#' #---------------------------------------------------------------------------------#
#' # Example 3: Ordinal Mediator
#' #---------------------------------------------------------------------------------#
#' \donttest{require(MASS)
#' fit.m3 <- polr(M.cat ~ R + C.num + C.bin, data = sdata.c)
#' fit.y3 <- lm(Y.num ~ R + M.cat + M.cat:R + X +
#' 	C.num + C.bin, data = sdata.c)
#' res3 <- pocr(fit.x = fit.x1, fit.m = fit.m3, fit.y = fit.y3,
#'        sims = 40, covariates = c("C.num", "C.bin"), treat = "R", seed = 111)
#' res3}
#'
#' #---------------------------------------------------------------------------------#
#' # Example 4: Nominal Mediator
#' #---------------------------------------------------------------------------------#
#' \donttest{require(nnet)
#' fit.m4 <- multinom(M.cat ~ R + C.num + C.bin, data = sdata.c)
#' res4 <- pocr(fit.x = fit.x1, fit.m = fit.m4, fit.y = fit.y3,
#'         sims = 40, covariates = c("C.num", "C.bin"), treat = "R", seed = 111)
#' res4}
pocr <- function(fit.x = NULL, fit.m, fit.y,
                 treat, covariates, sims = 100, conf.level = .95,
                 cluster = NULL, long = TRUE,
                 mc.cores = 1L, seed = NULL){

  if(!is.null(seed)){
    set.seed(seed)
  }
  
  fit.x <- fit.x

  # Warning for inappropriate settings
  if(inherits(fit.m, "list")){
    stop("There must be only one mediator model")
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

  isMultiConfounders <- FALSE

  # Model type indicators
  isLm.x <- inherits(fit.x, "lm") #fit.s[[1]]
  isLm.y <- inherits(fit.y, "lm")

  isGlm.m <- inherits(fit.m, "glm")
  isLm.m <- inherits(fit.m, "lm")
  isNominal.m <- inherits(fit.m, "multinom")
  isOrdinal.m <- inherits(fit.m, "polr")

  # Record family and link of fit.m if glm
  if(isGlm.m){
    FamilyM <- fit.m$family$family
    LinkM <- fit.m$family$link
    if(FamilyM == "gaussian" && LinkM == "identity"){
      M.fun <- gaussian(link = "identity")
    } else if(FamilyM == "binomial" && LinkM == "logit"){
      M.fun <- binomial(link = "logit")
    } else if(FamilyM == "binomial" && LinkM == "probit"){
      M.fun <- binomial(link = "probit")
    } else if(FamilyM == "binomial" && LinkM == "cloglog"){
      M.fun <- binomial(link = "cloglog")
    } else if(FamilyM == "poisson" && LinkM == "log"){
      M.fun <- poisson(link = "log")
    } else if(FamilyM == "poisson" && LinkM == "identity"){
      M.fun <- poisson(link = "identity")
    } else if(FamilyM == "poisson" && LinkM == "sqrt"){
      M.fun <- poisson(link = "sqrt")
    } else if(FamilyM == "Gamma" && LinkM == "inverse"){
      M.fun <- Gamma(link = "inverse")
    } else if(FamilyM == "Gamma" && LinkM == "identity"){
      M.fun <- Gamma(link = "identity")
    } else if(FamilyM == "Gamma" && LinkM == "log"){
      M.fun <- Gamma(link = "log")
    } else if(FamilyM == "inverse.gaussian" && LinkM == "1/mu^2"){
      M.fun <- inverse.gaussian(link = "1/mu^2")
    } else {
      stop("glm family for the mediator model not supported")
    } ### others excluded
  } else if (isOrdinal.m){
    pfun <- switch(fit.m$method, logistic = plogis, probit = pnorm,
                   loglog = pgumbel, cloglog = pGumbel, cauchit = pcauchy)
  }

  if(!is.null(fit.x) && !isLm.x){
    stop("unsupported model")
  }
  #fit.x1 <- fit.x[[1]]
  #fit.x2 <- fit.x[[2]]

  if(!isLm.y){
    stop("unsupported outcome model")
  }

  if(!isGlm.m && !isLm.m && !isNominal.m && !isOrdinal.m){
    stop("unsupported mediator model")
  }

  # Numbers of observations and model frame
  if(!is.null(fit.x) && isLm.x){
    n.x <- nrow(model.frame(fit.x))
  } else {
    n.x <- NULL
  }
  n.m <- nrow(model.frame(fit.m))
  y.data <- model.frame(fit.y)
  n.y <- nrow(y.data)
  if(is.null(n.x)){
    if(n.m != n.y){
      stop("number of observations do not match between mediator
           and outcome models")
    }
  } else if (!is.null(n.x)){
    if(n.m != n.y | n.x != n.y | n.m != n.x){
      stop("number of observations do not match between intermediate confounder,
           mediator and outcome models")
    }
  }

  if(isGlm.m){
    if(FamilyM == "binomial" & is.null(fit.x)){
      stop("'fit.x' must be specified for a binary mediator")
    } else if(FamilyM != "binomial" & !is.null(fit.x)){
      fit.x <- NULL
      message("'fit.x' must be NULL for a quantitative mediator. fit.x = NULL forced")
    }
  } else if (isLm.m & !isGlm.m){
    if(!is.null(fit.x)){
      fit.x <- NULL
      message("'fit.x' must be NULL for a quantitative mediator. fit.x = NULL forced")
    }
  } else if (isNominal.m & is.null(fit.x)){
    stop("'fit.x' must be specified for a categorical mediator")
  } else if (isOrdinal.m & is.null(fit.x)){
    stop("'fit.x' must be specified for a categorical mediator")
  }
  
  # Model frames for M and Y models (ver 0.1.0)
  m.data <- model.frame(fit.m)
  y.data <- model.frame(fit.y)
  
  if(!is.null(cluster)){
    row.names(m.data) <- 1:nrow(m.data)
    row.names(y.data) <- 1:nrow(y.data)
    
    if(!is.null(fit.m$weights)){
      m.weights <- as.data.frame(fit.m$weights)
      m.name <- as.character(fit.m$call$weights)  
      names(m.weights) <- m.name
      m.data <- cbind(m.data, m.weights)
    }
    
    if(!is.null(fit.y$weights)){
      y.weights <- as.data.frame(fit.y$weights)
      y.name <- as.character(fit.y$call$weights)  
      names(y.weights) <- y.name
      y.data <- cbind(y.data, y.weights)
    }
  }
  
  # Extracting weights from models (ver 0.1.0)
  weights.m <- model.weights(m.data)
  weights.y <- model.weights(y.data)
  
  if(!is.null(weights.m) && isGlm.m && FamilyM == "binomial"){
    message("weights taken as sampling weights, not total number of trials")
  }
  if(is.null(weights.m)){
    weights.m <- rep(1, nrow(m.data))
  }
  if(is.null(weights.y)){
    weights.y <- rep(1, nrow(y.data))
  }
  weights <- weights.y

  # Extract treatment and outcome variable names from models
  if(!is.null(fit.x)){
    X <- all.vars(formula(fit.x))[[1]]
  }
  M <- all.vars(formula(fit.m))[[1]]
  r <- length(levels(y.data[, treat]))
  outcome <- all.vars(formula(fit.y))[[1]]

  names <- c(treat, covariates)
  for(nn in 1:length(names)){
    if(!(names[nn] %in% attr(terms(fit.m),"term.labels"))){
      stop(paste(names[nn] , " must be in 'fit.m'", sep = ""))
    }
    if (!is.null(fit.x) && !(names[nn] %in% attr(terms(fit.x),"term.labels"))){
      stop(paste(names[nn] , " must be in 'fit.x'", sep = ""))
    }
  }
  if(!is.null(fit.x)){
    names <- c(treat, covariates, M, X)
  } else {
    names <- c(treat, covariates, M)
  }
  for(nn in 1:length(names)){
    if(!(names[nn] %in% attr(terms(fit.y),"term.labels"))){
      stop(paste(names[nn] , " must be in 'fit.y'", sep = ""))
    }
    if(!(paste(treat, M, sep = ":") %in% attr(terms(fit.y),"term.labels")) &&
       !(paste(M, treat, sep = ":") %in% attr(terms(fit.y),"term.labels"))){
      stop("interaction term between treatment and mediator must be in 'fit.y'")
    }
  }

  # Bootstrap QoI
  message("Running nonparametric bootstrap\n")

  # Make objects available to all functions
  environment(combine.b) <- environment(combine) <- environment(select.treat.pocr) <- environment()

  out <- matrix(NA, 3 * (r - 1), 3)
  all.out <- matrix(NA, 3 * (r - 1), sims)
  colnames(out) <- c("estimate", paste(conf.level * 100, "% CI Lower", sep = ""),
                     paste(conf.level * 100, "% CI Upper", sep = ""))
  rn <- rn2 <- NULL
  for(rr in 2:r){
    rn <- c(rn, c(paste("Initial Disparity   (", levels(y.data[, treat])[1],
                        " vs ", levels(y.data[, treat])[rr], ")", sep = ""),
                  paste("Disparity Remaining (", levels(y.data[, treat])[1],
                        " vs ", levels(y.data[, treat])[rr], ")", sep = ""),
                  paste("Disparity Reduction (", levels(y.data[, treat])[1],
                        " vs ", levels(y.data[, treat])[rr], ")", sep = "")))
    rn2 <- c(rn2, levels(y.data[, treat])[rr])
  }

  summary0 <- mclapply(1:sims, combine.b, fit.r = NULL, fit.x = fit.x, fit.m = fit.m, fit.y = fit.y,
                       cluster = cluster, func = "pocr", mc.cores = mc.cores)
  out.full <- simplify2array(summary0)
  ind.alphar <- seq(4, 5 * (r - 1), 5)
  ind.segamma <- seq(5, 5 * (r - 1), 5)
  ind.exclude <- c(ind.alphar, ind.segamma)
  
  out[, 1] <- combine(fit.r = NULL, fit.x = fit.x, fit.m = fit.m, fit.y = fit.y, func = "pocr", weights = weights)[- ind.exclude]
  out[, 2] <- apply(out.full[- ind.exclude, ], 1, quantile, prob = (1 - conf.level)/2, na.rm = TRUE)
  out[, 3] <- apply(out.full[- ind.exclude, ], 1, quantile, prob = 1/2 + conf.level/2, na.rm = TRUE)
  all.out <- out.full[- ind.exclude, ]
  rownames(out) <- rownames(all.out) <- rn
  
  alpha.r <- t(matrix(out.full[ind.alphar, ], ncol = sims))
  se.gamma <- t(matrix(out.full[ind.segamma, ], ncol = sims))
  colnames(alpha.r) <- colnames(se.gamma) <- rn2
  
  if(long){
    out <- list(result = out, all.result = all.out)
  } else {
    out <- list(result = out)
  }
  class(out) <- "pocr"
  return(out)

}

pgumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE){

  q <- (q - loc)/scale
  p <- exp(- exp(- q))
  if(!lower.tail){
    return(1 - p)
  } else {
    return(p)
  }
}

pGumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE){

  q <- (q - loc)/scale
  p <- exp(- exp(q))
  if(lower.tail){
    return(1 - p)
  } else {
    return(p)
  }
}
