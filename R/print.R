### objects from sensitivity.R ========================================================
##' @export
print.sensitivity <- function(x, ...){
  
  cat("\nDisparity Reduction:\n\n")
  redcolnames <- colnames(x$disparity.remaining)
  redout <- data.frame(x$disparity.remaining)
  colnames(redout) <- redcolnames
  print(redout)
  cat("\nDisparity Remaining:\n\n")
  remcolnames <- colnames(x$disparity.reduction)
  remout <- data.frame(x$disparity.reduction)
  colnames(remout) <- remcolnames
  print(remout)
  cat("\n")
  
  invisible(x)
  
}

### for objects from smi.R =====================================================
##' @export
print.smi <- function(x, ...){
  
  cat("\nResults:\n\n")
  print(x$result)
  cat("\n")
  
  invisible(x)
  
}

### for objects from mmi.R =====================================================
##' @export
print.mmi <- function(x, ...){
  
  cat("\nResults:\n\n")
  print(x$result)
  cat("\n")
  
  invisible(x)
  
}

### for objects from pocr.R =====================================================
##' @export
print.pocr <- function(x, ...){
  
  cat("\nResults:\n\n")
  print(x$result)
  cat("\n")
  
  invisible(x)
  
}