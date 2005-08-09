"print.summary.SIMEX" <-
function(x, digits = max(3, getOption("digits") - 3),...){
cat("Call:\n")
print(x$call)
  cat("\nNaive model: \n")
  print(x$naive.model)
  cat("\nSimex variable :\n")
  var.error <- t(matrix(x$measurement.error))
  dimnames(var.error) <- list("Measurement error :", x$SIMEXvariable )
  print(var.error)
  cat("\n\n\Number of iterations: ",x$B,"\n")
cat("\nResiduals: \n")
print(summary(x$residuals),digits)
cat("\n\Coefficients: \n")
  if(any(names(x$coefficients)=="asymptotic"))
  {
   cat("\n\Asymptotic variance: \n")
   printCoefmat(x$coefficients$asymptotic, digits =digits)
  }
  if(any(names(x$coefficients)=="jackknife"))
  {
   cat("\n\Jackknife variance: \n")
   printCoefmat(x$coefficients$jackknife, digits =digits)
  }
}

