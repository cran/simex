"print.summary.MCSIMEX" <-
function(x, digits = max(3, getOption("digits") - 3), ...){
cat("Call:\n")
print(x$call)
  cat("\nNaive model: \n")
  print(x$naive.model)
  cat("\nSimex variable :", x$SIMEXvariable, "\n")
  cat("Misclassification matrix: \n")
  if(is.character(x$mc.matrix)) dput(eval(my.mc)) else lapply(x$mc.matrix, print)
  cat("\n\Number of iterations: ",x$B,"\n")
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

