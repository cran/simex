"print.SIMEX" <-
function (
x# object to be printed
, digits = max(3, getOption("digits") - 3)# number of digits to be printed
, ...) 
{
     cat("\nNaive model:\n", deparse(x$model$call) , "\n", sep = "")
    cat("\nSIMEX-Variables: ")
    cat(x$SIMEXvariable, sep =", ")
    cat("\nNumber of Simulations: ", paste(x$B), "\n\n" , sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2,
            quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
}

