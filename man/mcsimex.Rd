\name{mcsimex}
\alias{mcsimex}
\title{ The Misclassification SIMEX }
\description{
  Implementation of the Misclassification SIMEX Algorithm as described by Küchenhoff, Mwalili and Lesaffre.
}
\usage{
mcsimex(model
	, SIMEXvariable
	, mc.matrix
	, lambda = c(0.5,1,1.5,2)
	, B = 100
	, jackknife.estimation = "quad"
	, asymptotic = TRUE
	, fitting.method = "quad")
}
\arguments{
  \item{model}{ The naive model, the misclassified variable must be a factor}
  \item{mc.matrix}{If one variable is misclassified it can be a matrix. If more than one variable is misclssified it must be a list of the misclassification matrices, names must match with the SIMEXvariabel names, column- and row-names must match with the factor levels. If a special misclssification is desired, the name of a function can be specified (see details)}
  \item{lambda}{ vector of exponents for the misclassification matrix (withou 0)}
  \item{SIMEXvariable}{ vector of names of the variables for which the MCSIMEX-method should be applied}
  \item{B}{ number of iterations for each lambda}
  \item{fitting.method}{ linear, quadratic and loglinear are implemented (first 4 letters are enough)}
  \item{jackknife.estimation}{ specifying the extrapolation method for jackknife variance estimation. Can be set to FALSE if it should not be performed}
  \item{asymptotic}{ logical, indicating if asymptotic variance estimation should be done, the option \code{x =TRUE} must be enabled in the naive model.}
  }
\details{
	if \code{mc.matrix} is a function the first argument of that function must 
	be the whole dataset used in the naive model, the second argument must be 
	the exponent (lambda) for the the misclassification. The function must return
	a \code{data.frame} containing the misclassified \code{SIMEXvariable}. An example can be found below.
	
	Asymptotic variance estimation is only implemented for \code{lm} and \code{glm}
	
	the loglinear fit has the form g(lambda,GAMMA) = exp(gamma0+gamma1*lambda). It is realized via the log() function. To avoid negaitve
	values the minimum +1 of the dataset is added and after the prediction later subtracted. exp(predict(...)) - min(data)-1
}
\value{
  object of class \code{MCSIMEX}
  \item{coefficients}{ corrected coefficients of the MCSIMEX-model }
  \item{SIMEX.estimates}{the MCSIMEX-estimates of the coefficients for each lambda}
  \item{lambda}{ the values of lambda}
  \item{model}{ naive model}
  \item{mc.matrix}{ the misclassification matrix}
  \item{B}{ the number of iterations}
  \item{extrapolation}{ model-object of the extrapolation step }
  \item{fitting.method}{ the fitting method used in the exrapolation step}
  \item{SIMEXvariable}{ name of the SIMEXvariables}
  \item{call}{ the function call,}
  \item{variance.asymptotic}{ the asymptotic variance estimates}
  \item{variance.jackknife}{ the jackknife variance estimates}
  \item{extrapolation.variance}{ the model-object of the variance extrapolation}
  \item{variance.jackknife.lambda}{ data set for the extrapolation}
  \item{theta}{all estimated coefficients for each lambda and B}
  ...
  }
\references{ Küchenhoff, H., Mwalili, S. M.  and Lesaffre (2005) E. A general method for dealing with misclassification in regression: the Misclassification SIMEX. \emph{Biometrics},in press}
\author{ Wolfgang Lederer, \email{wolfgang.lederer@stat.uni-muenchen.de}}
\seealso{ \code{\link{misclass}}, \code{\link{simex}},\code{\link{refit}}}
\examples{
x <- rnorm(200,0,1.142)
z <- rnorm(200,0,2)
y <- factor(rbinom(200,1,(1/(1+exp(-1*(-2 + 1.5*x -0.5*z))))))
Pi <- matrix(data = c(0.9,0.1,0.3,0.7), nrow =2, byrow =FALSE)
dimnames(Pi) <- list(levels(y),levels(y))
ystar <- misclass(data.frame(y), list(y = Pi), k=1)[,1]
naive.model <- glm(ystar ~ x + z, family = binomial, x=TRUE, y =TRUE)
true.model  <- glm(y ~ x + z, family = binomial)
simex.model <- mcsimex(naive.model, mc.matrix = Pi, SIMEXvariable = "ystar")

op <-par(mfrow = c(2,3))
invisible(lapply(simex.model$theta, boxplot, notch=TRUE, outline =FALSE, names=c(0.5,1,1.5,2)))
plot(simex.model)
par(op)


## example for a function which can be supplied to the function mcsimex()
## "xm" is the variable which is to be misclassified

my.mc <- function(datas,k){
	xm <- datas$"xm"
	p1 <- matrix(data = c(0.75,0.25,0.25,0.75), nrow =2, byrow = FALSE)
	colnames(p1) <- levels(xm)
	rownames(p1) <- levels(xm)
	p0 <- matrix(data = c(0.8,0.2,0.2,0.8), nrow =2, byrow =FALSE)
	colnames(p0) <- levels(xm)
	rownames(p0) <- levels(xm)
	xm[datas$y=="1"] <- misclass(data.frame(xm=xm[datas$y=="1"]),list(xm=p1), k=k)[,1]
	xm[datas$y=="0"] <- misclass(data.frame(xm=xm[datas$y=="0"]),list(xm=p0), k=k)[,1]
	xm <- factor(xm)
	return(data.frame(xm))
	}

}
\keyword{models}
