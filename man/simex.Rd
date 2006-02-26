\name{simex}
\alias{simex}
\title{Simulation Extrapolation}
\description{
Implementation of the SIMEX Algorithm for measurement error models according to Cook and Stefanski}
\usage{
simex(model
	, SIMEXvariable
	, measurement.error
	, lambda = c(0.5,1,1.5,2)
	, B = 100, fitting.method = "quadratic"
	, jackknife.estimation = "quad"
	, asymptotic = TRUE)
}
\arguments{
  \item{model}{the naive model}
  \item{SIMEXvariable}{character or vector of characters containing the names of the variables with measurement error}
  \item{measurement.error}{vector of standard deviations of the known measurement errors}
  \item{lambda}{vector of lambdas for which the simulation step should be done (without 0)}
  \item{B}{number of iterations for each lambda}
  \item{fitting.method}{fitting method \code{linear},\code{quadratic},\code{nonlinear} are implemented. (first 4 letters are enough)}
  \item{jackknife.estimation}{specifying the extrapolation method for jackknife variance estimation. Can be set to FALSE if it should not be performed}
  \item{asymptotic}{logical, indicating if asymptotic variance estimation should be done, in the Naive model the option \code{x = TRUE} have to be set.}
}
\details{
nonlinear is implemented as described in Cook and Stefanski, but is numerically not stable. It is not advisable to use this feature. See \code{fit.nls} for details.
If a nonlinear extrapolation is desired please use the \code{refit} function.

Asymptotic is only implemented for naive models of class \code{lm} or \code{glm}
}
\value{Returns an object of class SIMEX which contains:
   \item{coefficients}{the corrected coefficients of the SIMEX model,}
   \item{SIMEX.estimates}{the estimates for every lambda,}
   \item{model}{the naive model,}
   \item{measurement.error}{the known error variances,}
   \item{B}{the number of iterations,}
   \item{extrapolation}{the model object of the extrapolation step,}
   \item{fitting.method}{the fitting method of the extrapolation step,}
   \item{residuals}{residuals,}
   \item{fitted.values}{fitted values,}
   \item{call}{the function call,}
   \item{variance.jackknife}{the jackknife variance estimate,}
   \item{extrapolation.variance}{the model object of the variance extrapolation,}
   \item{variance.jackknife.lambda}{the data set for the extrapolation}
   \item{variance.asymptotic}{the asymptotic variance estimates}
   \item{theta}{estimates for every B and lambda}
   ...
   }
\references{
Cook, J.R. and Stefanski, L.A. (1994) Simulation-Extrapolation Estimation in Parametric Measurement error Models. \emph{Journal of American Statistical Assosiaction}, \bold{89}, 1314 -- 1328

Carroll, R.J., Küchenhoff,H., Lombard,F. and Stefanski L.A. (1996) Asymptotics for the SIMEX Estimator in Nonlinear Measurement Error Models. \emph{Journal of the American Statistical Association}, \bold{91}, 242 -- 250

Carrol, R.J., Ruppert, D. and Stefanski L.A. (1995). \emph{Measurement Error in Nonlinear Models.} London: Chapman and Hall.
}
\author{Wolfgang Lederer,\email{wolfgang.lederer@googlemail.com}}
\seealso{\code{\link{mcsimex}} for discreete data with misclassification, \code{\link[stats]{lm}},\code{\link[stats]{glm}}, \code{\link{refit}}}
\examples{
# to test nonlinear extrapolation set.seed(3)
x <- rnorm(200,0,100)
u <- rnorm(200,0,25)
w <- x+u
y <- x +rnorm(200,0,9)
true.model <- lm(y~x) # True model
naive.model <- lm(y~w, x=TRUE)
simex.model <- simex(model = naive.model
	, SIMEXvariable = "w"
	, measurement.error= 25)
plot(x,y)
abline(true.model,col="darkblue")
abline(simex.model,col ="red")
abline(naive.model,col = "green")
legend(min(x),max(y),legend=c("True Model","SIMEX model","Naive Model")
	, col = c("darkblue","red","green"),lty=1)

plot(simex.model, mfrow = c(2,2))
}
\keyword{models}
