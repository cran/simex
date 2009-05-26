\name{simex}
\Rdversion{1.1}
\encoding{latin1}
\alias{simex}
\alias{print.simex}
\alias{summary.simex}
\alias{print.summary.simex}
\alias{plot.simex}
\alias{predict.simex}
\alias{refit}
\alias{refit.simex}

\title{Measurement error in models using SIMEX}

\description{
  Implementation of the SIMEX algorithm for measurement error models according
  to Cook and Stefanski}

\usage{
simex(model, SIMEXvariable, measurement.error, lambda = c(0.5, 1, 1.5, 2),
  B = 100, fitting.method = "quadratic", jackknife.estimation = "quadratic",
  asymptotic = TRUE)

\method{print}{simex}(x, digits = max(3, getOption("digits") - 3), ...)
\method{summary}{simex}(object,...)
\method{plot}{simex}(x, xlab = expression((1 + lambda)), ylab = colnames(b[, -1]),
  ask = FALSE, show = rep(TRUE, NCOL(b) - 1), ...)
\method{predict}{simex}(object, newdata, ...)

refit(object, ...)
\method{refit}{simex}(object, fitting.method = "quadratic",
  jackknife.estimation = "quadratic", asymptotic = TRUE, ...)
}

\arguments{
  \item{model}{the naive model}
  \item{SIMEXvariable}{character or vector of characters containing the names of
    the variables with measurement error}
  \item{measurement.error}{vector of standard deviations of the known measurement
    errors}
  \item{lambda}{vector of lambdas for which the simulation step should be done
    (without 0)}
  \item{B}{number of iterations for each lambda}
  \item{fitting.method}{fitting method \code{linear}, \code{quadratic},
    \code{nonlinear} are implemented. (first 4 letters are enough)}
  \item{jackknife.estimation}{specifying the extrapolation method for jackknife
    variance estimation. Can be set to \code{FALSE} if it should not be performed}
  \item{asymptotic}{logical, indicating if asymptotic variance estimation should
    be done, in the naive model the option \code{x = TRUE} has to be set}
  \item{x}{object of class 'simex'}
  \item{digits}{number of digits to be printed}
  \item{object}{object of class 'simex'}
  \item{xlab}{optional name for the X-Axis}
  \item{ylab}{vector containing the names for the Y-Axis}
  \item{ask}{logical. If \code{TRUE}, the user is asked for input, before a new
    figure is drawn}
  \item{show}{vector of logicals indicating for wich variables a plot should be
    produced}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the fitted linear predictors are used}
  \item{\dots}{arguments passed to other functions}
}

\details{
  Nonlinear is implemented as described in Cook and Stefanski, but is numerically
  instable. It is not advisable to use this feature. If a nonlinear extrapolation
  is desired please use the \code{refit()} method.

  Asymptotic is only implemented for naive models of class \code{lm} or \code{glm}.

  \code{refit()} refits the object with a different extrapolation function.
}

\value{
  An object of class 'simex' which contains:
   \item{coefficients}{the corrected coefficients of the SIMEX model,}
   \item{SIMEX.estimates}{the estimates for every lambda,}
   \item{model}{the naive model,}
   \item{measurement.error}{the known error standard deviations,}
   \item{B}{the number of iterations,}
   \item{extrapolation}{the model object of the extrapolation step,}
   \item{fitting.method}{the fitting method used in the extrapolation step,}
   \item{residuals}{the residuals,}
   \item{fitted.values}{the fitted values,}
   \item{call}{the function call,}
   \item{variance.jackknife}{the jackknife variance estimate,}
   \item{extrapolation.variance}{the model object of the variance extrapolation,}
   \item{variance.jackknife.lambda}{the data set for the extrapolation,}
   \item{variance.asymptotic}{the asymptotic variance estimates,}
   \item{theta}{the estimates for every B and lambda,}
   ...
   }

\references{
Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in
  parametric measurement error models. \emph{Journal of American Statistical
  Association}, \bold{89}, 1314 -- 1328

Carroll, R.J., Küchenhoff, H., Lombard, F. and Stefanski L.A. (1996) Asymptotics
  for the SIMEX estimator in nonlinear measurement error models. \emph{Journal
  of the American Statistical Association}, \bold{91}, 242 -- 250

Carrol, R.J., Ruppert, D., Stefanski, L.A. and Crainiceanu, C. (2006).
  \emph{Measurement error in nonlinear models: A modern perspective.}, Second
  Edition. London: Chapman and Hall.
 
Lederer, W. and Küchenhoff, H. (2006) A short introduction to the SIMEX and MCSIMEX. \emph{R News}, \bold{6(4)}, 26--31 
 
}

\author{Wolfgang Lederer,\email{wolfgang.lederer@googlemail.com}}

\seealso{\code{\link[simex]{mcsimex}} for discrete data with misclassification,
  \code{\link[stats]{lm}}, \code{\link[stats]{glm}}}

\examples{
# to test nonlinear extrapolation
#set.seed(3)
x <- rnorm(200, 0, 100)
u <- rnorm(200, 0, 25)
w <- x + u
y <- x + rnorm(200, 0, 9)
true.model <- lm(y ~ x)
naive.model <- lm(y ~ w, x = TRUE)
simex.model <- simex(model = naive.model, SIMEXvariable = "w",
    measurement.error = 25)
plot(x, y)
abline(true.model, col = "darkblue")
abline(simex.model, col = "red")
abline(naive.model, col = "green")
legend(min(x), max(y), legend = c("True Model", "SIMEX model", "Naive Model"),
    col = c("darkblue", "red", "green"), lty = 1)

plot(simex.model, mfrow = c(2, 2))

simex.model2 <-  refit(simex.model, "line")
plot(simex.model2, mfrow = c(2, 2))
}

\keyword{models}
