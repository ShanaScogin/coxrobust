\name{coxr}
\alias{coxr}
\alias{print.coxr}
\alias{predict.coxr}

\title{Fit Robustly Proportional Hazards Regression Model}

\description{
Fits efficiently and robustly Cox proportional hazards regression model in its
basic form, where explanatory variables are time independent with one event
per subject.  Method is based on a smooth modification of the partial
likelihood.
}

\usage{

coxr(formula, data, subset, na.action, trunc = 0.95,
     f.weight = c("linear", "quadratic", "exponential"),
     singular.ok = TRUE, model = FALSE)

}

\arguments{

\item{formula}{
a formula object, with the response on the left of a \code{~}
operator, and the terms on the right.  The response must be a
survival object as returned by the \code{\link{Surv}} function.}

\item{data}{a data frame in which to interpret the variables
named in the \code{formula}, or in the \code{subset}.}

\item{subset}{expression saying that only a subset of the rows of the data
should be used in the fit.}

\item{na.action}{a missing-data filter function, applied to the model.frame,
after any subset argument has been used.}

\item{trunc}{roughtly, quantile of the sample \eqn{T_i exp(\beta'Z_i)},
it determines the trimming level for the robust estimator}

\item{f.weight}{type of weighting function, default is \code{"quadratic"}}

\item{singular.ok}{logical value indicating how to handle collinearity in the
model matrix. If \code{TRUE}, the program will automatically skip
over columns of the X matrix that are linear combinations of
earlier columns.  In this case the coefficients for such
columns will be \code{NA}, and the variance matrix will contain
zeros.  For ancillary calculations, such as the linear
predictor, the missing coefficients are treated as zeros.}

\item{model}{a logical value indicating whether model frame should be
included as a component of the returned value.}

}

\details{
The method consists in maximization of an objective function which is
a smooth modification of the partial likelihood. Observations with excessive
values of \eqn{\Lambda(T) exp(\beta'Z)}, where \eqn{\Lambda} is the cumulated
hazard, \eqn{\beta} vector of parameters, \eqn{Z} explanatory variables and
\eqn{T} possibly censored survival time, are down-weighted.  Both \eqn{\Lambda}
and \eqn{\beta} are iteratively robustly estimated.

Numerical results are supported by a graphical tool \code{plot}, which in
a series of 5 graphs let us compare how well data are explained by the estimated
proportional hazards model with non-robust (black color) and robust method
(green color).  The first graph shows standardized difference of two estimated
survival functions; one via the Cox model and the other via Kaplan Meier
estimator.  The following four graphs show the same differences for four
strata, defined by the quartiles of the estimated linear predictor.
Comparison of estimation results along with analysis of the graphs leads
frequently to a very detailed information about the model fit (see examples).
}

\value{
An object of class \code{coxr}. See \code{\link{coxr.object}} for details.
}

\references{

Bednarski, T. (1993). Robust estimation in Cox's regression model.
Scandinavian Journal of Statistics. Vol. 20, 213--225.

Bednarski, T. (1989). On sensitivity of Cox's estimator.
Statistics and Decisions. 7, 215--228.

Grzegorek, K.(1993). On robust estimation of baseline hazard under the Cox model and via
Frechet differentiability. Preprint of the Institute of Mathematics of the Polish Academy
of Sciences.518.

Minder, C.E. & Bednarski, T. (1996). A robust method for proportional hazards
regression. Statistics in Medicine Vol. 15, 1033--1047.

}

\seealso{
\code{\link{coxph}}, \code{\link{coxr.object}}, \code{\link{gen_data}},
\code{\link{plot.coxr}}, \code{\link{Surv}}
}

\examples{
# Create a simple test data set using the attached function gen_data
a <- gen_data(200, c(1, 0.1, 2), cont = 0.05, p.censor = 0.30)
result <- coxr(Surv(time, status) ~ X1 + X2 + X3, data = a , trunc = 0.9)
result
plot(result)

#use the lung cancer data at Mayo Clinic to
#compare results of non-robust and robust estimation
result <- coxr(Surv(time, status) ~ age + sex + ph.karno + meal.cal + wt.loss, data = lung)
result
plot(result)

#use the Veteran's Administration Lung Cancer Data
#to compare results of non-robust and robust estimation
result <- coxr(Surv(time,status) ~ age + trt + celltype + karno + diagtime + prior, data = veteran)
result
plot(result)
}

\keyword{robust}
\keyword{survival}
