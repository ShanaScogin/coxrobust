\name{coxr.object}
\alias{coxr.object}

\title{Fit Robustly Proportional Hazards Regression Object}

\description{
This class of objects is returned by \code{coxr} function to represent
efficiently and robustly fitted proportional hazards regression model.
Objects of this class have methods for the functions \code{print},
\code{plot} and \code{predict}.
}

\value{

The following components must be included in a legitimate \code{coxr} object.

\item{coefficients}{robust estimate of the regression parameter.}

\item{ple.coefficients}{non-robust (efficient) estimate of the regression
parameter.}

\item{var}{
an approximate variance matrix of the coefficients (estimated robustly). 
Rows and columns corresponding to any missing coefficients are set to zero.}

\item{ple.var}{
an approximate variance matrix of the coefficients (estimated non-robustly).
Rows and columns corresponding to any missing coefficients are set to zero.}

\item{lambda}{cumulated hazard (estimated robustly).}

\item{lambda.ple}{cumulated hazard (estimated non-robustly).}

\item{wald.test}{the value of Wald test.}

\item{ewald.test}{the value of extended Wald test.}

\item{skip}{skipped columns.}

\item{na.action}{the \code{na.action} attribute, if any, that was returned
by the \code{na.action} routine.}

The object also contain the following, for documentation see the
\code{\link{lm}} object: \code{terms}, \code{call}, \code{x}, \code{y}
and optionally \code{model}.

}

\seealso{
\code{\link{coxr}}
}

\keyword{robust}
\keyword{survival}
