\name{gen_data}
\alias{gen_data}

\title{Generate Data from the Proportional Hazards Regression Model}

\description{
Generates data set from the proportional hazards regression model
without or with contamination.
}

\usage{gen_data(n, beta, cont = 0, p.censor = 0)}

\arguments{

\item{n}{number of observations.}
\item{beta}{vector of regression coefficients.}
\item{cont}{fraction of contaminated observations.}
\item{p.censor}{probability of censoring.}

}

\details{
Covariates are generated independently, each from the standard normal
distribution.  Baseline hazard is equal to 1.  After generation of survival
times, covariates are replaced by independent \eqn{2N(0,1)+1} variables
in fraction \code{cont} of observations.
}

\value{

Data frame containing the following variables:

\item{time}{vector of survival times.}
\item{status}{vector of censoring status.}
\item{X1, X2, \ldots}{explanatory variables (their number is determined by the
dimension of vector of regression coefficients).}

}

\seealso{
\code{\link{coxph}}, \code{\link{coxr}}
}

\examples{
gen_data(50, c(2,-2), cont = 0.05)
}

\keyword{datagen}
\keyword{survival}
\keyword{robust}
