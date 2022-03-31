\name{plot.coxr}
\alias{plot.coxr}

\title{Plot Diagnostics for a coxr Object}

\description{
Graphical tool which in a series of 5 graphs let us compare how well data are
explained by the estimated proportional hazards model with non-robust
(black color) and robust method (green color).  The first graph gives
standardized difference of two estimated survival functions; one via the Cox
model and the other via Kaplan Meier estimator.  The following four graphs
show the same differences for four strata, defined by the quartiles of the
estimated linear predictor.  Comparison of estimation results along with
analysis of the graphs leads frequently to a very detailed information about
the model fit (see examples).
}

\usage{
\method{plot}{coxr}(x, caption = c("Full data set", "First quartile",
    "Second quartile", "Third quartile", "Fourth quartile"), main = NULL,
    xlab = "log time", ylab = "standardized survival differences", \dots,
    color = TRUE)
}

\arguments{

\item{x}{\code{coxr} object, typically result of \code{\link{coxr}}.}
\item{caption}{captions to appear above the plots.}
\item{xlab}{title for the x axis.}
\item{ylab}{title for the y axis.}
\item{main}{overall title for the plot.}
\item{\dots}{other parameters to be passed through to plotting functions.}
\item{color}{if \code{FALSE} grayscale mode is used.}

}

\seealso{
\code{\link{coxr}}
}

\examples{
#use the lung cancer data at Mayo Clinic 
#to compare results of non-robust and robust estimation
result <- coxr(Surv(time, status) ~ age + sex + ph.karno + meal.cal + wt.loss, data = lung)
plot(result, main = "Mayo Clinic Lung Cancer Data")
}

\keyword{hplot}
\keyword{robust}
\keyword{survival}
