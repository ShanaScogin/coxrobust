# SRS note 3/31/2022: add in the following:
# Covariates are generated independently, each from the standard normal
# distribution.  Baseline hazard is equal to 1.  After generation of survival
# times, covariates are replaced by independent \eqn{2N(0,1)+1} variables
# in fraction \code{cont} of observations.

#' Generate Data from the Proportional Hazards Regression Model
#'
#' Generates data set from the proportional hazards regression model
#' without or with contamination.
#'
#' @param n number of observations.
#' @param beta vector of regression coefficients.
#' @param cont fraction of contaminated observations.
#' @param p.censor probability of censoring.
#'
#' @return Data frame containing the following variables:
#' \itemize{
#' \item{time}{vector of survival times.}
#' \item{status}{vector of censoring status.}
#' \item{X1, X2, ...}{explanatory variables (their number is determined by the
#' dimension of vector of regression coefficients).}
#' }
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' \donttest{
#' if (interactive()) {
#' gen_data(50, c(2,-2), cont = 0.05)
#' }
#' }
#'
#' \dontshow{setwd(.old_wd)}
#' @export

gen_data <- function(n, beta, cont = 0, p.censor = 0) {

	beta <- as.double(beta)
	m <- length(beta)

    z <- array(rnorm(n*m), c(n,m))

    time <- rexp(n) / exp( z %*% beta )

	status <- sample(c(0,1), n, replace = TRUE, prob = c(p.censor, 1-p.censor))

	ncont <- floor(m*n*cont)
	if ( ncont > 0 ) {
		z <- as.double(t(z))
    	z[1:ncont] <- 2*rnorm(ncont) + 1
		z <- matrix(z, n, m, byrow = TRUE)
	}

    if ( ncol(z) == 1 ) {
        gdata <- data.frame(time, status, X1 = z)
    } else {
        gdata <- data.frame(time, status, z)
    }

    return(gdata)

}
