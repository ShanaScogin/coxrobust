# SRS note 3/31/2022: check the example in this and add or call
# data if need to. Also add in the following from original man file:
#
# The method consists in maximization of an objective function which is
# a smooth modification of the partial likelihood. Observations with excessive
# values of \eqn{\Lambda(T) exp(\beta'Z)}, where \eqn{\Lambda} is the cumulated
# hazard, \eqn{\beta} vector of parameters, \eqn{Z} explanatory variables and
# \eqn{T} possibly censored survival time, are down-weighted.  Both \eqn{\Lambda}
# and \eqn{\beta} are iteratively robustly estimated.
#
# Numerical results are supported by a graphical tool \code{plot}, which in
# a series of 5 graphs let us compare how well data are explained by the estimated
# proportional hazards model with non-robust (black color) and robust method
# (green color).  The first graph shows standardized difference of two estimated
# survival functions; one via the Cox model and the other via Kaplan Meier
# estimator.  The following four graphs show the same differences for four
# strata, defined by the quartiles of the estimated linear predictor.
# Comparison of estimation results along with analysis of the graphs leads
# frequently to a very detailed information about the model fit (see examples).


#' Fit Robustly Proportional Hazards Regression Model
#'
#' Fits efficiently and robustly Cox proportional hazards regression model in its
#' basic form, where explanatory variables are time independent with one event
#' per subject.  Method is based on a smooth modification of the partial
#' likelihood.
#'
#' @param formula a formula object, with the response on the left of a \code{~}
#' operator, and the terms on the right.  The response must be a
#' survival object as returned by the \code{\link{Surv}} function.
#' @param data a data frame in which to interpret the variables
#' named in the \code{formula}, or in the \code{subset}.
#' @param subset expression saying that only a subset of the rows of the data
#' should be used in the fit.
#' @param na.action a missing-data filter function, applied to the model.frame,
#' after any subset argument has been used.
#' @param trunc roughly, quantile of the sample \eqn{T_i exp(\beta'Z_i)},
#' it determines the trimming level for the robust estimator.
#' @param f.weight type of weighting function, default is \code{"quadratic"}
#' @param singular.ok logical value indicating how to handle collinearity in the
#' model matrix. If \code{TRUE}, the program will automatically skip
#' over columns of the X matrix that are linear combinations of
#' earlier columns.  In this case the coefficients for such
#' columns will be \code{NA}, and the variance matrix will contain
#' zeros.  For ancillary calculations, such as the linear
#' predictor, the missing coefficients are treated as zeros.
#' @param model a logical value indicating whether model frame should be
#' included as a component of the returned value.
#'
#' @references Bednarski, T. (1993). Robust estimation in Cox's regression model.
#' Scandinavian Journal of Statistics. Vol. 20, 213--225.
#' @references Bednarski, T. (1989). On sensitivity of Cox's estimator.
#' Statistics and Decisions. 7, 215--228.
#' @references Grzegorek, K.(1993). On robust estimation of baseline hazard under
#' the Cox model and via Frechet differentiability. Preprint of the Institute of
#' Mathematics of the Polish Academy of Sciences.518.
#' @references Minder, C.E. & Bednarski, T. (1996). A robust method for
#' proportional hazards regression. Statistics in Medicine Vol. 15, 1033--1047.
#'
#' @return a data frame containing MCMC summary statistics.An object of class
#' \code{coxr}. See \code{\link{coxr.object}} for details.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' \donttest{
#' if (interactive()) {
#' # Create a simple test data set using the attached function gen_data
#' a <- gen_data(200, c(1, 0.1, 2), cont = 0.05, p.censor = 0.30)
#' result <- coxr(Surv(time, status) ~ X1 + X2 + X3, data = a , trunc = 0.9)
#' result
#' plot(result)
#' }
#' }
#'
#' \dontshow{setwd(.old_wd)}
#' @export

coxr <- function(formula, data, subset, na.action, trunc = 0.95,
                f.weight = c("linear", "quadratic", "exponential"),
                singular.ok = TRUE, model = FALSE) {

    call <- match.call()

    mf <- match.call(expand.dots = FALSE)
    m  <- match(c("formula", "data", "subset", "na.action"), names(mf),
                nomatch = 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if ( NROW(mf) == 0 ) {
        stop("no (non-missing) observations")
    }

    mterms <- attr(mf, "terms")

    y <- model.extract(mf, "response")
    if ( !inherits(y, "Surv") ) {
        stop("response must be a \"Surv\" object")
    } else {
        type <- attr(y, "type")
        if ( type != "right" ) {
           stop(sprintf("\"%s\" type of survival data is not supported", type))
        }
    }

    x <- model.matrix(mterms, mf)
    x <- x[, -1, drop = FALSE]

    if ( missing(trunc) ) {
       trunc <- 0.95
    }

    if ( missing(f.weight) ) {
        f.weight <- "quadratic"
    } else {
        f.weight <- match.arg(f.weight)
    }
    f.weight = which(f.weight == c("linear", "quadratic", "exponential"))

    init <- rep(0,ncol(x))
    fit <- coxr.fit(x, y, trunc, init, f.weight, singular.ok)

    if ( length(fit$skip) > 0 ) {

        skipcol <- paste(fit$skip, collapse = ", ")

        msg <- sprintf("X matrix deemed to be singular; variable %s",
                        skipcol)

        if ( !singular.ok || length(fit$skip) == ncol(x) ) {
            stop(msg)
        } else {
            warning(msg)
        }

    }

    class(fit) <- "coxr"

    fit$na.action <- attr(mf, "na.action")
    fit$call <- call
    fit$terms <- mterms

    if (model) {
        fit$model <- mf
    }

    fit$x <- x
    fit$y <- y

    fit

}
