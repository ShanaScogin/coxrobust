#' @useDynLib coxrobust ple
#' @useDynLib coxrobust lambda
#' @useDynLib coxrobust lin
#' @useDynLib coxrobust re

coxr.ple <- function(beta, time, status, z, nrow, ncol) {

  eps <- 1e-06
  iter.max <- 100
  l <- 0.8;

  for (i in 1:iter.max) {

    res <- .C("ple", as.double(beta), as.double(time),
                 as.integer(status), as.double(z), as.integer(nrow),
                 as.integer(ncol), res = double(1), gradient = double(ncol),
                 hessian = double(ncol*ncol),
                 PACKAGE = "coxrobust")

    hess <- matrix(res$hessian, ncol, ncol)
    hessinv <- chol2inv(chol(hess))
    beta_new <- beta - l*(res$gradient %*% hessinv)
    minimum_new <- res$res

    norm <- sqrt(sum((beta - beta_new)^2))

    if ( norm <= eps ) {
      beta <- beta_new
      res <- .C("ple", as.double(beta), as.double(time),
                   as.integer(status), as.double(z), as.integer(nrow),
                   as.integer(ncol), res = double(1), gradient = double(ncol),
                   hessian = double(ncol*ncol),
                   PACKAGE = "coxrobust")

      hess <- matrix(res$hessian, ncol, ncol)
      hessinv <- chol2inv(chol(hess))
      break
    } else if ( i>1 && minimum_new > minimum ) {
      beta <- (beta + beta_new) / 2
    } else {
      beta <- beta_new
      minimum <- minimum_new
    }

  }

  return(list(beta = as.double(beta), hess = hess, hessinv = hessinv))

}

coxr.re <- function(beta, time, status, z, prev_ezbeta, M, nrow, ncol, f.weight) {

  eps <- 1e-06
  iter.max <- 100
  l <- 0.8;

  for (i in 1:iter.max) {

    res <- .C("re", as.double(beta), as.double(time), as.integer(status),
                 as.double(z), as.double(prev_ezbeta), as.double(M),
                 as.integer(nrow),  as.integer(ncol), as.integer(f.weight),
                 res = double(1),  gradient = double(ncol),
                 hessian  = double(ncol*ncol),
                 PACKAGE = "coxrobust")

    hess <- matrix(res$hessian, ncol, ncol)
    hessinv <- chol2inv(chol(hess))
    beta_new    <- beta - l*(res$gradient %*% hessinv)
    minimum_new <- res$res

    norm_re <- sqrt(sum((beta - beta_new)^2))

    if ( norm_re <= eps ) {
      beta <- beta_new
      res <- .C("re", as.double(beta), as.double(time), as.integer(status),
                   as.double(z), as.double(prev_ezbeta), as.double(M),
                   as.integer(nrow),  as.integer(ncol), as.integer(f.weight),
                   res = double(1),  gradient = double(ncol),
                   hessian  = double(ncol*ncol),
                   PACKAGE = "coxrobust")

      hess <- matrix(res$hessian, ncol, ncol)
      hessinv <- chol2inv(chol(hess))
      break
    } else if ( i>1 && minimum_new > minimum ) {
      beta <- (beta + beta_new) / 2
    } else {
      beta <- beta_new
      minimum <- minimum_new
    }

  }

  return(list(beta = as.double(beta), hess = hess, hessinv = hessinv))

}

coxr.lambda <- function(beta, time, status, z, prev_exp_zbeta, M, nrow, f.weight) {

  cres <- .C("lambda", as.double(exp(z %*% beta)), as.double(time),
                as.integer(status), as.double(prev_exp_zbeta), as.double(M),
                as.integer(nrow), as.integer(f.weight), lmb = double(nrow),
                PACKAGE = "coxrobust")

  return(cres$lmb)

}

coxr.covar <- function(beta, time, status, z, prev_exp_zbeta, M, nrow, ncol,
                       f.weight, hessinv) {

  cres <- .C("lin", as.double( exp(z %*% beta) ), as.double(time),
                as.integer(status), as.double(z), as.double(prev_exp_zbeta),
                as.double(M), as.integer(nrow), as.integer(ncol),
                as.integer(f.weight), lin = double(nrow * ncol),
                PACKAGE = "coxrobust")

  lin <- matrix(cres$lin, nrow, ncol)

  IF <- lin %*% -hessinv

  cv <- matrix(0, nrow = ncol, ncol = ncol)

  for (i in 1:nrow) {
    cv <- cv + IF[i,] %*% t(IF[i,])
  }

  return( list(covar = cv) )

}
