#' Get posterior::rvar object from output of multi_bipps
#'
#' @param out multi_bipps output
#' @param variable string of variable to extract
#'
#' @return posterio::rvar of variable
#' @export
get_rvars <- function(out,variable) {
  tryCatch({
    posterior::rvar(aperm(out[[paste0(variable,"_mcmc")]],perm=c(3,1,2)))
  },error=\(e) {
    stop(paste0("Variable ",variable," not found or can't convert to rvar"))
  })
}

#' Get cross-covariance matrix rvar at length-scale h
#'
#' @param out output from multi_bipps
#' @param h length scale (in [0,1])
#'
#' @return cross-covariance rvar object
#' @export
cross_covar <- function(out,h) {
  lambda <- get_rvars(out,"lambda")
  theta <- get_rvars(out,"theta")
  phi <- theta[1,]
  W <- diag_rfun(exp(-as.vector(phi)*h))
  covs <- lambda %*% W %*% t(lambda)
}

#' Get cross-correlation matrix rvar at length-scale h
#'
#' @param out output from multi_bipps
#' @param h length scale (in [0,1])
#'
#' @return cross-correlation rvar object
#' @export
cross_cor <- function(out,h) {
  covs <- cross_covar(out,h)
  cov2cor_rfun(covs)
}

#' List of cross-cor/covar matrices at multiple length scales
#'
#' @param out output from multi_bipps
#' @param hs vector of length scales (in [0,1])
#' @param method either "cor" or "covar"
#'
#' @return list of cross matrices, one per value in hs
#' @export
cross_list <- function(out,hs,method="cor") {
  if(method == "cor") {
    f <- cross_cor
  } else {
    f <- cross_covar
  }

  l <- lapply(hs,\(h) {
    f(out,h)
  })
}



diag_rfun <- posterior::rfun(diag)
cov2cor_rfun <- posterior::rfun(cov2cor)

