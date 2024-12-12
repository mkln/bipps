#' Get posterior::rvar object from output of multi_bipps
#'
#' @param out_list list of outputs from multi_bipps (one entry per chain)
#' @param variable string of variable to extract
#'
#' @return posterio::rvar of variable
#' @export
get_rvars <- function(out_list,variable,thin=1) {
  tryCatch({
    arr_list <- lapply(out_list,\(out) {
      arr <- out[[paste0(variable,"_mcmc")]]
      arr <- arr[,,seq(1,dim(arr)[3],by=thin),drop=FALSE]
      arr <- aperm(arr,perm=c(3,1,2))
      dims <- dim(arr)
      dims <- c(dims[1],1,dims[2:length(dims)])
      dim(arr) <- dims
      arr
    })
    arr <- abind::abind(arr_list,along=2)
    posterior::rvar(arr,with_chains = TRUE)
  },error=\(e) {
    stop(paste0("Variable ",variable," not found or can't convert to rvar"))
  })
}

#' Get cross-covariance matrix rvar at length-scale h
#'
#' Calculates the cross-covariance, e.g.,
#' \deqn{C{ij}(h) = \Lambda W \Lambda'}
#' with \eqn{\Lambda} estimated from MCMC and \eqn{W} the diagonal matrix with elements
#' \eqn{\exp{(-\phi_jh})} for \eqn{j=1,\ldots,k} on the diagonal.
#' @param out list of outputs from multi_bipps (one entry per chain)
#' @param h length scale (in [0,1])
#'
#' @return cross-covariance rvar object
#' @export
cross_covar <- function(out,h,...) {
  lambda <- get_rvars(out,"lambda",...)
  theta <- get_rvars(out,"theta",...)
  phi <- theta[1,]
  W <- diag_rfun(exp(-h*as.vector(phi)))
  covs <- lambda %*% W %*% t(lambda)
}

#' Get cross-correlation matrix rvar at length-scale h
#'
#' Calculates the cross-correlation matrix, e.g.,
#' \deqn{\rho_{ij}(h) = \frac{C_{ij}(h)}{\sqrt{C_{ii}(0)C_{jj}(0)}}}
#' with \eqn{C_{ij}(h)} the cross-covariance as calculated in [bipps::cross_covar].
#' @param out list of outputs from multi_bipps (one entry per chain)
#' @param h length scale (in \[0,1\])
#' @param cov0 the cross covariance calculated at h=0
#'
#' @return cross-correlation rvar object
#' @export
cross_cor <- function(out,h,cov0,...) {
  covs <- cross_covar(out,h,...)
  std_devs <- sqrt(diag(cov0))
  cross_cor <- sweep(sweep(covs, 1, std_devs, "/"), 2, std_devs, "/")
}

#' List of cross-cor/covar matrices at multiple length scales
#'
#' @param out list of outputs from multi_bipps (one entry per chain)
#' @param hs vector of length scales (in [0,1])
#' @param method either "cor" or "covar"
#'
#' @return list of cross matrices, one per value in hs
#' @export
cross_list <- function(out,hs,method="cor",...) {
  if(method == "cor") {
    cov0 <- cross_covar(out,0,...)
    l <- lapply(hs,\(h) {
      cross_cor(out,h,cov0,...)
    })
  } else {
    l <- lapply(hs,\(h) {
      cross_covar(out,h,...)
    })
  }


}



diag_rfun <- posterior::rfun(diag)
cov2cor_rfun <- posterior::rfun(cov2cor)

