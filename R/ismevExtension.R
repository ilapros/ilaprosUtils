
#' nicer print of gev.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param fitobj a fitted object of the class GEV
#' @keywords gev.fit
#' @export 
#' @examples
#' library(ismev)
#' print(gev.fit(c(50,45,65,78,12,23),show=FALSE))

print.gev.fit<-function(fitobj){
  zz<-list(mle=fitobj$mle,se=fitobj$se,conv=fitobj$conv,nllh=fitobj$nllh)
  print(zz)
}

#' nicer print of glo.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param fitobj a fitted object of the class GEV
#' @keywords glo.fit
#' @export 
#' @examples
#' library(ismev)
#' print(glo.fit(c(50,45,65,78,12,23),show=FALSE))
print.glo.fit<-function(fitobj){
  zz<-list(mle=fitobj$mle,se=fitobj$se,conv=fitobj$conv,nllh=fitobj$nllh)
  print(zz)
}

