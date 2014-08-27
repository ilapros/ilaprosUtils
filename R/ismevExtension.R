
#' nicer print of gev.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param fitobj a fitted object of the class gev.fit
#' @keywords gev.fit
#' @examples
#' library(ismev)
#' print(gev.fit(c(50,45,65,78,12,23),show=FALSE))
#' @export
print.gev.fit<-function(fitobj){
  zz<-list(mle=fitobj$mle,se=fitobj$se,conv=fitobj$conv,nllh=fitobj$nllh)
  print(zz)
}

#' nicer print of glo.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param fitobj a fitted object of the class glo.fit
#' @keywords glo.fit
#' @examples
#' library(ismev)
#' print(glo.fit(c(50,45,65,78,12,23),show=FALSE))
#' @export
print.glo.fit<-function(fitobj){
  zz<-list(mle=fitobj$mle,se=fitobj$se,conv=fitobj$conv,nllh=fitobj$nllh)
  print(zz)
}




#' nicer print of pp.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param fitobj a fitted object of the class pp.fit
#' @keywords pp.fit
#' @examples
#' data(rain)
#' a <- pp.fit(rain, 10)
#' a
#' @export
print.pp.fit<-function(obj){
  res<-list(mle=obj$mle,se=obj$se,nllh=obj$nllh)
  print(res)
}




#' @title Fitting the GLo distribution
#' @description This function fits the GLo distribution. It's built mirroring the \code{gev.fit} function of \code{ismev}
#' @param xdat the object to be fitted
#' @keywords glo.fit
#' @export
#' @examples
#' set.seed(5846)
#' print(glo.fit(rglo(n = 80, 50, 6, -0.2), show=FALSE))
glo.fit<-function(xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, mulink = identity, siglink = identity, shlink = identity, muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, method = "Nelder-Mead", maxit = 10000, ...){
    #
    # obtains mles etc for glo distn
    #
    ##### this is based on the gev.fit function - can not actually come with covariates yet
    z <- list()
    npmu <- length(mul) + 1
    npsc <- length(sigl) + 1
    npsh <- length(shl) + 1
    z$trans <- FALSE # if maximization fails, could try
    # changing in1 and in2 which are
    # initial values for minimization routine
    in2 <- sqrt(6 * var(xdat))/pi
    in1 <- mean(xdat) - 0.57722 * in2
    if(is.null(mul)) {
      mumat <- as.matrix(rep(1, length(xdat)))
      if( is.null( muinit)) muinit <- in1
    }
    else {
      z$trans <- TRUE
      mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
      if( is.null( muinit)) muinit <- c(in1, rep(0, length(mul)))
    }
    if(is.null(sigl)) {
      sigmat <- as.matrix(rep(1, length(xdat)))
      if( is.null( siginit)) siginit <- in2
    }
    else {
      z$trans <- TRUE
      sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
      if( is.null( siginit)) siginit <- c(in2, rep(0, length(sigl)))
    }
    if(is.null(shl)) {
      shmat <- as.matrix(rep(1, length(xdat)))
      if( is.null( shinit)) shinit <- 0.1
    }
    else {
      z$trans <- TRUE
      shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
      if( is.null( shinit)) shinit <- c(0.1, rep(0, length(shl)))
    }
    z$model <- list(mul, sigl, shl)
    z$link <- deparse(substitute(c(mulink, siglink, shlink)))
    a<- init <- c(muinit, siginit, shinit)
    glo.lik <- function(a) {
      # computes neg log lik of glo model
      mu <- mulink(mumat %*% (a[1:npmu]))
      sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
      xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
      y <- ((1-xi*(xdat - mu)/sc))
      if(any(!y>0) || any(sc <= 0)) return(10^6)
      y<- -log(y)/xi
      (sum(log(sc)) + sum(y*(1-xi)) + 2*sum(log(1+exp(-y) )))
    }
    x <- optim(init, glo.lik, hessian = TRUE, method = method,
               control = list(maxit = maxit, ...))
    z$conv <- x$convergence
    mu <- mulink(mumat %*% (x$par[1:npmu]))
    sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
    xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
    z$nllh <- x$value
    z$data <- xdat
    if(z$trans) {
      z$data <- - log(as.vector((1 + (xi * (xdat - mu))/sc)^( -1/xi)))
    }
    z$mle <- x$par
    z$cov <- solve(x$hessian)
    z$se <- sqrt(diag(z$cov))
    z$vals <- cbind(mu, sc, xi)
    if(show) {
      if(z$trans)
        print(z[c(2, 3, 4)])
      else print(z[4])
      if(!z$conv)
        print(z[c(5, 7, 9)])
    }
    class(z) <- "glo.fit"
    invisible(z)
}






#' @title GLO quantile function
#' @description This function fits the GLo distribution. It's built mirroring the \code{gev.fit} function of \code{ismev}. This is the same as qglo, just used to mimic ismev
#' @param loc,scale,sh the locatio, scale, shape parameters
#' @param p the non-exceedance probabiliy
gloq<-function (p, loc, scale, sh) {
  if (sh != 0) loc + (scale * (1-((1 - p)/p)^(-sh)))/sh
  else loc - scale * (log((1 - p)/p))
}


glo.rl.gradient<-function (loc, scale, sh, p) {
  scale <- a[2]
  shape <- a[3]
  out <- matrix(NA, nrow = 3, ncol = length(p))
  out[1, ] <- 1
  
  yp <- p/(1 - p)
  out[2, ] <- shape^(-1) * (1 - yp^(shape))
  out[3, ] <- scale * (shape^(-2)) * (1 - yp^(-shape)) -
    scale * shape^(-1) * yp^(-shape) * log(yp)
  
  return(out)
}

