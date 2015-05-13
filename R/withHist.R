



#' @title Maximum Likelihood Fitting for the GEV distibution in the presence of Historical Data
#' @description ML fitting procedures for sample of both systematic and historical data. Function structure inspired by \code{gev.fit} in ismev. 
#' Arguments naming and likelihood implementation follow Stedinger and Cohn (1986). 
#' @param xdat vector of historical and systematic/observed data - the first k elements of the vector should be the historical events
#' @param k number of historical events available. These events should be stored as the first observations of the xdat vector
#' @param h length of years covered by the information of the historical period
#' @param X0 the perception threshold which is exceeded by the historical events. This is most likely different from the lowest historical value. It should indicate a value after which we are confident the event would have been recorded or left traces 
#' @param binomialcens a logical value. Indicates whether the actual k values are to be used, or if only the information that the threshold X0 has been exceeded is used.
#' @param mulink the link function for the location parameter - default to identity
#' @param siglink the link function for the scale parameter - default to identity
#' @param shlink the link function for the shape parameter - default to identity
#' @param muinit initial values for the location parameter
#' @param siginit initial values for the scale parameter
#' @param shinit initial values for the shape parameter
#' @param method The optimization method (see \code{optim} for details)
#' @param maxit  The maximum number of iterations 
#' @param ...   Other control parameters for the optimization. These are passed to components of the control argument of optim.
#' @details The form of the GEV used is that of Coles (2001) Eq (3.2). Specifically, positive values of the shape parameter imply a heavy tail, and negative values imply a bounded upper tail.
#' @return An object of the class \code{gev.fit}. 
#' @export
#' @examples 
#' set.seed(7821567)
#' xx <- rglo(500, 40, 6, -0.2)
#' xxsist <- xx[471:500]; xxhist <- xx[1:470][xx[1:470] > 80]
#' glo.hist.fit(c(xxhist,xxsist), k = length(xxhist), h = 470, X0 = 80)
#' glo.hist.fit(c(xxhist,xxsist), k = length(xxhist), h = 470, X0 = 80, binomialcens = TRUE)
#' glo.fit(xxsist) ## notice the higher standard errors
glo.hist.fit<-
  function(xdat,k=0,h=NULL,X0=NULL,binomialcens=FALSE,ydat = NULL, 
           mul = NULL, sigl = NULL, shl = NULL, 
           mulink = identity, siglink = identity, shlink = identity, 
           muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, 
           method = "Nelder-Mead", maxit = 10000, ...)
  {
    #
    # obtains mles etc for glo dist when historical data are present
    #
    ### can not actually come with covariates yet
    ## binomialcens=TRUE if a binomial type of censoring is applied, 
    ## eg you only consider that some extremes happened in the past, not their dimension
    ## h is the number of past years you are covering with your historical data
    ## k is the number of event you see in the past
    ## X0 is the perception threshold
    if(k == 0 & !(is.null(h) & is.null(X0))) {warning("no historical data given, dropping X0 and h"); h=0; X0 = NULL} 
    z <- list()
    npmu <- length(mul) + 1
    npsc <- length(sigl) + 1
    npsh <- length(shl) + 1
    z$trans <- FALSE  # if maximization fails, could try
    # changing in1 and in2 which are 
    # initial values for minimization routine
    in2 <- sqrt(6 * var(xdat[(k+1):(length(xdat))]))/pi
    in1 <- mean(xdat[(k+1):(length(xdat))]) - 0.57722 * in2
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
    if(is.null(X0)) {
      X0<-min(xdat[1:k])
    }  
    z$model <- list(mul, sigl, shl)
    z$link <- deparse(substitute(c(mulink, siglink, shlink)))
    a<- init <- c(muinit, siginit, shinit)
    glo.hist.lik <- function(a,h,k,binomialcens) {
      # computes neg log lik of gev model
      mu <- mulink(mumat %*% (a[1:npmu]))
      sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
      xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
      y <- ((1-xi*(xdat - mu)/sc))
      if(any(!y>0) || any(sc <= 0)) return(10^6)
      y <- -log(y)/xi;y0<- -log(((1-xi[1]*(X0 - mu[1])/sc[1])))/xi[1]
      ## shape = 0 case
      if(any(abs(xi) <10^-5)) {y <- (xdat - mu[1])/sc[1]; y0<- (X0 - mu[1])/sc[1]}
      if(k==0) return(sum(log(sc)) + sum(y*(1-xi)) + 2*sum(log(1+exp(-y))))
      if(!binomialcens & !k==0) return(sum(log(sc)) + sum(y*(1-xi)) + 2*sum(log(1+exp(-y)))-log(choose(h,k))-(h-k)*log((1/(1+exp(-y0)))) )
      # If only historical exceedances are known, then below
      #      if(!binomialcens) return(-log(choose(h,k))-(h-k)*log((1/(1+exp(-y0)))) -k*log((1-(1/(1+exp(-y0)))) ))
      if(binomialcens & !k==0) return(sum(log(sc[-(1:k)])) + sum(y[-(1:k)]*(1-xi[-(1:k)])) + 2*sum(log(1+exp(-y[-(1:k)]))) -log(choose(h,k))-(h-k)*log((1/(1+exp(-y0)))) -k*log((1-(1/(1+exp(-y0)))))   )
    }
    x <- optim(init, glo.hist.lik, k=k, h=h, binomialcens=binomialcens, hessian = TRUE, method = method,
               control = list(maxit = maxit, ...))
    z$conv <- x$convergence
    mu <- mulink(mumat %*% (x$par[1:npmu]))
    sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
    xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
    z$nllh <- x$value
    z$data <- xdat
    if(z$trans) {
      z$data <-  - log(as.vector((1 + (xi * (xdat - mu))/sc)^(
        -1/xi)))
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
    z$k<-k;z$h<-h;z$X0<-X0
    class(z) <- "glo.fit"
    invisible(z)
}


#' @title Maximum Likelihood Fitting for the GEV distibution in the presence of Historical Data
#' @description ML fitting procedures for sample of both systematic and historical data. Function structure inspired by \code{gev.fit} in ismev. 
#' Arguments naming and likelihood implementation follow Stedinger and Cohn (1986). 
#' @param xdat vector of historical and systematic/observed data - the first k elements of the vector should be the historical events
#' @param k number of historical events available. These events should be stored as the first observations of the xdat vector
#' @param h length of years covered by the information of the historical period
#' @param X0 the perception threshold which is exceeded by the historical events. This is most likely different from the lowest historical value. It should indicate a value after which we are confident the event would have been recorded or left traces 
#' @param binomialcens a logical value. Indicates whether the actual k values are to be used, or if only the information that the threshold X0 has been exceeded is used.
#' @param mulink the link function for the location parameter - default to identity
#' @param siglink the link function for the scale parameter - default to identity
#' @param shlink the link function for the shape parameter - default to identity
#' @param muinit initial values for the location parameter
#' @param siginit initial values for the scale parameter
#' @param shinit initial values for the shape parameter
#' @param method The optimization method (see \code{optim} for details)
#' @param maxit  The maximum number of iterations 
#' @param ...	 Other control parameters for the optimization. These are passed to components of the control argument of optim.
#' @details The form of the GEV used is that of Coles (2001) Eq (3.2). Specifically, positive values of the shape parameter imply a heavy tail, and negative values imply a bounded upper tail.
#' @return An object of the class \code{gev.fit}. 
#' @export
#' @examples 
#' # library(ismev)
#' set.seed(5416574)
#' xx <- rgev(500, 40, 6, -0.2)
#' xxsist <- xx[471:500]; xxhist <- xx[1:470][xx[1:470] > 80]
#' gev.hist.fit(c(xxhist,xxsist), k = length(xxhist), h = 470, X0 = 80)
#' gev.hist.fit(c(xxhist,xxsist), k = length(xxhist), h = 470, X0 = 80, binomialcens = TRUE)
#' ismev::gev.fit(xxsist) # note the higher standard errors  
gev.hist.fit <-
  function(xdat, k=0, h=NULL, X0=NULL, binomialcens=FALSE, 
           mulink = identity, siglink = identity, shlink = identity, 
           muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, 
           method = "Nelder-Mead", maxit = 10000, ...){
  #
  # obtains mles etc for gev dist when historical data are present
  ## builds on ismev::gev.fit
  ### can not actually come with covariates yet
  ## binomialcens=TRUE if a binomial type of censoring is applied, 
  ## eg you only consider that some extremes happened in the past, not their dimension
  ## h is the number of past years you are covering with your historical data
  ## k is the number of event you see in the past
  ## X0 is the perception threshold 
  z <- list()
  if(k == 0 & !(is.null(h) & is.null(X0))) {warning("no historical data given, dropping X0 and h"); h=0; X0 = NULL} 
  mul = NULL; sigl = NULL; shl = NULL
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  z$trans <- FALSE
  in2 <- sqrt(6 * var(xdat[(k+1):(length(xdat))]))/pi
  in1 <- mean(xdat[(k+1):(length(xdat))]) - 0.57722 * in2
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, length(xdat)))
    if (is.null(muinit)) 
      muinit <- in1
  }
  else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
    if (is.null(muinit)) 
      muinit <- c(in1, rep(0, length(mul)))
  }
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, length(xdat)))
    if (is.null(siginit)) 
      siginit <- in2
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
    if (is.null(siginit)) 
      siginit <- c(in2, rep(0, length(sigl)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdat)))
    if (is.null(shinit)) 
      shinit <- 0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
    if (is.null(shinit)) 
      shinit <- c(0.1, rep(0, length(shl)))
  }
  z$model <- list(mul, sigl, shl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink)))
  init <- c(muinit, siginit, shinit)
  
  gev.hist.lik <- function(a,h,k,binomialcens) {
    # computes neg log lik of gev model
    # uses COles' formulation of the F(x) and f(x) of a GEV
    # the shape parameter is the negative of the one used in Hosking and Wallis 
    # so it must be interpreted in the opposite way than for the GLO
    # this is a mess, I am sorry... 
    mu <- mulink(mumat %*% (a[1:npmu]))
    sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
  y <- (xdat - mu)/sc
  y <- 1 + xi * y
  if(any(y <= 0) || any(sc <= 0)) 
      return(10^6)
  if(k == 0) return(sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi +  1)))
  y0 <- 1 + xi[1] * (X0 - mu[1])/sc[1]
  if(!binomialcens & !k==0) return(sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi +  1)) -
                                   log(choose(h,k))-(h-k)*(-y0^(-1/xi[1])))
  if(binomialcens & !k==0) return(sum(log(sc[-(1:k)])) + sum(y[-(1:k)]^(-1/xi[-(1:k)])) + sum(log(y[-(1:k)]) * (1/xi[-(1:k)] +  1)) -
                                     log(choose(h,k))-(h-k)*(-y0^(-1/xi[1])) - k*log(1-exp(-y0^(-1/xi[1]))))
  }
  x <- optim(init, gev.hist.lik, k=k, h=h, binomialcens=binomialcens, hessian = TRUE, method = method,
             control = list(maxit = maxit, ...))
  z$conv <- x$convergence
  mu <- mulink(mumat %*% (x$par[1:npmu]))
  sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
  xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
  z$nllh <- x$value
  z$data <- xdat
  if (z$trans) {
    z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
  }
  z$mle <- x$par
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  z$vals <- cbind(mu, sc, xi)
  z$k <- z$k; z$h <- z$h; z$X0 <- z$X0
  if (show) {
    if (z$trans) 
      print(z[c(2, 3, 4)])
    else print(z[4])
    if (!z$conv) 
      print(z[c(5, 7, 9)])
  }
  z$k <- z$k; z$h <- z$h; z$X0 <- z$X0
  class(z) <- "gev.fit"
  invisible(z)
}
