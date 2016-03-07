



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







#' @title L-moment estimation in the presence of Historical Data
#' @description L-moment fitting procedures for sample of both systematic and historical data via Partial Probability Weighted Moments (Wang, 1990).
#' @param xdat vector of historical and systematic/observed data - the first k elements of the vector should be the historical events
#' @param k number of historical events available. These events should be stored as the first observations of the xdat vector
#' @param h length of years covered by the information of the historical period
#' @param X0 the perception threshold which is exceeded by the historical events. This is most likely different from the lowest historical value. It should indicate a value after which we are confident the event would have been recorded or left traces 
#' @param nmom the number of moments to estimate
#' @return A vector of L-moment, as lmom::samlmu.
#' @references Wang, Q.J.(1990). Unbiased estimation of probability weighted moments and partial probability weighted moments from systematic and historical flood information and their application. Journal of hydrology, 120, 115--124.
#' @export
#' @examples 
#' # library(ismev)
#' set.seed(54165784)
#' xx <- rgev(500, 40, 6, -0.2)
#' xxsist <- xx[471:500]; xxhist <- xx[1:470][xx[1:470] > 80]
#' lmh <- lmom.hist.fit(c(xxhist,xxsist), k = length(xxhist), h = 470, X0 = 80, nmom = 5)
#' lmh
#' lmom::pelgev(lmh)
#' gev.hist.fit(c(xxhist,xxsist), k = length(xxhist), h = 470, X0 = 80, show = FALSE)$mle
#' ### note the different parametrization for the shape parameter
lmom.hist.fit <- function(xdat, k = NULL, h = NULL, X0 = NULL, nmom = 3){
  if((is.null(k) + is.null(X0) + is.null(h)) %in% c(1,2) ) stop("I need all information (X0, k, h) on historical records")
  if(all(is.null(k), is.null(X0), is.null(h))) {
    lmom::samlmu(xdat, nmom = nmom)
  }
  else {
    extSamlmu(wangPPWM(xcont = xdat[(k+1):(length(xdat))], xhist = xdat[1:k], 
                       nbans = h, seuil = X0, nmom = nmom))
  }
}



wangPPWM <- function(xcont, xhist, nbans, seuil, nmom = 3){ 
  ## naming after the BayesMCMC function (which doesn't make sense)
  n <- length(xcont) ## length of systematic record
  m <- n + nbans     ## length of systematic+historical record 
  l <- sum(xcont > seuil) + sum(xhist>seuil) ## number of events above threshold 
  xcont <- sort(xcont)
  xhist <- sort(xhist)
  F0hat <- 1-l/m
  aboveThresh <- sort(c(xcont[xcont > seuil], xhist[xhist > seuil]))
  if(any(xhist <= seuil)) stop("all historical events should be above the threshold")
  if(l != length(aboveThresh)) stop("something non-matching in length of obs. above threshold")
  ### the betas
  beta.second <- rep(0,nmom)
  for(r in 1:nmom){
    for(i in 1:n){
      beta.second[r] <- beta.second[r] +
        choose(i-1, (r-1)) * ifelse(xcont[i] <= seuil, xcont[i], 0)
    } 
    beta.second[r] <- ((n * choose(n-1, (r-1)))^(-1)) * beta.second[r]
  } 
  beta.prime <- rep(0,nmom)
  for(r in 1:nmom){
    for(i in seq(from = (m-l+1), to = m)){
      ## i goes to the high number -
      ## we do not have the values below the threshold for the m historical years  
      beta.prime[r] <- beta.prime[r] +
        choose(i-1, (r-1)) * aboveThresh[i-(m-l)]
    }
    beta.prime[r] <- ((m * choose(m-1, (r-1)))^(-1)) * beta.prime[r]
  }
  #   print(c(F0hat, 1/(1-F0hat), beta.second,beta.prime))
  #   beta.prime <- c(1/(1-F0hat),1/(1-F0hat^2), 1/(1-F0hat^3)) * beta.prime
  #   beta.second <- c(1/(1-F0hat),1/(1-F0hat^2), 1/(1-F0hat^3)) * beta.second
  list(beta.prime = beta.prime, beta.second = beta.second, 
       betas = beta.prime+beta.second)
}


extSamlmu <- function(x) {
  firstOut <- lmomcopwm2lmom(x)
  out <- c(firstOut$lambdas[1:2],firstOut$ratios[3:length(firstOut$ratios)])
  names(out) <- c("l_1","l_2", paste("t",seq(3,length(firstOut$ratios), by = 1),sep = "_")) 
  out
}



### copy of lmomco::pwm2lmom - avoid having to depend on the package
lmomcopwm2lmom <- function (pwm){
  if (is.list(pwm)) {
    if (!is.null(pwm$BETA0)) {
      z <- list(L1 = NULL, L2 = NULL, TAU3 = NULL, TAU4 = NULL, 
                TAU5 = NULL, LCV = NULL, L3 = NULL, L4 = NULL, 
                L5 = NULL)
      z$L1 <- pwm$BETA0
      z$L2 <- 2 * pwm$BETA1 - pwm$BETA0
      z$L3 <- 6 * pwm$BETA2 - 6 * pwm$BETA1 + pwm$BETA0
      z$L4 <- 20 * pwm$BETA3 - 30 * pwm$BETA2 + 12 * pwm$BETA1 - 
        pwm$BETA0
      z$L5 <- 70 * pwm$BETA4 - 140 * pwm$BETA3 + 90 * pwm$BETA2 - 
        20 * pwm$BETA1 + pwm$BETA0
      z$LCV <- z$L2/z$L1
      z$TAU3 <- z$L3/z$L2
      z$TAU4 <- z$L4/z$L2
      z$TAU5 <- z$L5/z$L2
      return(z)
    }
    if (!is.null(pwm$betas)) {
      pwm <- pwm$betas
    }
    else {
      warning("ambiguous call, do not find Betas for processing")
      return(NULL)
    }
  }
  nmom <- length(pwm)
  L <- vector(mode = "numeric", length = nmom)
  R <- vector(mode = "numeric", length = nmom)
  for (i in seq(1, nmom)) {
    r <- i - 1
    sum <- 0
    for (k in seq(0, r)) {
      weight <- (-1)^(r - k) * choose(r, k) * choose(r + 
                                                       k, k)
      sum <- sum + weight * pwm[k + 1]
    }
    L[i] <- sum
  }
  if (nmom >= 2) {
    R[2] <- L[2]/L[1]
  }
  if (nmom >= 3) {
    for (r in seq(3, nmom)) {
      R[r] <- L[r]/L[2]
    }
  }
  R[1] <- NA
  z <- list(lambdas = L, ratios = R)
  return(z)
}






