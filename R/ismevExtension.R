
#' nicer print of gev.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param fitobj a fitted object of the class gev.fit
#' @keywords gev.fit
#' @examples
#' # library(ismev)
#' print(ismev::gev.fit(c(50,45,65,78,12,23),show=FALSE))
#' @export
#' @import ismev
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
#' a <- ismev::pp.fit(rain, 10)
#' a
#' @export
print.pp.fit<-function(obj){
  res<-list(mle=obj$mle,se=obj$se,nllh=obj$nllh)
  print(res)
}





#' @title Maximum Likelihood Fitting for the GLO distibution
#' @description This function has the same structure as the \code{gev.fit} from \code{ismev}
#' @param xdat A numeric vector of data to be fitted
#' @param ydat A matrix of covariates for generalized linear modelling of the parameters (or NULL (the default) for stationary fitting). The number of rows should be the same as the length of xdat
#' @param mul  Numeric vectors of integers, giving the columns of ydat that contain covariates for generalized linear modelling of the location parameter (or NULL (the default) if the corresponding parameter is stationary)
#' @param sigl As \code{mul} for the scale parameter
#' @param shl As \code{mul} for the shape parameter
#' @param mulink the link function for the location parameter - default to identity
#' @param siglink the link function for the scale parameter - default to identity
#' @param shlink the link function for the shape parameter - default to identity
#' @param muinit initial values for the location parameter
#' @param siginit initial values for the scale parameter
#' @param shinit initial values for the shape parameter
#' @param method The optimization method (see \code{optim} for details)
#' @param maxit  The maximum number of iterations 
#' @param ...   Other control parameters for the optimization. These are passed to components of the control argument of optim.
#' @return an object of the glo.fit class - with values which mirror the ones of the gev.fit class
#' @keywords glo.fit
#' @export
#' @examples
#' set.seed(5846)
#' print(glo.fit(rglo(n = 80, 50, 6, -0.2), show=FALSE))
glo.fit<-function(xdat, ydat = NULL, 
                  mul = NULL, sigl = NULL, shl = NULL, 
                  mulink = identity, siglink = identity, shlink = identity, 
                  muinit = NULL, siginit = NULL, shinit = NULL, 
                  show = TRUE, method = "Nelder-Mead", maxit = 10000, ...){
    #
    # obtains mles etc for glo distn
    #
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
      if(any(abs(xi) <10^-5)) {y <- (xdat - mu[1])/sc[1]}
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




#' @title Fitted return curve for a glo.fit object
#' @description This function mimics the \code{gev.rl} function for GLO models.
#' It plots the flood frequency curve (return curve) based on the data and the fitted parameters for a \code{glo.fit} model, including ones with historical data. 
#' Very ad-hoc and working under assumption that flood data are plotted (hence the y-axis lab).
#' Also outputs useful informations. 
#' Mostly written by Thomas Kjeldsen 
#' @param a the mle estimates from a \code{glo.fit} object
#' @param mat the covariance function from a \code{glo.fit} object
#' @param dat data matrix from a \code{glo.fit} object
#' @param nk if historicla data are used, the length of historical record - default is 0, no historical data
#' @param X0 if historicla data are used, the perception threshold - default is NULL, no historical data
#' @param f frequencies for which return level (and 95\%) confidence intervals are calculated
#' @export
#' @examples 
#' set.seed(7821567)
#' xx <- rglo(500, 40, 6, -0.2)
#' xxsist <- xx[471:500]; xxhist <- xx[1:470][xx[1:470] > 80]
#' h1 <- glo.hist.fit(c(xxhist,xxsist), k = length(xxhist), h = 470, X0 = 80, show=FALSE)
#' s1 <- glo.fit(xxsist, show=FALSE) 
#' rls1 <- glo.rl(a=s1$mle,mat=s1$cov,dat=s1$data)
#' rlh1 <- glo.rl(h1$mle,h1$cov,h1$data,nh=h1$h,nk=h1$k,X0=h1$X0)
#' lines(log(rlh1$f/(1-rlh1$f)), rlh1$rl+1.96*sqrt(rlh1$var), lty = 2)
#' lines(log(rlh1$f/(1-rlh1$f)), rlh1$rl-1.96*sqrt(rlh1$var), lty = 2)
#' lines(log(rls1$f/(1-rls1$f)), rls1$rl, col = 2)
#' lines(log(rls1$f/(1-rls1$f)), rls1$rl-1.96*sqrt(rls1$var), lty = 2, col = 2)
#' lines(log(rls1$f/(1-rls1$f)), rls1$rl+1.96*sqrt(rls1$var), lty = 2, col = 2)
#' legend("topleft", col =c(1,2),legend = c("With historical","Systematic Only"), lty = 1)
#' ## similar fitted curve - but obvious reduction in uncertainty 
glo.rl <- function(a, mat,dat,nh=0,nk=0,X0=NULL, f=c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8, 0.9, 0.95, 0.99, 0.995, 0.999)) {
    #
    # produces return level curve and 95 % confidence intervals 
    # on usual scale
    #
    # TRKJ, 18 December
    # a: maximum likelihood estimates of 3 parameters
    # mat: covariance matrix
    # dat: amax data time series
    # nh:  period covered by historical data
    # nk: number of historic data, 0 = no historic data
    # X0: NULL no perception threshold !NULL no perception threshold (actual values)
    #
    # f non-exceedance probablity P(X<x), thus log(f/(1-f)) = log(T-1), where T is return period
    #
    #
    if(is.null(X0)  & nh>0) stop("historic data with no perception threshold")
    if(!is.null(X0) & nh<1) stop("perception threshold with no historic data")
    eps <- 1e-006
    a1 <- a
    a2 <- a
    a3 <- a
    a1[1] <- a[1] + eps
    a2[2] <- a[2] + eps
    a3[3] <- a[3] + eps
    ## quantile function
    gloqq<-
      function (a, p) 
      {
        if (a[3] != 0) 
          a[1] + (a[2] * (1-((1 - p)/p)^(-a[3])))/a[3]
        else a[1] - a[2] * (log((1 - p)/p))
      }
    # Points on the quantile function
    q <- gloqq(a, 1 - f)
    d <- t(glo.rl.gradient( a=a, p=1-f))
    v <- apply(d, 1, q.form, m = mat); print(v)
    
    # set-up FEH-style plotting surface diagram
    par(mgp=c(2.75,1,0))
    plot(log(f/(1-f)), q,type = "n", 
         xlab = paste(paste("Logistic reduced variate","\n"),"Ln(T-1)"), 
         ylab = expression(paste("Peak flow ",(m^3/s))),
         bty="l",xaxt = "n")
    axis(1)
    # Plot return period axis
    # Axis is plotted at the level of the 1.11-year event (10th entry in the q-array)
    # Upper limit currently fixed at 200 years. Maybe should be a plotting option?
    #
    yaxis<-q[10]
    
    Taxis<-c(2,5,10,25,50,100,200)
    et_axis<-log(Taxis-1)
    paxis<-rep(yaxis,length(Taxis))
    
    #points(gLaxis,paxis,
    points(et_axis,paxis,
           pch="|",cex=0.75)
    lines(x=c(log(2-1),log(200-1)),y=c(yaxis,yaxis),
          lwd=1)
    
    rtp<-c("2","5","10","25","50","100","200")
    #text(gLaxis,paxis-5,rtp,cex=0.5)
    text(et_axis,paxis-0.20*yaxis,rtp,cex=0.75)
    text(4.0,1.5*yaxis,"Return period [years]",cex=0.75)
    
    
    # Plot FFC curves
    #
    lines(log(f/(1-f)), q)
    #   x<-log(f/(1-f))
    #   curve((a1+(a2/a3)*(1-exp(x)^(-a3))),from=1,to=499,n=500,add=TRUE)
    #     points(-1/(log((1:length(dat))/(length(dat) + 1))), sort(dat))
    #     lines(1/(1-f), q,col=2)
    #     lines(-1/log(f), q)
    #     lines(1/(1-f), q + 1.96 * sqrt(v), col = 4)
    #     lines(1/(1-f), q - 1.96 * sqrt(v), col = 4)
    #### change from ismev: use Gringorten plotting position
    # points(-1/(log((1:length(dat)-0.44)/(length(dat) + 1 - 0.88))), sort(dat),col=1)
    nex<-ifelse(is.null(X0),0,length(dat[dat>X0]))   # Number of observations > X0
    n<-length(dat)-nk+nh 
    fplot1<-NULL 
    # fplot1 is the nex points that exceed the perception threshold 
    # (both in systematic and hist data)
    if(nex!=0) fplot1<-((1:nex-0.44)/(nex + 1 - 0.88))*(nex/n)
    # flpot2 is the points that do not exceed the perception threshold in systematic data
    fplot2<- 
      (nex/n)+((n-nex)/n)*((seq((nex+1),length(dat))-nex-0.44)/(n-nh-length(dat[-seq(0,nk)][dat[-seq(0,nk)]>X0])+1-0.88))
    #    points(-1/log(1-c(fplot1,fplot2)),sort(dat,decreasing = TRUE),col=2,pch=19)
    if(is.null(X0)) X0<-max(dat)+1
    ##    points((-1/log(1-fplot2)),sort(dat[dat<X0],decreasing=TRUE),col=1,pch=1,cex=0.9)
    points(log((1-fplot2)/(fplot2)),sort(dat[dat<X0],decreasing=TRUE),col=1,pch=1,cex=0.9)
    ### we use a different symbol for the historical data 
    ### and it gets messy
    dat[1:nk]<-(dat[1:nk]+0.00001) ### to make sure we identify the right ones!!!
    histpoints<-seq(1,length(dat[dat>X0]))[(sort(dat[dat>X0]) %in% dat[1:nk])]
    if(nh>0){  ## if historical data exists
      #### this is not very pretty, but it is correct,it's all because of sorting issues
      ##      points(-1/log(1-sort(fplot1[length(dat[dat>X0])+1-histpoints])),
      fplot3<-(1-sort(fplot1[length(dat[dat>X0])+1-histpoints]))
      points(log(fplot3/(1-fplot3)),
             sort((sort(dat[dat>X0]))[histpoints],decreasing = TRUE),
             col=1,pch=17,cex=0.9)
      ##      points(-1/log(1-sort(fplot1[-(length(dat[dat>X0])+1-histpoints)])),
      fplot4<-1-sort(fplot1[-(length(dat[dat>X0])+1-histpoints)])
      cat(log(fplot4/(1-fplot4)))
      points(log(fplot4/(1-fplot4)),
             sort((sort(dat[dat>X0]))[-histpoints],decreasing = TRUE),
             col=1,pch=1,cex=0.9)
    }
    #    points(-1/(log((1:length(dat)-0.44)/(length(dat) + 1 - 0.88))), sort(dat),col=1)   
    ##    title("Return Level Plot")        
    ## add an output, so we can use the rl again
    res<-list(f=f,rl=q,var=v,dat=dat,n=length(dat),mle=a)
    res
}




#' @title GLO quantile function
#' @description This function gives quantiles the GLo distribution. It's built mirroring the \code{gev.fit} function of \code{ismev}. This is the same as qglo, just used to mimic ismev
#' @param loc,scale,sh the locatio, scale, shape parameters
#' @param p the non-exceedance probabiliy
gloq<-function (p, loc, scale, sh) {
  if (sh != 0) loc + (scale * (1-((1 - p)/p)^(-sh)))/sh
  else loc - scale * (log((1 - p)/p))
}


glo.rl.gradient<-function (a, p) {
  scale <- a[2]
  shape <- a[3]
  out <- matrix(NA, nrow = 3, ncol = length(p))
  out[1, ] <- 1
  
  yp <- p/(1 - p)
  out[2, ] <- shape^(-1) * (1 - yp^(shape))
  # out[3, ] <- scale * (shape^(-2)) * (1 - yp^(-shape)) -
  #   scale * shape^(-1) * yp^(-shape) * log(yp)
  out[3, ] <- - scale * (shape^(-2)) * (1 - yp^(shape)) -
    scale * shape^(-1) * yp^(shape) * log(yp)
  return(out)
}

