
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



#' nicer print of gpd.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param fitobj a fitted object of the class pp.fit
#' @keywords gpd.fit
#' @examples
#' a <- ismev::gpd.fit(c(53, 52, 49, 58, 50, 48, 47, 50, 46, 46, 49, 51, 47, 49, 50), threshold = 46, show=FALSE)
#' a
#' @export
print.gpd.fit<-function(fitobj){
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
print.pp.fit<-function(fitobj){
  res<-list(mle=fitobj$mle,se=fitobj$se,conv=fitobj$conv,nllh=fitobj$nllh)
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
      if( is.null(siginit)) siginit <- c(in2, rep(0, length(sigl)))
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




gumbX <- function(x, emp=FALSE) {
  if(emp) x <- (1:length(x)-0.44)/(length(x) + 1 - 0.88)
  -log(-log(x))
}

redVarX <- function(x, emp=FALSE) {
  if(emp) x <- (1:length(x)-0.44)/(length(x) + 1 - 0.88)
  log((x)/(1-x))
}


#' @title Plots of return curves and confidence intervals estimated by the delta method
#' @description This function is based on the \code{ismev::gev.rl} function but also alows the case in which historical data are added to the systematic record.
#' It plots the return curve (flood frequency curve in hydrology) based on the data and the fitted parameters for a \code{glo.fit} or \code{gev.fit} object, including ones with historical data. 
#' The data points are plotted using the Gringorten plotting positions. 
#' The function also outputs useful the estimated return levels and the corresponding standard error as estimated via the delta method (see Coles, 2001). 
#' @param obj a \code{glo.fit} or \code{gev.fit} object
#' @param p non-exceedance probabilities for which return level and standard errors are calculated. If the given non-exceedance probability don't cover the whole range of empirical non-exceedance probailities of the data the figure will automatically draw a line covering the whole range, but the output dataframe will only contain the specified non-exceedance probabilities.
#' @param vtype a character specifying the reduced variate type to be used in the figure. The types allowed are "Gumbel", corresponding to -log(-log(p)), and "redVar", corresponding to log(p/(1-p)). The dafault is set to "redVar".
#' @param sign.alpha significance level required for the confidence intervals - default to 0.05
#' @param pchHist pch parameter to be used for the historical data (if present) - default to 15
#' @param plot.out logical, indicating whether the plot should actually be displayed; set to FALSE to only compute the return levels and standard errors
#' @param ...	Arguments to be passed to methods, such as graphical parameters (see par)
#' @return a return levels figure and a data.frame containing the estimated return levels and corresponding standard errors for the specified exceedance probabilities
#' @export
#' @examples 
#' set.seed(7821567)
#' xx <- rglo(500, 40, 6, -0.2)
#' xxsist <- xx[471:500]; xxhist <- xx[1:470][xx[1:470] > 80]
#' s1 <- glo.fit(xxsist, show=FALSE) 
#' rls1 <- retPlot(s1, sign.alpha = 0.1, col = 4, p = c(seq(0.01,0.991,length=45),seq(0.992,0.9992,length=120)))
#' h1 <- glo.hist.fit(c(xxhist,xxsist), k = length(xxhist), h = 470, X0 = 80, show=FALSE)
#' rlh1 <- retPlot(h1, vtype = "Gumbel", col = 1, sign.alpha = 0.05, p = rls1$p, xlab = "Gumbel reduced variate (-log(-log(1-1/T)))")
#' lines(-log(-log(rls1$p)), rls1$retLev, col = 2)
#' lines(-log(-log(rls1$p)), rls1$retLev+qnorm(0.025)*rls1$se, lty = 2, col = 2)
#' lines(-log(-log(rls1$p)), rls1$retLev-qnorm(0.025)*rls1$se, lty = 2, col = 2)
#' legend("topleft", col =c(1,2),legend = c("With historical","Systematic Only"), lty = 1)
#' ## similar fitted curve - but large reduction in uncertainty for rare events
retPlot <- function(obj, p=NULL, 
                    sign.alpha = 0.05,
                    plot.out = TRUE,
                    vtype = "redVar", 
                    ylim = NULL, xlim = NULL, pch = par()$pch, pchHist = 15, ...) { # , gridProb = NULL, gridLab,
  a <- as.vector(obj$mle); mat <- as.matrix(obj$cov); dat <- as.vector(obj$data)
  nh <- 0; nk <- 0; X0 <- 0
  if(exists("X0",obj)){
    nh <- obj$h; nk <- obj$k; X0 <- obj$X0
  }
  a <- a + 1e-006
  if(is.null(p)) p <- c(0.01, 0.02, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8, 0.9, 0.95, 0.99, 0.991, seq(0.992, 0.995, by = 0.001))
  f <- 1-p
  fline <- f
  ### function to calculate the right type of reduced variate
  xcalc <- function(x, vtype, emp = FALSE) switch(vtype,"Gumbel" = gumbX(x, emp = emp),"redVar" = redVarX(x, emp = emp))
  
  xemp <- xcalc(c(dat,rep(mean(dat),nh)), vtype = vtype, emp = TRUE)
  if(max(xcalc(1-f, vtype = vtype)) < max(xemp)) f <- sort(c(f, seq(1/((length(dat)+nh)), 1/(2*(length(dat)+nh)),l=120)))
  xv <- xcalc(1-f, vtype = vtype)

  ## quantiles
  y <- switch(class(obj),
    "glo.fit" = gloq(f,loc = a[1],scale = a[2],sh = a[3]),
    "gev.fit" = ismev::gevq(p = f, a = a))

  # Confidence intervals by means of delta method
  v <- switch (class(obj),
                  "glo.fit" = apply(t(glo.rl.gradient(a=a, p=f)), 1, ismev::q.form, m = mat),
                  "gev.fit" = apply(t(ismev::gev.rl.gradient(a=a, p=f)), 1, ismev::q.form, m = mat))
  # v <- v[seq(length(v),1)]
  if(plot.out){
  if(is.null(ylim)) ylim <- range(c(y, dat))
  if(is.null(xlim)) xlim <- range(c(xv, xemp))
  plot(xv, y, type = "n", ylim = ylim, xlim = xlim, ...)
  lines(xv, y, ...)
  lines(xv, y + qnorm(1-sign.alpha/2) * sqrt(v), ...)
  lines(xv, y - qnorm(1-sign.alpha/2) * sqrt(v), ...)
  if(!(nh > 0)){
    points(xemp, sort(dat), pch = pch, ...)
  }
  if(nh > 0){
    ### following Bayliss and Reed
    nex <- length(dat[dat>X0])   # Number of observations > X0 (both in systematic and hist data when present)
    n<-length(dat)-nk+nh 
    fplot1<-NULL 
    # fplot1 is the nex points that exceed the perception threshold 
    # (both in systematic and hist data)
    if(nex!=0) fplot1<-((1:nex-0.44)/(nex + 1 - 0.88))*(nex/n)
    # flpot2 is the points that do not exceed the perception threshold in systematic data
    fplot2 <- 
      (nex/n)+((n-nex)/n)*((seq((nex+1),length(dat))-nex-0.44)/(n-nh-length(dat[-seq(0,nk)][dat[-seq(0,nk)]>X0])+1-0.88))
    points(xcalc(1-fplot2, vtype = vtype),sort(dat[dat<X0], decreasing=TRUE), pch = pch, ...)
    ### we use a different symbol for the historical data 
    ### and it gets messy
    dat[1:nk]<-(dat[1:nk]+0.00001) ### to make sure we identify the right ones!!!
    histpoints<-seq(1,length(dat[dat>X0]))[(sort(dat[dat>X0]) %in% dat[1:nk])]
    fplot3 <- sort(fplot1[length(dat[dat>X0])+1-histpoints])
    points(xcalc(1-fplot3, vtype = vtype),
           sort((sort(dat[dat>X0]))[histpoints],decreasing = TRUE), pch=pchHist, ... )
    fplot4 <- sort(fplot1[-(length(dat[dat>X0])+1-histpoints)])
    points(xcalc(1-fplot4, vtype = vtype),
           sort((sort(dat[dat>X0]))[-histpoints],decreasing = TRUE),pch = pch, ...)
  }
  }
  invisible(data.frame(p = 1-f, retLev = y, se = sqrt(v))[f %in% fline, ])
}


# xl <- ifelse(any(names(list(...)) == "xlab"), 
#              list(...)["xlab"],  
#              switch(vtype,
#                     "Gumbel" = "Gumbel variate: -log(-log(p))",
#                     "redVar" = "Reduced variate: log(p/(1-p))"))


#' @title Fitted return curve for a glo.fit object
#' @description This function mimics the \code{gev.rl} function for GLO models.
#' It plots the flood frequency curve (return curve) based on the data and the fitted parameters for a \code{glo.fit} model, including ones with historical data. 
#' Very ad-hoc and working under assumption that flood data are plotted (hence the y-axis lab).
#' Also outputs useful informations. 
#' Mostly written by Thomas Kjeldsen. 
#' This function is deprecated and has been replaced by the generic retPlot function.
#' @param a the mle estimates from a \code{glo.fit} object
#' @param mat the covariance function from a \code{glo.fit} object
#' @param dat data matrix from a \code{glo.fit} object
#' @param nk if historical data are used, the length of historical record - default is 0, no historical data
#' @param X0 if historical data are used, the perception threshold - default is NULL, no historical data
#' @param f frequencies for which return level (and 95\%) confidence intervals are calculated
#' @export
#' @seealso \code{\link{retPlot}}
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

