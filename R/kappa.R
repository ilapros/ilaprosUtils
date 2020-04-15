#### kappa as in Hosking and Wallis 


#----- the KAPPA distributon -------
#' @title The Kappa distribution
#' @description Density, distribution function, quantile function and random generation 
#' for the Kappa distribution with location parameter equal to \code{loc}, scale parameter equals 
#' to \code{scale}, first shape parameter equal to \code{sh} and second shape parameter equal to \code{sh2}.
#' @param loc location parameter
#' @param scale scale parameter
#' @param sh first shape parameter
#' @param sh2 second shape parameter
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}
#' @keywords Kappa distribution
#' @aliases pkappa dkappa rkappa qkappa
#' @name Kappa distribution 
#' @details#' The distribution is described in the Hosking and Wallis book and can be seen as a generelisation of several distributions used in extreme value modelling. 
#' Depending on the sh2 parameter value the Kappa distribution reduces to others commonly used distributions For example: 
#' when sh2 is equal to -1 the distribution reduces to a GLO, when sh2 is equal to 0 the distribution reduces to a GEV, 
#' when sh2 is equal to +1 the distribution reduces to a Generalised Pareto distribution (GPA).  
#' @return dkappa gives the density, pkappa gives the distribution function, 
#' qkappa gives the quantile function, and rkappa generates random deviates.
#' The length of the result is determined by n for rkappa, and is the maximum of the lengths of the numerical arguments for the other functions.
#' The numerical arguments are recycled to the length of the result. 
#' Only the first elements of the logical arguments are used.
#' @family kappa distribution
#' @export
#' @references 
#' Hosking, J.R.M. and Wallis, J.R., 2005. Regional frequency analysis: an approach based on L-moments. Cambridge university press. 
#' @examples 
#' plot(seq(-26,40,by=0.2),dkappa(seq(-26,40,by=0.2),4,6,0.2,-0.4),type="l", ylab = "density")
#' lines(seq(-26,40,by=0.2),dkappa(seq(-26,40,by=0.2),4,6,0.2,-1),type="l", col = 2)
#' set.seed(123)
#' plot(ecdf(rkappa(100,4,6,0.2,-0.4)))
#' lines(seq(-20,30,by=0.5),pkappa(seq(-20,30,by=0.5),4,6,0.2,-0.4),col=2)
#' ## notable quantiles 
#' qkappa(c(0.5,0.99,0.995,0.995,0.999),4,6,0.2, -0.4) 
pkappa <- function(q, loc, scale, sh, sh2, lower.tail = TRUE) {
  ## in HW notation 
  allTogether <- cbind(q,loc,scale,sh,sh2)
  apply(allTogether,1,
        pkappa_int,lower.tail=lower.tail)
}


pkappa_int <- function(x, lower.tail = TRUE) {
  q = x[1]; loc = x[2]; scale= x[3]; sh= x[4]; sh2 = x[5]
  ## in HW notation 
  xi = loc;  alpha = scale; k = sh; h = sh2
  if (abs(sh) < 10^-7) {
    sh <- 10^(-12)
  }
  isGEV <- abs(sh2) < 10^-7
  isGLO <- (sh2 < -1+10^-7 & sh2 > -1-10^-7)
  if (isGEV) {
    Fd <- pgev(q = q, loc = loc, scale = scale, sh = sh, lower.tail = lower.tail)
  }
  if (isGLO) {
    Fd <- pglo(q = q,loc = loc, scale = scale, sh = sh, lower.tail = lower.tail)
  }
  if(!(isGEV | isGLO)){
     Fd <- (1 - h * (1-k * (q - xi)/alpha)^(1/k))^(1/h)
  }
  if(lower.tail) return(Fd)
  else return(1-Fd)
}

#' @name Kappa distribution
#' @family kappa distribution
#' @export
rkappa <- function(n, loc, scale, sh, sh2) {
  ## in HW notation 
  allTogether <- cbind(rep(loc, length.out=n), 
                       rep(scale, length.out=n),
                       rep(sh, length.out=n),
                       rep(sh2, length.out=n))
  ## add random percentiles
  allTogether <- cbind(runif(n), allTogether)
  apply(allTogether, 1, qkappa_int, lower.tail=FALSE)
  # apply(allTogether, 1, rkappa_int)
}

# rkappa_int <- function(x){
#   loc = x[1]; scale= x[2]; sh= x[3]; sh2 = x[4]
#   qkappa(p = runif(1), loc = loc, scale = scale, sh = sh, sh2=sh2)
# }


#' @name Kappa distribution
#' @family kappa distribution
#' @export
qkappa <- function(p, loc, scale, sh, sh2, lower.tail = TRUE) {
  ## in HW notation 
  allTogether <- cbind(p,loc,scale,sh,sh2)
  apply(allTogether,1,
        qkappa_int,lower.tail=lower.tail)
}

qkappa_int <- function(x, lower.tail = TRUE) {
  p = x[1]; loc = x[2]; scale= x[3]; sh= x[4]; sh2 = x[5]
  ## in HW notation 
  xi = loc;  alpha = scale; k = sh; h = sh2
  if (abs(sh) < 10^-7) {
    sh <- 10^(-12)
  }
  isGEV <- abs(sh2) < 10^-7
  isGLO <- (sh2 < -1+10^-7 & sh2 > -1-10^-7)

  if (isGEV) {
    qf <- qgev(p = p, loc = loc, scale = scale, sh = sh, lower.tail = lower.tail)
  }
  if (isGLO) {
    qf <- qglo(p = p, loc = loc, scale = scale, sh = sh, lower.tail = lower.tail)
  }
  if(!(isGEV | isGLO)){
    if(lower.tail) qf <- xi + (alpha/k) * (1 - ((1 -p^h)/h)^k) 
    if(!lower.tail) qf <- xi + (alpha/k) * (1 - ((1 -(1-p)^h)/h)^k) 
  }
qf
}

#' @name Kappa distribution
#' @family kappa distribution
#' @export
dkappa <- function(x, loc, scale, sh, sh2) {
  ## in HW notation 
  allTogether <- cbind(x,loc,scale,sh,sh2)
  apply(allTogether,1,
        dkappa_int)
}

dkappa_int <- function(x) {
  loc = x[2]; scale= x[3]; sh= x[4]; sh2 = x[5]; x = x[1]
  ## in HW notation 
  xi = loc;  alpha = scale; k = sh; h = sh2
  if (abs(sh) < 10^-7) {
    sh <- 10^(-12)
  }
  isGEV <- abs(sh2) < 10^-7
  isGLO <- (sh2 < -1+10^-7 & sh2 > -1-10^-7)
  
  if (isGEV) {
    f <- dgev(x = x, loc = loc, scale = scale, sh = sh)
  }
  if (isGLO) {
    f <- dglo(x = x, loc = loc, scale = scale, sh = sh)
  }
  if(!(isGEV | isGLO)){
    f <- (alpha^(-1)) * (1 - k * (x - xi)/alpha)^((1/k) - 1) * 
      (pkappa(q = x, loc = xi, scale = alpha, sh = k, sh2 = h))^(1 - h)
  }
  return(f)
}


#### Kappa distribution fitting ----- 


#' @title Maximum Likelihood Fitting for the Kappa distribution - CV model
#' @description Maximum-likelihood fitting for the Kappa distribution, 
#' including generalized linear modelling of each parameter.  The function differs from 
#' \code{kappad.fit} because it uses a different parametrisation of the distribution based on 
#' the \eqn{\tau}, the ratio of the scale parameter and the location parameter. This means that when regression models are applied for the location,
#' these also affect the scale. 
#' The function allows any parameter to be kept fixed and to not be estimated. 
#' @param xdat A numeric vector of data to be fitted
#' @param ydat A matrix of covariates for generalized linear modelling of the parameters (or NULL (the default) for stationary fitting). The number of rows should be the same as the length of xdat
#' @param mul  Numeric vectors of integers, giving the columns of ydat that contain covariates for generalized linear modelling of the location parameter (or NULL (the default) if the corresponding parameter is stationary)
#' @param taul As \code{mul} for the tau parameter
#' @param shl As \code{mul} for the shape parameter
#' @param sh2l As \code{mul} for the second shape parameter
#' @param mulink the link function for the location parameter - default to identity
#' @param taulink the link function for the tau parameter - default to identity
#' @param shlink the link function for the shape parameter - default to identity
#' @param sh2link the link function for the second shape parameter - default to identity
#' @param muinit initial values for the location parameter
#' @param tauinit initial values for the tau parameter
#' @param shinit initial values for the shape parameter
#' @param sh2init initial values for the second shape parameter
#' @param show Logical; if \code{TRUE} (the default), print details of the fit.
#' @param method The optimization method (see \code{optim} for details)
#' @param optimPars A string with other parameters to pass into \code{optim}. For example, depending on \code{method}, one could have "lower = 10, upper = 20"
#' @param maxit  The maximum number of iterations 
#' @param fixedPars a named list to fix any of the distribution parameter to a given value. When the named parameter is set to \code{NULL} its value is estimated.
#' @param ...   Other control parameters for the optimization. These are passed to components of the control argument of optim.
#' @return An object of the kappacv.fit class - with values which mirror the ones of the gev.fit class in \code{ismev}. 
#' 
#' In the output the \code{vals} matrix gives the location and scale values obtained as scale = \eqn{\tau} *  location.  
#' @seealso \link{kappad.fit}
#' @references 
#' Hosking, J.R.M. and Wallis, J.R., 2005. Regional frequency analysis: an approach based on L-moments. Cambridge university press. 
#' @export
#' @examples
#' set.seed(12)
#' x <- runif(120)
#' y <- rkappa(120,loc = 40+3*x,
#'    scale = 0.2*(40+3*x), sh = -0.2, sh2=-0.4)
#' fit1 <- kappacvd.fit(y, show=FALSE)
#' fit1
#' ## now add a regression model for the location
#' fit2 <- kappacvd.fit(y, ydat = cbind(x), mul=1, show=FALSE)
#' fit2
#' ## now a fit with a fixed shape parameter 
#' fitf2 <- kappacvd.fit(y, show=FALSE, fixedPars = list(sh2 = -0.4))
#' fitf2 ## only three parameters are estimated 
kappacvd.fit <- function(xdat, ydat = NULL, mul = NULL, taul = NULL, shl = NULL, sh2l = NULL,
                         mulink = identity, taulink = identity, shlink = identity, sh2link = identity,
                         muinit = NULL, tauinit = NULL, shinit = NULL, sh2init = NULL, 
                         show = TRUE, 
                         method = "Nelder-Mead", optimPars = NULL, maxit = 10000, 
                         fixedPars = list(mu = NULL, sig = NULL, sh = NULL, sh2 = NULL), ...) {
  z <- list()
  npmu <- ifelse(is.null(fixedPars[["mu"]]), length(mul) + 1, 0)
  npcv <- ifelse(is.null(fixedPars[["tau"]]), length(taul) + 1, 0)
  npsh<- ifelse(is.null(fixedPars[["sh"]]), length(shl) + 1, 0)
  npsh2<- ifelse(is.null(fixedPars[["sh2"]]), length(sh2l) + 1, 0)
  z$trans <- FALSE
  in2 <- sqrt(6 * var(xdat))/pi
  in1 <- mean(xdat) - 0.57722 * in2
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
  if (is.null(taul)) {
    taumat <- as.matrix(rep(1, length(xdat)))
    if (is.null(tauinit)) 
      tauinit <- in2/in1
  }
  else {
    z$trans <- TRUE
    taumat <- cbind(rep(1, length(xdat)), ydat[, taul])
    if (is.null(tauinit)) 
      tauinit <- c(in2/in1, rep(0, length(taul)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdat)))
    if (is.null(shinit)) 
      shinit <- -0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
    if (is.null(shinit)) 
      shinit <- c(-0.1, rep(0, length(shl)))
  }
  
  if (is.null(sh2l)) {
    sh2mat <- as.matrix(rep(1, length(xdat)))
    if (is.null(sh2init)) 
      sh2init <- -0.1
  }
  else {
    z$trans <- TRUE
    sh2mat <- cbind(rep(1, length(xdat)), ydat[, sh2l])
    if (is.null(shinit))
      sh2init <- c(-0.1, rep(0, length(sh2l)))
  }
  z$model <- list(mul, taul, shl,sh2l)
  z$link <- deparse(substitute(c(mulink, taulink, shlink, sh2link)))
  if(is.null(fixedPars[["mu"]])) muinit <- muinit
  if(!is.null(fixedPars[["mu"]])) muinit <- NULL
  if(is.null(fixedPars[["tau"]])) tauinit <- tauinit
  if(!is.null(fixedPars[["tau"]])) tauinit <-NULL
  if(is.null(fixedPars[["sh"]])) shinit <- shinit
  if(!is.null(fixedPars[["sh"]])) shinit <- NULL
  if(is.null(fixedPars[["sh2"]])) sh2init <- sh2init
  if(!is.null(fixedPars[["sh2"]])) sh2init <- NULL
  init <- c(muinit, tauinit, shinit, sh2init)
  # kappa.lik <- function(xdat, mu , sc, sh, sh2val) {
  kappa.lik <- function(a) {
    if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (a[1:npmu]))
    if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
    if(is.null(fixedPars[["tau"]])) cv <- taulink(taumat %*% (a[seq(npmu + 1, length = npcv)]))
    if(!is.null(fixedPars[["tau"]])) cv <-fixedPars[["tau"]]
    if(is.null(fixedPars[["sh"]])) sh <-  shlink(shmat %*% (a[seq(npmu + npcv + 1, length = npsh)]))
    if(!is.null(fixedPars[["sh"]])) sh <- fixedPars[["sh"]]
    if(is.null(fixedPars[["sh2"]])) sh2 <- sh2link(shmat %*% (a[seq(npmu + npcv + npsh + 1, length = npsh2)]))
    if(!is.null(fixedPars[["sh2"]])) sh2 <- fixedPars[["sh2"]]
    sc <- cv*mu
    y <- (xdat - mu)/sc
    y <- 1 - sh * y
    if (any(y <= 0) || any(sc <= 0))
      return(10^6)
    sum(log(sc)) + 
      sum((1-(1/sh))*log(y)) + 
      sum((sh2-1) * log(pkappa(q = xdat, loc = mu, scale = sc, sh=sh, sh2=sh2)))
    # 
    # f <- (alpha^(-1)) * (1 - k * (x - xi)/alpha)^((1/k) - 1) * 
    #   (pkappa(q = x, loc = xi, scale = alpha, sh = k, sh2 = h))^(1 - h)
    # 
  }
  x <- eval(parse(text = paste0(
    "optim(init, kappa.lik, ",
    optimPars,
    ",hessian = TRUE, method = method, control = list(maxit = maxit, ...))")))
  # x <- optim(init, kappa.lik, hessian = TRUE, method = method, 
  #            control = list(maxit = maxit, ...))
  z$conv <- x$convergence
  if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (x$par[1:npmu]))
  if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
  if(is.null(fixedPars[["tau"]])) cv <- taulink(taumat %*% (x$par[seq(npmu + 1, length = npcv)]))
  if(!is.null(fixedPars[["tau"]])) cv <-fixedPars[["tau"]]
  if(is.null(fixedPars[["sh"]])) sh <-  shlink(shmat %*% (x$par[seq(npmu + npcv + 1, length = npsh)]))
  if(!is.null(fixedPars[["sh"]])) sh <- fixedPars[["sh"]]
  if(is.null(fixedPars[["sh2"]])) sh2 <- sh2link(shmat %*% (x$par[seq(npmu + npcv + npsh + 1, length = npsh2)]))
  if(!is.null(fixedPars[["sh2"]])) sh2 <- fixedPars[["sh2"]]
  z$nllh <- x$value
  z$data <- xdat
  # if (z$trans) {
  #   z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
  # }
  z$mle <- x$par
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  vals <- cbind(mu, cv*mu, sh, sh2)
  colnames(vals) <- c("location","scale","shape1","shape2")
  z$vals <- vals 
  if (show) {
    if (z$trans) 
      print(z[c(2, 3, 4)])
    else print(z[4])
    if (!z$conv) 
      print(z[c(5, 7, 9)])
  }
  class(z) <- "kappacv.fit"
  invisible(z)
}



#' @title Maximum Likelihood Fitting for the Kappa distribution
#' @description Maximum-likelihood fitting for the Kappa distribution, 
#' including generalized linear modelling of each parameter. 
#' This function has the same structure as the \code{gevd.fit} and is inspired by \code{ismev::gev.fit}. 
#' The function allows any parameter to be kept fixed and to not be estimated. 
#' @param xdat A numeric vector of data to be fitted
#' @param ydat A matrix of covariates for generalized linear modelling of the parameters (or NULL (the default) for stationary fitting). The number of rows should be the same as the length of xdat
#' @param mul  Numeric vectors of integers, giving the columns of ydat that contain covariates for generalized linear modelling of the location parameter (or NULL (the default) if the corresponding parameter is stationary)
#' @param sigl As \code{mul} for the scale parameter
#' @param shl As \code{mul} for the shape parameter
#' @param sh2l As \code{mul} for the second shape parameter
#' @param mulink the link function for the location parameter - default to identity
#' @param siglink the link function for the scale parameter - default to identity
#' @param shlink the link function for the shape parameter - default to identity
#' @param sh2link the link function for the second shape parameter - default to identity
#' @param muinit initial values for the location parameter
#' @param siginit initial values for the scale parameter
#' @param shinit initial values for the shape parameter
#' @param sh2init initial values for the second shape parameter
#' @param show Logical; if \code{TRUE} (the default), print details of the fit.
#' @param method The optimization method (see \code{optim} for details).
#' @param optimPars A string with other parameters to pass into \code{optim}. For example, depending on \code{method}, one could have "lower = 10, upper = 20"
#' @param maxit  The maximum number of iterations.
#' @param fixedPars a named list to fix any of the distribution parameter to a given value. When the named parameter is set to \code{NULL} its value is estimated.
#' @param ...   Other control parameters for the optimization. These are passed to components of the control argument of optim.
#' @return An object of the kappa.fit class - with values which mirror the ones of the gev.fit class in \code{ismev}
#' @details The distribution is discussed in the Hosking and Wallis book and can be seen as a generelisation of several distributions used in extreme value modelling. 
#' Depending on the sh2 parameter value the Kappa distribution reduces to others commonly used distributions For example: 
#' when sh2 is equal to -1 the distribution reduces to a GLO, when sh2 is equal to 0 the distribution reduces to a GEV, 
#' when sh2 is equal to +1 the distribution reduces to a Generalised Pareto distribution (GPA). 
#' @seealso \link{dkappa}
#' @export
#' @examples
#' set.seed(12)
#' x <- runif(500)
#' y <- rkappa(500,loc = 40+4*x,scale = 6, sh = 0.2, sh2=-0.4)
#' fit1 <- kappad.fit(y, show=FALSE)
#' fit1
#' ## now add a regression model for the location
#' fit2 <- kappad.fit(y, ydat = cbind(x), mul=1, show=FALSE)
#' fit2
#' ## now a fit with a fixed shape parameter 
#' fitf <- kappad.fit(y, show=FALSE, fixedPars = list(sh = 0.2))
#' fitf ## only three parameters are estimated 
#' ## now fix the second shape parameter
#' fitf2 <- kappad.fit(y, show=FALSE, fixedPars = list(sh2 = -0.4))
#' fitf2 ## only three parameters are estimated 
kappad.fit <- function(xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, sh2l = NULL,
                       mulink = identity, siglink = identity, shlink = identity, sh2link = identity,
                       muinit = NULL, siginit = NULL, shinit = NULL, sh2init = NULL, 
                       show = TRUE, method = "Nelder-Mead", optimPars=NULL, maxit = 10000, 
                       fixedPars = list(mu = NULL, sig = NULL, sh = NULL, sh2 = NULL), ...) {
  z <- list()
  npmu <- ifelse(is.null(fixedPars[["mu"]]), length(mul) + 1, 0)
  npsc <- ifelse(is.null(fixedPars[["sig"]]), length(sigl) + 1, 0)
  npsh<- ifelse(is.null(fixedPars[["sh"]]), length(shl) + 1, 0)
  npsh2<- ifelse(is.null(fixedPars[["sh2"]]), length(sh2l) + 1, 0)
  z$trans <- FALSE
  in2 <- sqrt(6 * var(xdat))/pi
  in1 <- mean(xdat) - 0.57722 * in2
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
      shinit <- -0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
    if (is.null(shinit)) 
      shinit <- c(-0.1, rep(0, length(shl)))
  }
  
  if (is.null(sh2l)) {
    sh2mat <- as.matrix(rep(1, length(xdat)))
    if (is.null(sh2init)) 
      sh2init <- -0.1
  }
  else {
    z$trans <- TRUE
    sh2mat <- cbind(rep(1, length(xdat)), ydat[, sh2l])
    if (is.null(shinit))
      sh2init <- c(-0.1, rep(0, length(sh2l)))
  }
  z$model <- list(mul, sigl, shl,sh2l)
  z$link <- deparse(substitute(c(mulink, siglink, shlink, sh2link)))
  if(is.null(fixedPars[["mu"]])) muinit <- muinit
  if(!is.null(fixedPars[["mu"]])) muinit <- NULL
  if(is.null(fixedPars[["sig"]])) siginit <- siginit
  if(!is.null(fixedPars[["sig"]])) siginit <-NULL
  if(is.null(fixedPars[["sh"]])) shinit <- shinit
  if(!is.null(fixedPars[["sh"]])) shinit <- NULL
  if(is.null(fixedPars[["sh2"]])) sh2init <- sh2init
  if(!is.null(fixedPars[["sh2"]])) sh2init <- NULL
  init <- c(muinit, siginit, shinit, sh2init)
  # kappa.lik <- function(xdat, mu , sc, sh, sh2val) {
  kappa.lik <- function(a) {
    if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (a[1:npmu]))
    if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
    if(is.null(fixedPars[["sig"]])) sc <-  siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    if(!is.null(fixedPars[["sig"]])) sc <-fixedPars[["sig"]]
    if(is.null(fixedPars[["sh"]])) sh <-  shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    if(!is.null(fixedPars[["sh"]])) sh <- fixedPars[["sh"]]
    if(is.null(fixedPars[["sh2"]])) sh2 <- sh2link(shmat %*% (a[seq(npmu + npsc + npsh + 1, length = npsh2)]))
    if(!is.null(fixedPars[["sh2"]])) sh2 <- fixedPars[["sh2"]]
    y <- (xdat - mu)/sc
    y <- 1 - sh * y
    if (any(y <= 0) || any(sc <= 0))
      return(10^6)
    sum(log(sc)) + 
      sum((1-(1/sh))*log(y)) + 
      sum((sh2-1) * log(pkappa(q = xdat, loc = mu, scale = sc, sh=sh, sh2=sh2)))
    # 
    # f <- (alpha^(-1)) * (1 - k * (x - xi)/alpha)^((1/k) - 1) * 
    #   (pkappa(q = x, loc = xi, scale = alpha, sh = k, sh2 = h))^(1 - h)
    # 
  }
  x <- eval(parse(text = paste0(
    "optim(init, kappa.lik, ",
    optimPars,
    ",hessian = TRUE, method = method, control = list(maxit = maxit, ...))")))
  # x <- optim(init, kappa.lik, hessian = TRUE, method = method, 
  #            control = list(maxit = maxit, ...))
  z$conv <- x$convergence
  if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (x$par[1:npmu]))
  if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
  if(is.null(fixedPars[["sig"]])) sc <-  siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
  if(!is.null(fixedPars[["sig"]])) sc <-fixedPars[["sig"]]
  if(is.null(fixedPars[["sh"]])) sh <-  shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
  if(!is.null(fixedPars[["sh"]])) sh <- fixedPars[["sh"]]
  if(is.null(fixedPars[["sh2"]])) sh2 <- sh2link(shmat %*% (x$par[seq(npmu + npsc + npsh + 1, length = npsh2)]))
  if(!is.null(fixedPars[["sh2"]])) sh2 <- fixedPars[["sh2"]]
  z$nllh <- x$value
  z$data <- xdat
  # if (z$trans) {
  #   z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
  # }
  z$mle <- x$par
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  vals <- cbind(mu, sc, sh, sh2)
  names(vals) <- c("location","scale","shape1","shape2")
  z$vals <- vals 
  if (show) {
    if (z$trans) 
      print(z[c(2, 3, 4)])
    else print(z[4])
    if (!z$conv) 
      print(z[c(5, 7, 9)])
  }
  class(z) <- "kappa.fit"
  invisible(z)
}



# 
# kappagivensh2.fit <- function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, 
#           mulink = identity, siglink = identity, shlink = identity, 
#           muinit = NULL, siginit = NULL, shinit = NULL, sh2val = -0.5, 
#           show = TRUE, method = "Nelder-Mead", maxit = 10000, ...) {
#   z <- list()
#   npmu <- length(mul) + 1
#   npsc <- length(sigl) + 1
#   npsh <- length(shl) + 1
#   z$trans <- FALSE
#   in2 <- sqrt(6 * var(xdat))/pi
#   in1 <- mean(xdat) - 0.57722 * in2
#   if (is.null(mul)) {
#     mumat <- as.matrix(rep(1, length(xdat)))
#     if (is.null(muinit)) 
#       muinit <- in1
#   }
#   else {
#     z$trans <- TRUE
#     mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
#     if (is.null(muinit)) 
#       muinit <- c(in1, rep(0, length(mul)))
#   }
#   if (is.null(sigl)) {
#     sigmat <- as.matrix(rep(1, length(xdat)))
#     if (is.null(siginit)) 
#       siginit <- in2
#   }
#   else {
#     z$trans <- TRUE
#     sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
#     if (is.null(siginit)) 
#       siginit <- c(in2, rep(0, length(sigl)))
#   }
#   if (is.null(shl)) {
#     shmat <- as.matrix(rep(1, length(xdat)))
#     if (is.null(shinit)) 
#       shinit <- -0.1
#   }
#   else {
#     z$trans <- TRUE
#     shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
#     if (is.null(shinit)) 
#       shinit <- c(-0.1, rep(0, length(shl)))
#   }
#   z$model <- list(mul, sigl, shl)
#   z$link <- deparse(substitute(c(mulink, siglink, shlink)))
#   init <- c(muinit, siginit, shinit)
#   # kappa.lik <- function(xdat, mu , sc, sh, sh2val) {
# 
#   kappa.lik <- function(a) {
#     mu <- mulink(mumat %*% (a[1:npmu]))
#     sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
#     sh <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
#     y <- (xdat - mu)/sc
#     y <- 1 - sh * y
#  
#     if (any(y <= 0) || any(sc <= 0)) 
#       return(10^6)
#     sum(log(sc)) + 
#       sum((1-(1/sh))*log(y)) + 
#       sum((sh2val-1) * log(pkappa(q = xdat, loc = mu, scale = sc, sh=sh, sh2=sh2val)))
#     # 
#     # f <- (alpha^(-1)) * (1 - k * (x - xi)/alpha)^((1/k) - 1) * 
#     #   (pkappa(q = x, loc = xi, scale = alpha, sh = k, sh2 = h))^(1 - h)
#     # 
#   }
#   x <- optim(init, kappa.lik, hessian = TRUE, method = method, 
#              control = list(maxit = maxit, ...))
#   z$conv <- x$convergence
#   mu <- mulink(mumat %*% (x$par[1:npmu]))
#   sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
#   xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
#   z$nllh <- x$value
#   z$data <- xdat
#   # if (z$trans) {
#   #   z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
#   # }
#   z$mle <- x$par
#   z$cov <- solve(x$hessian)
#   z$se <- sqrt(diag(z$cov))
#   z$vals <- cbind(mu, sc, xi)
#   if (show) {
#     if (z$trans) 
#       print(z[c(2, 3, 4)])
#     else print(z[4])
#     if (!z$conv) 
#       print(z[c(5, 7, 9)])
#   }
#   class(z) <- "kappa.fit"
#   invisible(z)
# }
# 
# y <- rkappa(500,30,6,-0.3,-1)
# plot(y)
# kappagivensh2.fit(xdat = y, sh2val = -1)
# ilaprosUtils::glo.fit(y)
# 

# kappagivensh2_cv.fit <- function(xdat, ydat = NULL, mul = NULL, taul = NULL, shl = NULL, 
#                         mulink = identity, taulink = identity, shlink = identity, 
#                         muinit = NULL, tauinit = NULL, shinit = NULL, 
#                         hinit = NULL, show = TRUE, 
#                         method = "Nelder-Mead", maxit = 10000, 
#                         fixedPars = list(mu = NULL, tau = NULL, sh = NULL, h = NULL),...) {
#   z <- list()
#   npmu <- length(mul) + 1
#   npcv <- length(taul) + 1
#   npsh <- length(shl) + 1
#   z$trans <- FALSE
#   in2 <- sqrt(6 * var(xdat))/pi
#   in1 <- mean(xdat) - 0.57722 * in2
#   ins <- in2/in1
#   if (is.null(mul)) {
#     mumat <- as.matrix(rep(1, length(xdat)))
#     if (is.null(muinit)) 
#       muinit <- in1
#   }
#   else {
#     z$trans <- TRUE
#     mumat <- cbind(rep(1, length(xdat)), ydat[, mul])
#     if (is.null(muinit)) 
#       muinit <- c(in1, rep(0, length(mul)))
#   }
#   if (is.null(taul)) {
#     taumat <- as.matrix(rep(1, length(xdat)))
#     if (is.null(tauinit)) 
#       tauinit <- in2
#   }
#   else {
#     z$trans <- TRUE
#     taumat <- cbind(rep(1, length(xdat)), ydat[, taul])
#     if (is.null(tauinit)) 
#       tauinit <- c(in2, rep(0, length(taul)))
#   }
#   if (is.null(shl)) {
#     shmat <- as.matrix(rep(1, length(xdat)))
#     if (is.null(shinit)) 
#       shinit <- 0.1
#   }
#   else {
#     z$trans <- TRUE
#     shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
#     if (is.null(shinit)) 
#       shinit <- c(0.1, rep(0, length(shl)))
#   }
#   if (is.null(hinit)) hinit <- 0.1
#   z$model <- list(mul, taul, shl)
#   z$link <- deparse(substitute(c(mulink, taulink, shlink)))
#   if(is.null(fixedPars[["mu"]])) muinit <- muinit
#   if(!is.null(fixedPars[["mu"]])) muinit <- NULL
#   if(is.null(fixedPars[["tau"]])) tauinit <- tauinit
#   if(!is.null(fixedPars[["tau"]])) tauinit <-NULL
#   if(is.null(fixedPars[["sh"]])) shinit <- shinit
#   if(!is.null(fixedPars[["sh"]])) shinit <- NULL
#   if(is.null(fixedPars[["h"]])) hinit <- hinit
#   if(!is.null(fixedPars[["h"]])) hinit <- NULL
#   init <- c(muinit, tauinit, shinit, hinit)
#   print(length(init))
#   kappa.lik <- function(a) {
#     if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (a[1:npmu]))
#     if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
#     if(is.null(fixedPars[["tau"]])) cv <- taulink(taumat %*% (a[seq(npmu + 1, length = npcv)]))
#     if(!is.null(fixedPars[["tau"]])) cv <-fixedPars[["tau"]]
#     if(is.null(fixedPars[["sh"]])) xi <- shlink(shmat %*% (a[seq(npmu + npcv + 1, length = npsh)]))
#     if(!is.null(fixedPars[["sh"]])) xi <- fixedPars[["sh"]]
#     if(is.null(fixedPars[["h"]])) h <- rep(a[npmu + npcv + npsh + 1], length(xdat))
#     if(!is.null(fixedPars[["h"]])) h <- rep(fixedPars[["h"]], length(xdat))
#     sc <- cv*mu
#     ### lots of conditions 
#     if(any(sc <= 0)) return(10^6)
#     if(any(h <= -1)) return(10^6)
#     if(any(xi <= -1)) return(10^6)
#     if(h > 0 & ((xdat - mu)/((1-h^(-xi))*sc)) < 0) return(10^6) 
#     if(h > 0 & any(xi > 0) & ((xdat - mu)/sc) > 0) return(10^6) 
#     if(h > 0 & abs(k) < 10^-12 & ((xdat - mu)/(sc*log(h)) < 0)) return(10^6) 
#    #### negative cases to be done 
#     
#     # y <- 1 + xi * y
#     print(c(mu[1],cv[1],xi[1],h[1],y))
#     else -sum(apply(cbind(xdat, mu, sc, xi, h), 1, 
#               function(x) log(dkappa(x = x[1], xi = x[2], alpha = x[3], k = x[4], h = x[5]))))
#   }
#   x <- optim(init, kappa.lik, hessian = TRUE, method = method, 
#              control = list(maxit = maxit, ...))
#   z$conv <- x$convergence
#   if(is.null(fixedPars[["mu"]])) mu <- mulink(mumat %*% (x$par[1:npmu])) 
#   if(!is.null(fixedPars[["mu"]])) mu <- fixedPars[["mu"]]
#   if(is.null(fixedPars[["tau"]])) cv <- taulink(taumat %*% (x$par[seq(npmu + 1, length = npcv)]))
#   if(!is.null(fixedPars[["tau"]])) cv <-fixedPars[["tau"]]
#   if(is.null(fixedPars[["sh"]])) xi <- shlink(shmat %*% (x$par[seq(npmu + npcv + 1, length = npsh)]))
#   if(!is.null(fixedPars[["sh"]])) xi <- fixedPars[["sh"]]
#   if(is.null(fixedPars[["h"]])) h <- x$par[length(x$par)]
#   if(!is.null(fixedPars[["h"]])) h <- fixedPars[["h"]]
#   sc <- mu * cv
#   z$nllh <- x$value
#   z$data <- xdat
#   # if (z$trans) {
#   #   z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
#   # }
#   z$mle <- x$par
#   z$cov <- solve(x$hessian)
#   z$se <- sqrt(diag(z$cov))
#   z$vals <- cbind(mu, cv, xi, h)
#   if (show) {
#     if (z$trans) 
#       print(z[c(2, 3, 4)])
#     else print(z[4])
#     if (!z$conv) 
#       print(z[c(5, 7, 9)])
#   }
#   class(z) <- "kappa.fit"
#   invisible(z)
# }
# 



#' Nicer print of kappa.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param x a fitted object of the class kappa.fit.
#' @param ... further arguments passed to \code{print}.
#' @seealso \code{\link{kappad.fit}}
#' @export
print.kappa.fit<-function(x, ...){
  zz<-list(mle=x$mle,se=x$se,conv=x$conv,nllh=x$nllh)
  print(zz, ...)
}

#' Nicer print of kappa.fit results
#'
#' This function prints the MLE, se, convergence info and negative log-likelihood value.
#' @param x a fitted object of the class kappa.fit.
#' @param ... further arguments passed to \code{print}.
#' @seealso \code{\link{kappacvd.fit}}
#' @export
print.kappacv.fit<-function(x, ...){
  zz<-list(mle=x$mle,se=x$se,conv=x$conv,nllh=x$nllh)
  print(zz, ...)
}

