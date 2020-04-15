# 
# 
# load_all("../ilaprosUtils/"); document("../ilaprosUtils");
# build("../ilaprosUtils");install("../ilaprosUtils");library(ilaprosUtils)
# ?PPmatG



#' @title Plotting positions of a matrix of data
#' @description This can be useful when building DDF models or when comparing seasonal maxima of the same record. 
#' It takes each columns of a matrix and plots them as Gringorten plotting positions. 
#' The x-axis is shown using the Gumbel variate \code{-log(-log(f))}, with f the non exceedance probability
#' If a distribution is given to fit.lines, it adds an estimated return level curve based on L-moments estimate
#' @param mat data matrix, each column will be analysed separately.
#' @param cols colour in which the data points of each column (and estimated line) should be displayed
#' @param fit.lines a character string indicating the distribution to be used to draw the fitted line. 
#' Current options are "gev" (default), "glo", "gamma", "gpa" and "none" which results in no lines
#' @param ff the non-exceedance probabilities for which the lines are estimated. 
#' It also affects the width of the x-axis. Default is \code{seq(0.005,0.995,by=0.0025)}.
#' @param yll ylimits (optional)
#' @param timeGrid logical. If \code{TRUE} a time grid with return periods expressed in years
#' corresponding to values in \code{vgrid} will be displayed.
#' @param vgrid return periods expressed in years shown as a grid
#' @param ... other parameters passed into \code{plot} 
#' @return A plot with the plotting positions of each column and a matrix of the estimated parameters  
#' @examples
#' x <- matrix(c(rgev(20, 10, 3,   -0.2), 
#'               rgev(20, 15, 4.5, -0.2), 
#'               rgev(20, 40, 6.5, -0.2)),byrow = FALSE, ncol=3)
#' PPmatG(x, cols = c(4,5,2), fit.lines = "gev")
#' @export
#' @import lmom
PPmatG <- function(mat,cols=NULL,fit.lines="gev",
                   ff=seq(0.005,0.995,by=0.0025),
                   yll=NULL,
                   timeGrid=TRUE, 
                   vgrid = c(1.2,2,5,10,50,100,200),
                   ...){
  if(is.null(cols) | length(cols)<ncol(mat)) cols<-seq(1,ncol(mat))
  if(is.null(yll)) yll<-c(0,1.2*max(mat[,1:ncol(mat)],na.rm=TRUE))
  prevPlot<-FALSE
  mat.pars<-matrix(NA,ncol=switch(fit.lines,
                                    "gev"=3,
                                    "glo"=3,
                                    "gpa"=3,
                                    "gamma"=2,
                                    "none"=4),nrow=ncol(mat))
  plot(range(-log(-log(ff))),yll,xlab=" ",ylab=" ",type="n",bty="l",...)  
  for(j in 1:ncol(mat)){
    x<-mat[,j];x<-x[!is.na(x)]
    prevPlot<-((all(is.na(x))))
    if(!all(is.na(x))) {
        axis(1)
        femp<-(1:length(x)-0.44)/(length(x) + 1 - 0.88)
        points(-log(-log(femp)),sort(x),col=cols[j],...)
        if(length(x)>3){
          lines(-log(-log(ff)),switch(fit.lines,
                "gev"=quagev(ff,pelgev(samlmu(x))),
                "gamma"=quagam(ff,pelgam(samlmu(x))),
                "gpa"=quagpa(samlmu(x)),
                "glo"=quaglo(ff,pelglo(samlmu(x))),
                "none"=rep(-500,l=length(ff))),col=cols[j])
          mat.pars[j,]<-switch(fit.lines,
                  "gev"=pelgev(samlmu(x)),
                  "gamma"=pelgam(samlmu(x)),
                  "glo"=pelglo(samlmu(x)),
                  "gpa"=pelgpa(samlmu(x)),
                  "none"=samlmu(x))
        }
    }    
 if(timeGrid){
       abline(v=-log(-log(1-1/vgrid)),lty=2,col=8)
       axis(3,col=0,at=-log(-log(1-1/vgrid)),
            labels=paste(vgrid,"yrs"),
            cex.axis=.75,line=-1,cex=0.8)
     } 
  }
  invisible(mat.pars)
}



#' @title Add return level curves based on a set of (estimated) parameters.
#' @description This can be useful in combination with PPmatG. 
#' The x-axis is shown using the Gumbel variate \code{-log(-log(f))}, with f the non exceedance probability. 
#' If a distribution is given to fit.lines, it adds an estimated return level curve based on L-moments estimate.
#' @param mat data matrix, each column will be analysed seperately.
#' @param pars matrix of parameter estimates. 
#' Each column should correspond to a parameter, each line to an observation (station).
#' @param cols colour in which the lines points of each row should be displayed
#' @param fit.lines a character string indicating the distribution to be used to draw the fitted line. 
#' Current options are "gev" (default), "glo", "gpa", "gamma"
#' @param ff the non-exceedance probabilities for which the lines are estimated
#' It also affects the width of the x-axis. Defaul is \code{seq(0.005,0.995,by=0.0025)}
#' @param ... additional graphical parameters
#' @return A set of estimated return curves added to the current plot
#' @examples
#' x <- matrix(c(rgev(20, 10, 3,   -0.2), 
#'               rgev(20, 15, 4.5, -0.2), 
#'               rgev(20, 40, 6.5, -0.2)), byrow = FALSE, ncol=3)
#' library(lmom)
#' PPmatG(x, cols = c(4,5,2), fit.lines = "gev")
#' pars <- t(apply(x,2,function(x) pelglo(samlmu(x))))
#' PPlinesG(pars = pars, cols= c(4,5,2), lty=2)
#' PPlinesG(mat  = x, cols= c(4,5,2), lty=3, fit.lines = "gamma")
#' @export
PPlinesG <- function(mat = NULL, pars = NULL, cols = NULL, fit.lines = "gev", ff = seq(0.005, 0.995, by = 0.0025), ...) {
  if(is.null(pars) & is.null(mat)) stop("Either a parameter or data matrix needed")
  if(!is.null(mat)) pars <- t(switch(fit.lines,gev = apply(mat, 2, function(x) pelgev(samlmu(x))), 
                                  gamma = apply(mat, 2, function(x) pelgam(samlmu(x))), 
                                  glo = apply(mat, 2, function(x) pelglo(samlmu(x))), 
                                  gpa = apply(mat, 2, function(x) pelgpa(samlmu(x)))))
  if (is.null(cols) | length(cols) < nrow(pars)) 
    cols <- seq(1, ncol(pars))
  for (j in 1:nrow(pars)) {
    lines(-log(-log(ff)), switch(fit.lines, 
                                 gev = quagev(ff, as.numeric(pars[j, ])), 
                                 gamma = quagam(ff, as.numeric(pars[j, ])), 
                                 glo = quaglo(ff, as.numeric(pars[j, ])),
                                 gpa = quagpa(ff, as.numeric(pars[j, ]))), col = cols[j], ...)
  }
  invisible(pars)
}




