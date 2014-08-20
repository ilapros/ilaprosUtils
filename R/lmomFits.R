
PPmatG<-function(mat,cols=NULL,fit.lines="gev",ff=seq(0.005,0.995,by=0.0025),yll=NULL,timeGrid=TRUE,...){
  if(is.null(cols) | length(cols)<ncol(mat)) cols<-seq(1,ncol(mat))
  if(is.null(yll)) yll<-c(0,1.2*max(mat[,1:ncol(mat)],na.rm=TRUE))
  prevPlot<-FALSE
  mat.pars<-matrix(NA,ncol=switch(fit.lines,
                                    "gev"=3,
                                    "gamma"=2,
                                    "none"=4),nrow=ncol(mat))
  plot(range(-log(-log(ff))),yll,xlab=" ",ylab=" ",type="n",bty="l",...)  
  for(j in 1:ncol(mat)){
    x<-mat[,j];x<-x[!is.na(x)]
    prevPlot<-((all(is.na(x))))
    if(!all(is.na(x))) {
        axis(1)
        femp<-(1:length(x)-0.44)/(length(x) + 1 - 0.88)
#         print(femp)
        points(-log(-log(femp)),sort(x),col=cols[j],...)
        if(length(x)>3){
          lines(-log(-log(ff)),switch(fit.lines,
                "gev"=quagev(ff,pelgev(samlmu(x))),
                "gamma"=quagam(ff,pelgam(samlmu(x))),
                "glo"=quaglo(ff,pelglo(samlmu(x))),
                "none"=rep(-500,l=length(ff))),col=cols[j])
          mat.pars[j,]<-switch(fit.lines,
                  "gev"=pelgev(samlmu(x)),
                  "gamma"=pelgam(samlmu(x)),
                  "glo"=pelglo(samlmu(x)),
                  "none"=samlmu(x))
        }
    }    
 if(timeGrid){
       abline(v=-log(-log(1-1/c(1.2,2,5,10,50,100,200))),lty=2,col=8)
       axis(3,col=0,at=-log(-log(1-1/c(1.2,2,5,10,50,100,200))),labels=paste(c(1.2,2,5,10,50,100,200),"yrs"),cex.axis=.75,line=-1,cex=0.8)
     } 
  }
  invisible(mat.pars)
}




PPlinesG<-function(cols=NULL,fit.lines="gev",ff=seq(0.005,0.995,by=0.0025),pars=NULL,...){
  if(is.null(cols) | length(cols)<nrow(pars)) cols<-seq(1,ncol(pars))
  for(j in 1:nrow(pars)){
      lines(-log(-log(ff)),switch(fit.lines,
            "gev"=quagev(ff,as.numeric(pars[j,])),
            "gamma"=quagam(ff,as.numeric(pars[j,])),
            "glo"=quaglo(ff,as.numeric(pars[j,]))),col=cols[j], ...)
  }
}





ll<-seq(0.005,0.995,by=0.0025)
retP<-c(1.2,1.5,2,3.5,5,10,20,50,100,200)

## GB and NI
ukcoast<-read.table("P:\\HRRD\\NEC05031 EA Small Catchments Phase 2\\analysis\\ukcoast.txt",header=TRUE)
ukcoast<-ukcoast/1000

gbcoast <- read.table("P:\\HRRD\\NEC05031 EA Small Catchments Phase 2\\analysis\\gbcoast.txt",header=TRUE)


par.set.size <- function(x, y){
  limits <- c(range(x), range(y))
  foo.l <<- limits
  fin <- par()$fin
  #scale the plots to maximum size but with x and y scales equal
  xrange <- limits[2] - limits[1]
  yrange <- limits[4] - limits[3]
  xscale <- xrange/(fin[1] - par()$mai[2] - par()$mai[4])
  yscale <- yrange/(fin[2] - par()$mai[1] - par()$mai[3])
  samescale <- max(xscale, yscale)
  par(pin = (1/samescale) * c(xrange, yrange), usr = limits,mai=c(.60000,0.75,.4,0.560000))
  invisible()
}

do.uk<-function(mycol=1,xlim=c(0,670),ylim=c(0,1050)){
  par.set.size(xlim,ylim)
  plot(xlim, ylim, type = "n", 
       xlab = "Easting (km)", ylab = "Northing (km)",
       mgp=c(1.7,.6,0), axes = T, bty="l")
  #          mgp= c(0.8, 0.2, 0), axes = T)
  lines(ukcoast$east, ukcoast$north, lwd = 0.1, col = mycol)
}

