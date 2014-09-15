# Use of ilaprosUtils
*Ilaria Prosdocimi* (ilapro@ceh.ac.uk)


This document aims at showing the use of some of the functions contained in the `ilaprosUtils` package. I will use the data obtained via the package `rnrfa`, which retrieves data from the National River Flow Archive. Indeed, `ilaprosUtils` is mostly a combination of functions I developed (and continue to develop) in my years at the Centre for Ecology and Hydrology and can be useful when working with hydrological data, in particular high end extremes. 

This document shows how to use some functions in the `ilaprosUtils` package to perform an analysis on the annual and seasonal maxima of a daily river flow as gauged at one of the stations of the National River Flow Archive. 

First a number of libraries are needed for everything to work:


```r
suppressPackageStartupMessages(library(zoo))
## this actually only needs to be installed 
## only needed to extract the index later on
```


```r
suppressPackageStartupMessages(library(devtools))
library(rnrfa)
# install_github("ilapros/ilaprosUtils")
library(ilaprosUtils)
library(lmom)
library(ismev)
```

We download via `rnrfa` average daily flow data for the Conder at Galgate, as measured at station 72014. The result is a pretty complex object, we are only interested in the daily river flow stored in `wmlTS`:

```r
st72014 <- SearchNRFA(72014)
```

```
## http://www.ceh.ac.uk/nrfa/xml/waterml2?db=nrfa_public&stn=72014&dt=gdf
```

```
## Warning: some methods for "zoo" objects do not work if the index entries
## in 'order.by' are not unique
```

```r
str(st72014)
```

```
## List of 2
##  $ wmlInfo:'data.frame':	1 obs. of  6 variables:
##   ..$ stationName      :List of 1
##   .. ..$ : chr "Conder at Galgate"
##   ..$ Latitude         :List of 1
##   .. ..$ : num 54
##   ..$ Longitude        :List of 1
##   .. ..$ : num -2.79
##   ..$ typeOfMeasurement:List of 1
##   .. ..$ : chr "Flow, m3/s, Mean, Day"
##   ..$ timeZone         :List of 1
##   .. ..$ :List of 2
##   .. .. ..$ zoneOffset      : chr "+00:00"
##   .. .. ..$ zoneAbbreviation: chr "GMT"
##   ..$ remarks          :List of 1
##   .. ..$ : chr "\nFlat V Crump profile weir in confined concrete channel to 1.775m, with concrete wall above as a flood barrier. Weir operates "| __truncated__
##  $ wmlTS  :'zoo' series from 1976-02-01 to NA
##   Data: num [1:13168] 0.247 0.246 0.241 0.24 0.228 0.311 0.395 0.311 0.321 0.292 ...
##   Index:  POSIXlt[1:13168], format: "1976-02-01" "1976-02-02" ...
```

```r
st72014 <- st72014[["wmlTS"]]
head(st72014)
```

```
## 1976-02-01 1976-02-02 1976-02-03 1976-02-04 1976-02-05 1976-02-06 
##      0.247      0.246      0.241      0.240      0.228      0.311
```

We simplify the structure of the data and just keep the date and flow value recorded at each day


```r
class(st72014)
```

```
## [1] "zoo"
```

```r
st72014 <- data.frame(Date = zoo::index(st72014), Flow = as.numeric(st72014))
plot(st72014,type="l")
```

![plot of chunk unnamed-chunk-4](./README_files/figure-html/unnamed-chunk-4.png) 

We see that some values are missing, in particular, for some yers there seem to be no record at all. 

The water year in the UK is commonly taken to start in October. Summer is taken to be the period running between March and October. We identify the Water year with the `date2wy` function and the summer month with the `date2summer` function. 


```r
st72014$WaterYear <- date2wy(st72014$Date)
## prints the month of the first and last day of each WaterYear
head(tapply(st72014$Date, factor(st72014$WaterYear), function(x) date2month(x[c(1,length(x))])))
```

```
## $`1975`
## [1] 2 9
## 
## $`1976`
## [1] 10  9
## 
## $`1977`
## [1] 10  9
## 
## $`1978`
## [1] 10  9
## 
## $`1979`
## [1] 10  9
## 
## $`1980`
## [1] 10  9
```

```r
st72014$Summer <- date2summer(st72014$Date)

plot(st72014$Date, st72014$Flow, type="l")
lines(st72014$Date[st72014$Summer == 1], st72014$Flow[st72014$Summer == 1], col = 2)
```

![plot of chunk unnamed-chunk-5](./README_files/figure-html/unnamed-chunk-5.png) 


Years in which less than 80% of the flow data are recorded are to be considered incomplete. If the interest lies in some extremal part of the distribution, one should be careful when keeping records with missing data in the dataset: it could be that we do not record the part of the flow which are actually of interest to the analysis. 

```r
### delete Water Years for which less than 80% of the data is present (~292 days complete)
### crude way to be sure we are not missing some important events
tapply(st72014$Flow, factor(st72014$WaterYear), length)
```

```
## 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 
##  243  361  210  344  366  365  360  365  366  365  365  358  366   97  273 
## 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 
##  365  366  365  365  365  366  365  365  365  366  329  365  365  366  365 
## 2005 2006 2007 2008 2009 2010 2011 2012 
##  365  365  365  365  365  365  366  363
```

```r
st72014 <- st72014[!st72014$WaterYear %in% c(1975,1977,1988,1989), ]
```

To extract the annual maxima, we use the `inst2daily` function, which needs as input a matrix of data, and the information of which column stores the variable of interest (`indVec`) and which column stores the (factor) variable which should be used to divide the data (`indFac`):


```r
names(st72014)
```

```
## [1] "Date"      "Flow"      "WaterYear" "Summer"
```

```r
anMax <- inst2daily(zz = st72014, indFac = 3, indVec = 2)
## annual maxima of Flow (2nd column) extracted for each WaterYear (3rd column)
summary(anMax)
```

```
##       Date                          Flow         WaterYear 
##  Min.   :1976-11-29 00:00:00   Min.   : 2.62   1976   : 1  
##  1st Qu.:1986-03-17 18:00:00   1st Qu.: 6.73   1978   : 1  
##  Median :1996-12-29 23:00:00   Median : 8.34   1979   : 1  
##  Mean   :1995-12-03 18:42:21   Mean   : 8.63   1980   : 1  
##  3rd Qu.:2004-08-17 12:15:00   3rd Qu.: 9.50   1981   : 1  
##  Max.   :2012-12-06 00:00:00   Max.   :15.92   1982   : 1  
##                                                (Other):28  
##      Summer     
##  Min.   :0.000  
##  1st Qu.:0.000  
##  Median :0.000  
##  Mean   :0.265  
##  3rd Qu.:0.750  
##  Max.   :1.000  
## 
```

```r
### incidentally, the indFac variable doesn't need to be sorted
anMax2 <- inst2daily(zz = st72014[sample(seq_along(st72014$Date)),], indFac = 3, indVec = 2)
identical(anMax,anMax2)
```

```
## [1] TRUE
```
Note that the output of `inst2daily` is a matrix with the same columns of the original matrix. The function simply selects the day in which the maximum was recorded, so if any additional information is available in the matrix (e.g. daily rainfall, temperature, etc.) the information will be kept in the matrix. 

We now have a dataset with all annual maxima recorded at the station. We can easily extract the seasonal maxima as well and display them in a plot. 



```r
#### extract seasonal maxima
winMax <- inst2daily(zz = st72014[st72014$Summer != 1,], indFac = 3, indVec = 2)
sumMax <- inst2daily(zz = st72014[st72014$Summer == 1,], indFac = 3, indVec = 2)

plot(st72014$Date, st72014$Flow, type="l")
points(anMax$Date, anMax$Flow, pch=16)
points(winMax$Date, winMax$Flow, pch=4, col=4)
points(sumMax$Date, sumMax$Flow, pch=4, col=2)
```

![plot of chunk unnamed-chunk-8](./README_files/figure-html/unnamed-chunk-8.png) 

```r
### the outputs of inst2daily are ordered, no need to merge
allMax <- cbind(anMax$Flow, winMax$Flow, sumMax$Flow)
```
Different methods are used to fit a distribution to the data. The maximum likelihood and the L-moment approaches are probably the most frequently used ones. We use the `ismev` package to estimate the distribution parameters via maximum likelihood, and the `lmom`  package for L-moments. Other useful packages are `evd` or `nsRFA`, which could be also used to obtain Maximum likelihood or L-moment estimates of the parameters of a distribution.












