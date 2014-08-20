
#' @title Length without NAs
#' @name lengthNOTna
#' @param x a vector
#' @return The length of the vector when NAs are excluded 
#' @examples
#' x <- c(5,12,NA,4)
#' lengthNOTna(x)
#' @export
lengthNOTna <- function(x) length(x[!is.na(x)])


#' @title Mean without NAs
#' @name meanNOTna
#' @param x a vector
#' @return The mean of the vector when NAs are excluded, just \code{mean(x, na.rm=TRUE)}
#' @examples
#' x <- c(5,12,NA,4)
#' meanNOTna(x)
#' @export
meanNOTna <- function(x) mean(x, na.rm = TRUE)

