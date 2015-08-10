#' Get Convex Hull of Points
#' 
#' Returns the sequence of indexes from supplied points x and y, that describe the convex hull of those points. 
#' With the exception of a few brief checks, is almost a direct wrapper to the \code{\link[geometry]{convhulln}} 
#' function as part of the geometry package.
#' @inheritParams getDelaunayMesh
#' @param closeHull whether to close the hull or not, by closed, meaning that the last index in the hull ix exactly the
#' same as the first index
#' @rdname getConvexHull
#' @examples
#' #Generate the Convex Hull of a Series of Points
#' set.seed(1)
#' x  = runif(100)
#' y  = runif(100)
#' ch = getConvexHull(x,y,closeHull=TRUE)
#' 
#' #To demonstrate, Lets view the hull
#' library(ggplot2)
#' df = data.frame(x,y)
#' ggplot(data=df,aes(x,y)) + 
#'    geom_path(data=df[ch,]) + 
#'    geom_point()
getConvexHull <- function(x,y,closeHull=TRUE){
  if(!all(is.numeric(x),is.numeric(y))) stop('x and y must be numeric')
  if(length(x) != length(y)) stop('x and y must be the same length')
  ret = unique(as.vector(convhulln(data.frame(x,y))))
  ret = ret[orderPoints(x=x[ret],y=y[ret])]
  if(length(ret) > 1 & closeHull)
    ret = c(ret,ret[1])
  ret
}

