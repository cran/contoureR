% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/RcppExports-Doc.R
\name{convexHullAM}
\alias{convexHullAM}
\alias{convexHullAM_Indexes}
\alias{convexHullAM_Points}
\title{Convex Hull via Andrews Monotone, Rcpp Interface to C++ Routine}
\usage{
convexHullAM_Indexes(x, y, includeColinear=FALSE,zeroBased = TRUE)
convexHullAM_Points(x, y,includeColinear=FALSE)
}
\arguments{
\item{x}{NumericVector of x values}

\item{y}{NumericVector of y values}

\item{includeColinear}{whether to inlude points that line \strong{ON} the hull, by default this is set to FALSE, as this is
the true definition of the convex hull.}

\item{zeroBased}{Whether the return indexes should be zero based (true, for use in C++),
or One-Based (false, for use in R).}
}
\value{
\code{convexHullAM_Indexes} returns an integer vector of the indexes of the points,
whilst \code{convexHullAM_Points} returns an \code{n x 2} matrix of the points themselves.
}
\description{
This function is the R interface to the C++ implementation of Andrews Monotone, a well known algorithm for solving
the convex hull in \code{O(nlogn)} time complexity.
}
\examples{
library(contoureR)
library(ggplot2)
set.seed(1)
x  = runif(100)
y  = runif(100)
ch = convexHullAM_Indexes(x,y,includeColinear=FALSE,zeroBased = FALSE)
ggplot(data.frame(x,y),aes(x,y)) +
 geom_point() +
 geom_path(data=data.frame(x,y)[ch,],colour="red")
}

