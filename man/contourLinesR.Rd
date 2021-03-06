% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/contourLines.R
\name{contourLinesR}
\alias{contourLinesR}
\title{Get Contour Lines (list)}
\usage{
contourLinesR(x, y, z, nlevels = 10, levels = pretty(range(z, na.rm =
  TRUE)), ...)
}
\arguments{
\item{x}{Numeric data for x and y coordinate, a single matrix or data-frame object can be provided for
\code{x}, which will be used in preference to the y and z arguments. These do \strong{NOT} need to be in any particular order
nor do they need to be regular.}

\item{y}{Numeric data for x and y coordinate, a single matrix or data-frame object can be provided for
\code{x}, which will be used in preference to the y and z arguments. These do \strong{NOT} need to be in any particular order
nor do they need to be regular.}

\item{z}{numeric Data for z coordinate (the coordinate to model)}

\item{nlevels}{An integer number of bins to split the data into \strong{iff} \code{levels} or \code{binwidth}
have not been specified.}

\item{levels}{A numeric vector of the explicitly specified levels (z values) to contour, by specifying this argument,
it will override \code{nlevels} and/or \code{binwidth}. If this argument is provided, the stacking order of the
contours will be preserved in the order of first occurence within the supplied vector.}

\item{...}{any other parameters passed through to \code{\link{getContourLines}}}
}
\value{
A list of contours is returned, Each contour is a list with the elements:
\item{\code{level}}{The contour of the level}
\item{\code{x}}{The x-coordinates of the contour}
\item{\code{y}}{The y-coordinates of the contour}
}
\description{
A wrapper to the \code{\link{getContourLines}} function, provided to ease the transition from
the \code{\link[grDevices]{contourLines}} function as part of \code{\link{grDevices}}
}
\details{
This function returns data in the same format/structure as \code{\link[grDevices]{contourLines}},
ie, list of lists, which is different from the (preferred) \code{dataframe} object returned
by the \code{\link{getContourLines}} function as part of the present work.
}
\examples{
library(contoureR)
library(ggplot2)
x = runif(100)
y = runif(100)
df = expand.grid(x=x,y=y)
z  = with(df,x+y)
result = contourLinesR(df$x,df$y,z)
}
\seealso{
\code{\link{getContourLines}}
}

