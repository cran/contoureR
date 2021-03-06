% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/orderPoints.R
\name{orderPoints}
\alias{orderPoints}
\title{Order Points Clockwise or Counter-Clockwise}
\usage{
orderPoints(x, y, ..., xm = mean(range(x)), ym = mean(range(y)),
  clockwise = TRUE)
}
\arguments{
\item{x}{numeric vector of x values}

\item{y}{numeric vector of y values of same length as x}

\item{...}{not used}

\item{xm}{the x value of the reference point}

\item{ym}{the y value of the reference point}

\item{clockwise}{order in clockwise or anticlockwise manner}
}
\description{
Returns the indexes of supplied points, \code{x} and \code{y}, ordered either clockwise or
anticlockwise about another point, which by default is taken to be the non-weighted midpoint
of the supplied data
}
\examples{
#Generate a random set of points and put them clockwise order
set.seed(1)
x  = runif(100)
y  = runif(100)
op    = orderPoints(x,y)

#To demonstrate, Lets view the points in order
library(ggplot2)
df    = data.frame(x,y)
df    = df[op,];
df$id = 1:nrow(df)
ggplot(data=df,aes(x,y,colour=id)) +
    geom_path() + geom_point() +
    scale_colour_gradient(low="green",high="red")
}

