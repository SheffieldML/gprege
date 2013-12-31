\name{exhaustivePlot}
\Rdversion{1.0}
\alias{exhaustivePlot}
\title{Plot of the LML function by exhaustive search.}
\description{
  Exhaustively searches the hyperparameter space by a grid, whose resolution is passed as an argument, and plots the LML function for every point in the space.
}
\usage{
 exhaustivePlot(y, x, xstar, options, maxwidth, res, nlevels)
}
\arguments{
  \item{y}{the target (output) data.} 
  \item{x}{the input data matrix.}
  \item{xstar}{the points to predict function values.}
  \item{options}{options structure as defined by gpOptions.m.}
  \item{maxwidth}{ maximum lengthscale to search for.}
  \item{res}{The search resolution. Number of points to plot for in the search range.}
  \item{nlevels}{Number of contour levels.}
}
\value{
  \item{area}{Area under the ROC curve of method-A.}
}
\seealso{
  \code{ \link{rocStats} }
}
\examples{
  noiseLevel <- 0.2
  noiseVar <- noiseLevel^2
  options <- gpOptions()
  options$kern$comp <- list('rbf','white')
  ## Create data set
  l <- 9; x <- matrix(seq(0,240,by=20), ncol=1)
  trueKern <- kernCreate(x, 'rbf')
  trueKern$inverseWidth <- 1/(20^2) ## Characteristic inverse-width.
  K <- kernCompute(trueKern, x) + diag(dim(x)[1])*noiseVar
  ## Sample some true function values.
  y <- gaussSamp(Sigma=K, numSamps=1)
  xTest <- as.matrix(seq(0, 240, length=200))
  graphics.off(); dev.new(); plot.new(); dev.new(); plot.new()
  exhaustivePlot(y, x, xTest, options=options, maxwidth=100, res=50, nlevels=75)
}
