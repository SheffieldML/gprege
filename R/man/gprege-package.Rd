\name{gprege-package}
\alias{gprege-package}
\docType{package}
\title{
  gprege - Gaussian Process Ranking and Estimation of Gene Expression.
}
\description{
  This package implements the method of Kalaitzis and Lawrence (2011)
  for Gaussian process modelling gene expression time-series data. The method can be
  used to filter quiet genes and quantify differential expression in time-series expression ratios.
}
\details{
  \tabular{ll}{
    Package: \tab gprege\cr
    Type: \tab Package\cr
    Version: \tab 0.99.0\cr
    Date: \tab 2011-07-08\cr
    License: \tab A-GPL Version 3\cr
  }
  For details of using the package please refer to the Vignette.
}
\author{
  Alfredo A. Kalaitzis

  Maintainer: Alfredo A. Kalaitzis <alkalait@gmail.com>
}
\references{
  A.~A.~Kalaitzis and N.~D.~Lawrence.
  A Simple Approach to Ranking Differentially Expressed Gene Expression Time Courses through Gaussian Process Regression
  \emph{BMC Bioinformatics} 2011, 12:180.
  DOI:10.1186/1471-2105-12-180.
}
\keyword{ package }
\seealso{
  \code{
    \link{demGpCov2D}, \link{demGpSample}, \link{demInterpolation},
    \link{demOptimiseGp}, \link{demRegression}
  }
}
\examples{
## see demTp63Gp1.R
}
