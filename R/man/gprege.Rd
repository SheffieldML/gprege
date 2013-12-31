\name{gprege}
\Rdversion{1.0}
\alias{gprege}
\title{Gaussian process ranking and estimation of gene expression time-series}
\description{
  Fits two GPs with the an RBF (+ noise diagonal) kernel on each profile. One GP kernel is initialised wih a short lengthscale
  hyperparameter, signal variance as the observed variance and a zero noise variance. It is optimised via scaled conjugate gradients (netlab). The
  other GP has fixed hyperparameters with a zero inverse-width, zero signal variance and noise variance as the observed variance. The log-ratio of
  marginal likelihoods of the two hypotheses acts as a score of differential expression for the profile. Comparison via ROC curves is
  performed against BATS (Angelini et.al, 2007). See Kalaitzis & Lawrence (2011) for a detailed discussion of the ranking algorithm and dataset used.
}
\usage{
  gprege(data, inputs, gpregeOptions)
}
\arguments{
  \item{data}{The matrix of gene expression profiles; one profile per row.}
  \item{inputs}{Inputs (timepoints) to the GP.}
  \item{gpregeOptions}{Options list for gprege with fields
    \describe{
      \item{explore}{Logical. TRUE operates in a user interactive mode. Used for examining individual gene expression profiles.}
      \item{labels}{A binary vector. TRUE specifies whether the corresponding profile comes from a differentially expressed gene (usually from a ground truth).}
      \item{indexRange}{A numeric vector. Range of indices of profiles on which the function should operate. Useful for selective exploration of specific profiles, e.g. only genes marked as differentially expressed in a ground truth list.}
      \item{interpolatedT}{A numeric vector. New timepoints to interpolate for each profile, based on the estimated function values.}
      \item{iters}{A scalar. The number of iterations for scaled-conjugate gradients (SCG) optimisation.}
      \item{display}{Logical. Display gradient and LML information on each SCG iteration.}
      \item{inithypers}{The matrix of hyperparameter configurations as its rows. Each row has the following format:
	[inverse-lengthscale   percent-signal-variance   percent-noise-variance]
	The first row corresponds to a (practically constant) function with a very large lengthscale. Such a function will account for 0 percent of the observed variance in the expression profile (hence 0 for signal) and explain it as noise (hence 1 for noise). Subsequent rows (initialisations for SCG optimisation) correspond to functions of various lengthscales that explain all the observed variance as signal. A reasonable lengthscale would be roughly in line with the time-point sampling intervals.}
      \item{exhaustPlotRes}{A scalar. The search resolution. Used for interactive mode (explore == 1).}
      \item{exhaustPlotLevels}{A scalar. Number of contour levels in the exhaustive plot. Used for interactive mode (explore == 1).}
      \item{exhaustPlotMaxWidth}{A scalar. the maximum lengthscale to search for. Used for interactive mode (explore == 1).}
    }
  }
}
\value{
  \item{gpregeOutput}{Output list with fields:
    \describe{
      \item{signalvar}{A numeric vector of the vertical lengthscales of the optimised RBF kernel; one for each profile.}
      \item{noisevar}{A numeric vector. Similar to signalvar, but for the noise hyperparameter.}
      \item{width}{A numeric vector. Similar to signalvar and noisever, but for the horizontal lengthscales of the RBF.}
      \item{LMLs}{A numeric vector of log-marginal likelihoods of the GP; one for each profile.}
      \item{interpolatedData}{A matrix of the extended dataset with interpolated values as the augmenting columns.}
      \item{rankingScores}{A numeric vector of the ranking scores, based on the log-ratio of marginal likelihoods.}
    }
  }
}
\seealso{
  \code{
    \link{gpOptions}, \link{gpCreate}, \link{gpExpandParam}, \link{gpOptimise},
    \link{gpExtractParam}, \link{gpLogLikelihood}, \link{gpPosteriorMeanVar}.
  }
}
\examples{
  ## see demTp63Gp1.R
  data(FragmentDellaGattaData) ## Load demo data.
  ## Setup other gprege options.
  gpregeOptions = list(indexRange=(1:2), explore=TRUE, exhaustPlotRes=30, exhaustPlotLevels=10,
    exhaustPlotMaxWidth=100, iters=100, labels=DGatta_labels_byTSNI, display=FALSE)
  ## Matrix of different hyperparameter configurations as rows:
  ## [inverse-lengthscale   percent-signal-variance   percent-noise-variance].
  gpregeOptions$inithypers <- matrix( c(
	   1/1000,	1e-3,	0.999
	  ,1/20,	0.999,	1e-3
    ), ncol=3, byrow=TRUE)
  gpregeOutput <- gprege(data=exprs_tp63_RMA, inputs=matrix(seq(0,240,by=20), ncol=1), gpregeOptions=gpregeOptions)
}
