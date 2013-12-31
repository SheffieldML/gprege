demTp63Gp1 <- function(fulldataset=FALSE) {

## DEMTP63GP1 gprege on TP63 expression time-series.
## FORMAT 
## DESC Demo of Gaussian Process Regression and Estimation of Gene
## Expression on TP63 time-series data (see gprege.m).
## See Kalaitzis & Lawrence (2011) for a detailed discussion of the
## ranking algorithm and dataset used.
## ARG fulldataset : (LOGICAL) Download and use the full dataset.
##
## USAGE: demTp63Gp1(fulldataset=FALSE)
##
## COPYRIGHT: Alfredo A. Kalaitzis, 2011
##
## SEEALSO : gprege
##
## GPREGE

DGatta_labels_byTSNI <- NULL; DGatta_labels_byTSNItop100 <- NULL; exprs_tp63_RMA <- NULL
rm(DGatta_labels_byTSNI, DGatta_labels_byTSNItop100, exprs_tp63_RMA)

gpregeOptions <- list()

if (fulldataset) {
  con <- url('http://staffwww.dcs.shef.ac.uk/people/A.Kalaitzis/DellaGattaData.RData')
  load(con, envir=environment())   ## Download full Della Gatta dataset (in this function's invironment).
  close.connection(con)
  gpregeOptions$indexRange <- which(DGatta_labels_byTSNItop100)[1:2]
} else {
  data(FragmentDellaGattaData, envir=environment()) ## Load demo data.
  gpregeOptions$indexRange <- c(1:2)
}

## Download BATS rankings (Angelini, 2007)
## Case 1: Delta error prior, case 2: Inverse Gamma error prior, case 3: Double Exponential error prior
BATSranking = matrix(0, length(DGatta_labels_byTSNItop100), 3)
for (i in 1:3) { 
  tmp <- read.table(url(paste('http://arxiv.org/src/1106.4333v1/anc/DGdat_p63_case',i,'_GL.txt',sep='')), skip=1) ## Read the gene numbers
  genenumbers <- as.numeric(lapply( as.character(tmp[,2]), function(x) x=substr(x,2,nchar(x))))
  BATSranking[,i] <- tmp[sort(genenumbers, index.return=TRUE)$ix, 4] ## Sort rankings by gene numbers.
}
## The smaller the BATS ranking metric is, the better the rank that the gene reporter gets.
BATSranking = 1/BATSranking ## Invert those ranking metrics to compare on a common ground.

tTrue = matrix(seq(0,240,by=20), ncol=1)
## Setup other gprege options.
gpregeOptions$explore <- TRUE ## Explore individual profiles in interactive mode.
gpregeOptions$exhaustPlotRes <- 30 ## Exhaustive plot resolution of the LML function.
gpregeOptions$exhaustPlotLevels <- 10 ## Exhaustive plot contour levels.
gpregeOptions$exhaustPlotMaxWidth <- 100 ## Exhaustive plot maximum lengthscale.
gpregeOptions$iters <- 100
gpregeOptions$labels <- DGatta_labels_byTSNI ## Noisy ground truth labels (which genes are in the top 786 ranks of the TSNI ranking?).
gpregeOptions$display <- FALSE ## SCG optimisation: display messages.

## Matrix of different hyperparameter configurations as rows:
## [inverse-lengthscale   percent-signal-variance   percent-noise-variance].
gpregeOptions$inithypers <- matrix( c(
	1/1000,	1e-3,	0.999
#	,1/8,	0.999,	1e-3
	,1/20,	0.999,	1e-3
#	,1/30,	0.999,	1e-3
  ), ncol=3, byrow=TRUE)

graphics.off()

precalculated_rankingScores = gpregeOutput$rankingScores
gpregeOutput <- gprege(data=exprs_tp63_RMA, inputs=tTrue, gpregeOptions=gpregeOptions)
compareROC(output=precalculated_rankingScores, groundTruthLabels=DGatta_labels_byTSNItop100, compareToRanking=BATSranking)

}

