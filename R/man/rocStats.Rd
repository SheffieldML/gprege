\name{rocStats}
\Rdversion{1.0}
\alias{rocStats}
\title{Make ROC curve data.}
\description{
  Computes the points on an ROC curve by varying a threshold on the sorted outputs of the method in question.
}
\usage{
  rocStats(outputs, groundTruthLabels, decreasing = TRUE)
}
\arguments{
  \item{outputs}{A numeric vector with the outputs of the evaluated method (e.g. likelihoods from gprege).} 
  \item{groundTruthLabels}{A binary vector than contains the ground truth (e.g. which genes belong in the top-100 ground truth list).}
  \item{decreasing}{Logical. TRUE sorts outputs by decreazing order.}
}
\value{
  \item{stats}{A list of numeric variables with the necessary statistics to compute an ROC curve, a precision-recall curve, etc.}
}
\examples{
  ## see compareROC.R
  data(FragmentDellaGattaData) ## Load demo data.
  rocStats(gpregeOutput$rankingScores, DGatta_labels_byTSNItop100, decreasing=TRUE)
}
