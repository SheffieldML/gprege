\name{compareROC}
\Rdversion{1.0}
\alias{compareROC}
\title{Make ROC plots.}
\description{
  This rocStats wrapper superimposes ROC curves on a plot to analyse the output performance of a method-A, and optionally compare it with that
  of a method-B, based on some ground thruth labels.
}
\usage{
  compareROC(output, groundTruthLabels, compareToRanking)
}
\arguments{
  \item{output}{(vector) The output of ranking scores returned by method-A for each data-point.} 
  \item{groundTruthLabels}{(vector) A binary vector that contains the ground truth (e.g. which genes are members of the top-100 ground truth list).}
  \item{compareToRanking}{A matrix where each column is the output vector of ranking scores returned by another competing method.}
}
\value{
  \item{area}{A scalar. The area under the ROC curve of method-A.}
}
\seealso{
  \code{ \link{rocStats} }
}
\examples{
  data(FragmentDellaGattaData) ## Load demo data.
  compareROC(output= rnorm(length(DGatta_labels_byTSNI))>0, groundTruthLabels=DGatta_labels_byTSNI)
}
