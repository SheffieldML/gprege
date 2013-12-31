compareROC <- function(output, groundTruthLabels, compareToRanking) {

## COMPAREROC Wrapper for rocStats that prints ROC plots.
## FORMAT
## DESC This rocStats wrapper superimposes ROC curves on a plot to analyse
## the output performance of a method-A, and optionally compare it with that
## of a method-B, based on some ground thruth labels.
## ARG output : The output (or ranking score) returned by method-A for each
## data-point.
## ARG groundTruthLabels : Binary vector that contains the ground truth.
## ARG compareToRanking : The output (or ranking score) returned by method-B
## for each data-point.
## RETURN area : Area under the ROC curve of method-A.
##
## USAGE : compareROC(diffLMLs, DGatta_labels_byTSNI, compareToRanking)
##
## SEEALSO : rocStats
##
## COPYRIGHT : Alfredo A. Kalaitzis, 2010, 2011
##
## GPREGE

ArocStats <- rocStats(output, groundTruthLabels, decreasing=TRUE) ## Compute points of the ROC curve method-A performance.

#graphics.off(); dev.new();
par(ps=15)
plot(0:1, 0:1, type = "n", xlab='FPR', ylab='TPR', main='ROC curve') ## Close all devices, open new device, new plot.

ArocStats$FPR[is.nan(ArocStats$FPR)] <- 0
ArocStats$TPR[is.nan(ArocStats$TPR)] <- 0
x <- ArocStats$FPR
y <- ArocStats$TPR
lines(x, y, col='darkred', lwd=6) ## Plot ROC curve of method A.

rocArea = c()
rocArea[1] <- sum(diff(x)*(y[-length(y)]+y[-1])/2) ## Trapezoidal numerical integration.
area <- rocArea[1]

## Plot ROC curves for method-B 
if (!missing(compareToRanking)) {
  lstyle <- c(3, 2, 4) ## Dotted, dashed, dotdash.
  lcolors <- c('darkblue', 'darkgreen', 'darkmagenta')

  for (f in 1:dim(compareToRanking)[2]) {
    BrocStats <- rocStats(compareToRanking[,f], groundTruthLabels, decreasing=TRUE)
    BrocStats$FPR[is.nan(BrocStats$FPR)] <- 0; BrocStats$TPR[is.nan(BrocStats$TPR)] <- 0
    x <- BrocStats$FPR; y <- BrocStats$TPR
    rocArea[1+f] <- sum(diff(x)*(y[-length(y)]+y[-1])/2) ## Trapezoidal numerical integration.
    lines(x, y, lty=lstyle[f], col=lcolors[f], lwd='3')
  }

  legend('bottomright', legend = as.expression(c(
    paste('GP (auc=', as.character(signif(rocArea[1],3)), ')', sep=''),
    bquote(paste(BATS[G], ' (auc=', .(as.character(signif(rocArea[2],3))), ')', sep='')),
    bquote(paste(BATS[T], ' (auc=', .(as.character(signif(rocArea[3],3))), ')', sep='')),
    bquote(paste(BATS[DE], ' (auc=', .(as.character(signif(rocArea[4],3))), ')', sep=''))
    )), col=c('darkred', lcolors), pch=rep(-1,4), lty=c(1, lstyle), lwd=c(6,3,3,3))
} else {
  legend('bottomright', legend=c(paste('GP (auc=', as.character(signif(rocArea[1],3)), ')',sep='')),
    col='darkred', pch=-1, lty=1, lwd=6)
}
return(area)
}
