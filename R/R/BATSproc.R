library(affy)
library(R.matlab)
library(spam)
source("~/Documents/R/fitexprs_gpr_NDLgptoolkitver.R")

flags = read.csv('~/bats/UserData/SyntheticData1_T_1_1.txt',row.names=1,sep='')
BATSdat = read.csv('~/bats/UserData/SyntheticData1_D_1_1.txt',row.names=1,sep='')

## column: 1 2  3 4 5  6 7  8 9  10 11 12  13 14  15 16 17  18 19  20 21  22 23  24 25
## time/h: 1 1  2 2 2  4 4  6 6   8  8  8  12 12  16 16 16  20 20  24 24  28 28  32 32
## average out the replicated timepoints
# BATSdat = cbind(rowMeans(BATSdat[,1:2]), rowMeans(BATSdat[,3:5]), rowMeans(BATSdat[,6:7]),
# 		rowMeans(BATSdat[,8:9]), rowMeans(BATSdat[,10:12]),rowMeans(BATSdat[,13:14]),
# 		rowMeans(BATSdat[,15:17]),rowMeans(BATSdat[,18:19]),rowMeans(BATSdat[,20:21]),
# 		rowMeans(BATSdat[,22:23]),rowMeans(BATSdat[,24:25]))

gpdata = fitBATSExprsGPR(BATSdat, explore.mode=FALSE)

