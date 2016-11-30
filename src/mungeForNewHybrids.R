#Basically to prepare file types for newHybrid input

source("src/dart.r")
source("src/read.dart.r")

#library(devtools)
#dev_mode(on=T)
#dev_mode(on=F)

#Install non- cran versions (the older version of adgenet has a serious flaw and is not updated yet.)
#install_github("thibautjombart/adegenet")
#install_github("green-striped-gecko/PopGenReport")


all.dart <- read.dart("OtherData/allDArTsnps.csv", topskip = 5)
gl.dart <- dart2genlight(all.dart, covfilename = "OtherData/qslDartCovariatesAll.csv")


##Below is not required but does Fst and PCOA
library(StAMPP)
system.time(snpfst <-stamppFst(gl.dart,nboots=1, percent=95, nclusters=8 ))


#pcoa
system.time(pca1 <- glPca(gl.dart  , parallel=F, nf=3))
s.class(pca1$scores, pop(gl.dart), col=rainbow(nlevels(pop(gl.dart))))