#Load support functions

source("src/dart.r")
source("src/is.fixed.r")
source("src/gl.fixed.diff.r")
source("src/gl.filter.CallRate.r")
source("src/gl.filter.RepAvg.r")
source("src/gl2nhyb.r")

library(reshape2)
library(MASS)

#Note: gl is a genlight object. You will need to have created this using Bernd’s script
#Filter stringently on callrate (Optional)

gl <- gl.filter.CallRate(gl.dart,0.99)

#Filter stringently on reproducibility

gl <- gl.filter.RepAvg(gl,0.95)

#Define parental populations

p1<-c("Narranderah","Bendora")
p0<-c("Nerreman","Bullen")

#Filter out from the master file, those loci that are not fixed between parental, and output
#those that are fixed in NewHybrids format

gl.out <- gl2nhyb(gl, outfile="posteriors.txt", p0, p1, threshold=0)
