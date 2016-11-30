# This is based on Arthurs 'create a genlight object' script.

#Not yet working. Need to sort our correct file inputs and larv variables.

source("src/generic/dart.r")
source("src/generic/gl.read.dart.r")

trial<-gl.read.dart(datafile="data/Report-DMac15-1861.csv", topskip=5, nmetavar=17, nas="-", ind.metafile="data/qslAllLarvaInfo.csv")

gl1 <- gl.recode.pops(gl, pop.recode="Mygeo_pop_recode.csv")

gl1 <- gl.filter.callrate(gl1,by="loc", 0.975)
gl1 <- gl.filter.repavg(gl1,0.99)

df <- gl2related(gl1)
rtd <- coancestry(df, wang=1)

relations <- rtd$relatedness[,c(2,3,6)]
colnames(relations) <- c("from","to","weight")

entities <- cbind(gl1$ind.names, as.character(pop(gl1)))
colnames(entities) <- c("name","pop")
