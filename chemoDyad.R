require(dplyr)
library(tidyr)
library(reshape2)

dfChemoType<-read.csv(file="outData/coreChem.csv")

rownames(dfChemoType)<-dfChemoType$LarvalID
dfChemoType<-select(dfChemoType,Li:Delta15N)

toDrop <- c("221","248") #not Murray cod but got chem sampled.
dfChemoType<-dfChemoType[ !(rownames(dfChemoType) %in% toDrop), ] 

chemoType<-dist(dfChemoType)
summary(chemoType)

chemoTypelong <- melt(as.matrix(chemoType), varnames = c("row", "col"))
chemoTypelong<-chemoTypelong[chemoTypelong$row > chemoTypelong$col,]

names(chemoTypelong)[names(chemoTypelong)=="row"] <- "larvaA"
names(chemoTypelong)[names(chemoTypelong)=="col"] <- "larvaB"
names(chemoTypelong)[names(chemoTypelong)=="value"] <- "chemoType"

test<-merge(chemoTypelong,dyads, by.x=c("larvaA","larvaB"))





# test2<-merge(chemoTypelong,dyads, by.x=c("larvaA","larvaB"))
# test3<-rbind(test1,test2)

# library(clusterSim)
# 
# forMant<-select(test3,genDiff,chemoType)
# 
# 
# test4<-mantel.rtest(as.dist(data.Normalization(test3$genDiff,type="n4",normalization="column")),as.dist(data.Normalization(test3$chemoType,type="n4",normalization="column")), nrepet = 9999)
# # test4



# genDist<-dist(test[,c(1:21076)])
# chemDist<-dist(test[,c(21077:21086)])
