#Make dyad data set. This is to include:
# related probability, genetic distance, geographic distance, hatchDiff, Year difference, chemoDist.
# The final data frame will/may also have columns common mum and common dad (T/F). This is to be stored in Maccullochella database.
#It should be about 243^2 /2 - 243/2 (29403rows,8cols)

require(dplyr)
library(tidyr)
library(reshape2)

#related probability
load("output 3 Years Larval snps")
dfRelatedness<-data.frame(output$relatedness) 

#Give row names to then enable dit matrices with names.
dm<-qslMPeeliiForRelated
row.names(dm)<-dm[,1]

#genetic distance
df<-larvalPeeliSnps; main<-"3 Years Larval snps"
df<-df[-c(1:4),-c(1:17)]
df<-data.frame(t(df))
row.names(df)<-df[,1]
df[,1]<-NULL
genDiff<-dist(df)
summary(genDiff)

genDifflong <- melt(as.matrix(genDiff), varnames = c("row", "col"))
genDifflong<-genDifflong[genDifflong$row > genDifflong$col,]

names(genDifflong)[names(genDifflong)=="row"] <- "larvaA"
names(genDifflong)[names(genDifflong)=="col"] <- "larvaB"
names(genDifflong)[names(genDifflong)=="value"] <- "genDiff"
        
        
#geo distance difference nestEst
dm1<-select(dm,Distance.to.Angle.Crossing..m.,estimatedAge,Day.of.Year)
dm1[,4]<-dm$Distance.to.Angle.Crossing..m.-(dm$Day.of.Year-dm$estimatedAge*1200/2)
dm1<-select(dm1,V4)

geoDiff<-dist(dm1)
summary(geoDiff)

geoDifflong <- melt(as.matrix(geoDiff), varnames = c("row", "col"))
geoDifflong<-geoDifflong[geoDifflong$row > geoDifflong$col,]

names(geoDifflong)[names(geoDifflong)=="row"] <- "larvaA"
names(geoDifflong)[names(geoDifflong)=="col"] <- "larvaB"
names(geoDifflong)[names(geoDifflong)=="value"] <- "geoDiffMetres"

#hatch date difference
dm2<-select(dm,hatchedDoY) #should be days between dates proper.
hatchDiff<-dist(dm2)
summary(hatchDiff)

hatchDifflong <- melt(as.matrix(hatchDiff), varnames = c("row", "col"))
hatchDifflong<-hatchDifflong[hatchDifflong$row > hatchDifflong$col,]

names(hatchDifflong)[names(hatchDifflong)=="row"] <- "larvaA"
names(hatchDifflong)[names(hatchDifflong)=="col"] <- "larvaB"
names(hatchDifflong)[names(hatchDifflong)=="value"] <- "hatchDiffDays"

#Make Year Difference (should be 0,1 or 2)
dm1<-select(dm,YearOnly)
yearDiff<-dist(dm1)
summary(yearDiff)

yearDifflong <- melt(as.matrix(yearDiff), varnames = c("row", "col"))
yearDifflong<-yearDifflong[yearDifflong$row > yearDifflong$col,]

# require(dplyr)#NOT working
# dplyr::rename(yearDifflong, row = larvaA)
# dplyr::rename(yearDifflong, col = larvaB)
#so old way below

names(yearDifflong)[names(yearDifflong)=="row"] <- "larvaA"
names(yearDifflong)[names(yearDifflong)=="col"] <- "larvaB"
names(yearDifflong)[names(yearDifflong)=="value"] <- "yearsDiff"

names(dfRelatedness)[names(dfRelatedness)=="ind1.id"] <- "larvaA"
names(dfRelatedness)[names(dfRelatedness)=="ind2.id"] <- "larvaB"

dyads<-merge(yearDifflong,geoDifflong, by.x=c("larvaA","larvaB"))
dyads<-merge(dyads,hatchDifflong, by.x=c("larvaA","larvaB"))
dyads<-merge(dyads,genDifflong, by.x=c("larvaA","larvaB"))

dyadsm1<-merge(dyads,dfRelatedness, by.x=c("larvaA","larvaB"))
dyadsm2<-merge(dyads,dfRelatedness, by.x=c("larvaB","larvaA"))
dyads<-rbind(dyadsm1,dyadsm2)


#Add in Mum and Dad by two details.

dyads<-dyads %>% left_join(dm, by =c("larvaA"="LarvaID"))
dyads<-dyads %>% left_join(dm, by =c("larvaB"="LarvaID"))

#eg: then
#select(dyadsm, larvaA,larvaB, mumPulldown.x, mumPulldown.y,yearsDiff,hatchDiffDays)
#Write the data frame out to file for import into the database.
look<-select(dyads, larvaA,larvaB,yearsDiff,mumPulldown.x,Fathers.x,mumPulldown.y,Fathers.y,hatchDiffDays,geoDiffMetres,genDiff,trioml,wang,lynchli,lynchrd,ritland,quellergt,dyadml,Delta13C.x,Delta15N.x,CNRatio.x,Delta13C.y,Delta15N.y,CNRatio.y)

write.csv(dyads,file = "outData/allTheDyads.csv")


#Also, in case it is needed
#ageDiff
dm1<-select(dm,estimatedAge)
ageDiff<-dist(dm1)


##########################trying to fix problem mentioned at line 11
#before
# setdiff(dyads$larvaA,dfRelatedness$larvaA)
# setdiff(dyads$larvaB,dfRelatedness$larvaB)
# setdiff(dfRelatedness$larvaA,dyads$larvaA)
# setdiff(dfRelatedness$larvaB,dyads$larvaB)
# 
# sort(unique(dfRelatedness$larvaA))
# sort(unique(dfRelatedness$larvaB))
# 
# hold<-filter(dyads,larvaA==351 & larvaB==185)
# holdA<-filter(dyads,larvaA==377 & larvaB==70)
# 
# holdOut<-rbind(hold,holdA)
# 
# hold2<-filter(dyads,!(larvaA==351 & larvaB==185))
# holdB<-filter(hold2,!(larvaA==377 & larvaB==70))
# 
# 
# head(holdOut)
# names(holdOut)[names(holdOut)=="larvaA"] <- "larvaB1"
# names(holdOut)[names(holdOut)=="larvaB"] <- "larvaA"
# names(holdOut)[names(holdOut)=="larvaB1"] <- "larvaB"
# holdOut<-holdOut[,c(2,1,3:ncol(hold))] # reorder correctly
# head(holdOut)
# 
# 
# #second pair
# head(holdA)
# names(holdA)[names(holdA)=="larvaA"] <- "larvaB1"
# names(holdA)[names(holdA)=="larvaB"] <- "larvaA"
# names(holdA)[names(holdA)=="larvaB1"] <- "larvaB"
# holdA<-holdA[,c(2,1,3:ncol(hold))] # reorder correctly
# head(holdA)
# 
# letsSee<-rbind(holdOut,holdB)
# 
# 
# #after
# # setdiff(letsSee$larvaA,dfRelatedness$larvaA)
# # setdiff(letsSee$larvaB,dfRelatedness$larvaB)
# # setdiff(letsSee$larvaA,dyads$larvaA)
# # setdiff(letsSee$larvaB,dyads$larvaB)
# # 
# # setdiff(letsSeeAB$larvaA,dfRelatedness$larvaA)
# # setdiff(letsSeeAB$larvaB,dfRelatedness$larvaB)
# # setdiff(letsSeeAB$larvaA,dyads$larvaA)
# # setdiff(letsSeeAB$larvaB,dyads$larvaB)
# 
# 
# #letsSee<-letsSee[-c(which(duplicated(letsSee)==TRUE)), ]
# dyadsm<-merge(letsSee,dfRelatedness, by.x=c("larvaA","larvaB"))
# dyadsm2<-merge(dfRelatedness,letsSee, by.x=c("larvaB","larvaA"))
# dyadsm3<-rbind(dyadsm,dyadsm2)
# 
# dyadsm4<-merge(dyadsm3,dfRelatedness, by=c("larvaA","larvaB"))
# 
# dfRelatedness$larvaA<-as.integer(dfRelatedness$larvaA)
# dyadsm<-inner_join(letsSee, dm, by =c("larvaA"="larvaID"))
# dyads<-merge(yearDifflong,geoDifflong, by.x=c("larvaA","larvaB"))
# dyads<-merge(dyads,hatchDifflong, by.x=c("larvaA","larvaB"))
