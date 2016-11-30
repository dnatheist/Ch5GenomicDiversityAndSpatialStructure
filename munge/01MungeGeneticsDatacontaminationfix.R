## Load larv data frame
require(dplyr)
allIdentifiers<-dplyr::mutate(DARTallIdentifiers, identifier=paste(DARTallIdentifiers$PlateID,DARTallIdentifiers$Row,DARTallIdentifiers$Column,sep=""))

larv<-qslAllLarvaInfo #Make dataframe from larval data

## Extract more meaningful names (Labels) in place of Larval ID numbers
row.names(larv)<-larv$Label # #Make Labels row names too old line :rownames(larv) <- larv[,1]
#rm(qslAllLarvaInfo2015) #Tidy up

##Load DArT 2 row snps data frame
DArTsnps<-Report.DMac15.1861
#rm(Report.DMac15.1861)# Tidy up
##291 fish from DArT. 2x94, 1x93(plate 3) plus 8 plus 2 Trout cod controls from trial 1)

DArTsnps<-t(DArTsnps) # Just transpose this way larvae become observations and snps become the variables.
DArTsnps<-data.frame(DArTsnps)

## Create a unique identifier as DArT did not import comments with larval fish ID
DArTsnps$X1<-paste(DArTsnps$X2,DArTsnps$X3,DArTsnps$X4,sep="")

vec<-DArTsnps$X1# List of DArT 'larval identifiers' for use recoding names with recoderFunc which translates the DaRT identifiers into something more meaningful.

recoderFunc <- function(data, oldvalue, newvalue) {
  
  # convert any factors to characters
  
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  
  # create the return vector
  
  newvec <- data
  
  # put recoded values into the correct position in the return vector
  
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  
  newvec
  
}

vec<-recoderFunc(vec, allIdentifiers$identifier, allIdentifiers$Comments)

DArTsnps[,1]<-vec # Load more meaningful names (Labels) in place of Larval ID numbers
DArTsnps$X5<-as.character(DArTsnps$X5)
DArTsnps[18:308,5]<-DArTsnps[18:308,1]
rm(vec)# Tidy up and now new DArTsnps data frame now has correct fish labels 

##Create List with Trial Fish Samples
trialFish<-c(subset(DArTsnps, grepl("T", X1),X1))
trialFish<-trialFish$X1

DArTsnps<-DArTsnps[-which(DArTsnps$X1%in% trialFish),]

degradedFish<-c("221","243","248","249","251","252","253")
goodDArTsnps<-DArTsnps[-which(DArTsnps$X1%in% degradedFish),] #larvae were rotting as a result of a fresh in the river and inability to retrieve nets for a few days. Turns out 251, 252,253 were carp larvae (mitochondrial sequencing). 243, 249 have Murray cod mitos but sequenced very poorly so are also excluded. 221 and 248 also poor but unknown mitos.

#Need to remove hybrids and trout cod for 'related' BUT NOT for dart.analysis.
# Will also need to heit covariates - another reason why I need to make related a seperate project.

#contaminatedFish<- c("177","283","314","332","366")#see lab book p36 potentially but probably not contaminated.
#goodDArTsnps<-goodDArTsnps[-which(goodDArTsnps$X1%in% contaminatedFish),]

hybridFish<-c("102","106","141","145","178","262","269", "302") #from hybrids paper
allPeeliSnps<-goodDArTsnps[-which(goodDArTsnps$X1%in% hybridFish),]

TCcontrols<-c("Bendora20130408","Narranderah")
allPeeliSnps<-allPeeliSnps[-which(allPeeliSnps$X1%in% TCcontrols),]
ncol(allPeeliSnps)-17 #17 is columns up front of DArT file.

#Identify and remove adult fish
adultFish<-c(subset(DArTsnps, grepl("A", X1),X1))
adultFish<-adultFish$X1
larvalPeeliSnps<-allPeeliSnps[-which(allPeeliSnps$X1%in% adultFish),]
ncol(larvalPeeliSnps)-17 #17 is columns up front of DArT file.
adultPeeliSnps<-allPeeliSnps[which(allPeeliSnps$X1%in% adultFish),]
ncol(adultPeeliSnps)-17 #17 is columns up front of DArT file.

#Keep some files
DArTsnps<-as.data.frame(t(DArTsnps))
write.csv(DArTsnps,".\\otherData\\allDArTsnps.csv", row.names=FALSE)

goodDArTsnps<-as.data.frame(t(goodDArTsnps))
write.csv(goodDArTsnps,".\\otherData\\goodDArTsnps.csv", row.names=FALSE)

#write.csv(larvalPeeliiSnps,".\\otherData\\larvalPeeliiSnps.csv", row.names=FALSE)
#write.csv(adultPeeliiSnps,".\\otherData\\adultPeeliiSnps.csv", row.names=FALSE)
rm("Report.DMac15.1861") #tidy up
rm("DARTallIdentifiers")
rm("allIdentifiers")
rm("DArTsnps")
rm("goodDArTsnps")
rm("qslAllLarvaInfo")
rm("hybridFish")
rm("degradedFish")
rm("TCcontrols")
rm("trialFish")
