require(dplyr)

elLarvList<-filter(larv,YearOnly==2011)
elLarvList<-select(elLarvList, LarvalRecords_LarvaID)
elLarvList<-elLarvList$LarvalRecords_LarvaID
ellarvalPeeliSnps1<-larvalPeeliSnps[1:17,]
ellarvalPeeliSnps2<-larvalPeeliSnps[which(larvalPeeliSnps$X1%in% elLarvList),]
ellarvalPeeliSnps<-rbind(ellarvalPeeliSnps1,ellarvalPeeliSnps2)
ellarvalPeeliSnps<-data.frame(t(ellarvalPeeliSnps))

twLarvList<-filter(larv,YearOnly==2012)
twLarvList<-select(twLarvList, LarvalRecords_LarvaID)
twLarvList<-twLarvList$LarvalRecords_LarvaID
twlarvalPeeliSnps1<-larvalPeeliSnps[1:17,]
twlarvalPeeliSnps2<-larvalPeeliSnps[which(larvalPeeliSnps$X1%in% twLarvList),]
twlarvalPeeliSnps<-rbind(twlarvalPeeliSnps1,twlarvalPeeliSnps2)
twlarvalPeeliSnps<-data.frame(t(twlarvalPeeliSnps))

thLarvList<-filter(larv,YearOnly==2013)
thLarvList<-select(thLarvList, LarvalRecords_LarvaID)
thLarvList<-thLarvList$LarvalRecords_LarvaID
thlarvalPeeliSnps1<-larvalPeeliSnps[1:17,]
thlarvalPeeliSnps2<-larvalPeeliSnps[which(larvalPeeliSnps$X1%in% thLarvList),]
thlarvalPeeliSnps<-rbind(thlarvalPeeliSnps1,thlarvalPeeliSnps2)
thlarvalPeeliSnps<-data.frame(t(thlarvalPeeliSnps))

allPeeliSnps<-data.frame(t(allPeeliSnps))
adultPeeliSnps<-data.frame(t(adultPeeliSnps))
larvalPeeliSnps<-data.frame(t(larvalPeeliSnps))

#tidy up temporary files
rm("ellarvalPeeliSnps1","ellarvalPeeliSnps2","twlarvalPeeliSnps1","twlarvalPeeliSnps2","thlarvalPeeliSnps1","thlarvalPeeliSnps2")
