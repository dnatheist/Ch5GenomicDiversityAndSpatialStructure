require(dplyr)
larvalPeeliSnps<-data.frame(t(larvalPeeliSnps))

eltwLarvList<-filter(larv,YearOnly==2011|YearOnly==2012)
eltwLarvList<-select(eltwLarvList, LarvalRecords_LarvaID)
eltwLarvList<-eltwLarvList$LarvalRecords_LarvaID
eltwlarvalPeeliSnps1<-larvalPeeliSnps[1:17,]
eltwlarvalPeeliSnps2<-larvalPeeliSnps[which(larvalPeeliSnps$X1%in% eltwLarvList),]
eltwlarvalPeeliSnps<-rbind(eltwlarvalPeeliSnps1,eltwlarvalPeeliSnps2)
eltwlarvalPeeliSnps<-data.frame(t(eltwlarvalPeeliSnps))

twthLarvList<-filter(larv,YearOnly==2012|YearOnly==2013)
twthLarvList<-select(twthLarvList, LarvalRecords_LarvaID)
twthLarvList<-twthLarvList$LarvalRecords_LarvaID
twthlarvalPeeliSnps1<-larvalPeeliSnps[1:17,]
twthlarvalPeeliSnps2<-larvalPeeliSnps[which(larvalPeeliSnps$X1%in% twthLarvList),]
twthlarvalPeeliSnps<-rbind(twthlarvalPeeliSnps1,twthlarvalPeeliSnps2)
twthlarvalPeeliSnps<-data.frame(t(twthlarvalPeeliSnps))

elthLarvList<-filter(larv,YearOnly==2011|YearOnly==2013)
elthLarvList<-select(elthLarvList, LarvalRecords_LarvaID)
elthLarvList<-elthLarvList$LarvalRecords_LarvaID
elthlarvalPeeliSnps1<-larvalPeeliSnps[1:17,]
elthlarvalPeeliSnps2<-larvalPeeliSnps[which(larvalPeeliSnps$X1%in% elthLarvList),]
elthlarvalPeeliSnps<-rbind(elthlarvalPeeliSnps1,elthlarvalPeeliSnps2)
elthlarvalPeeliSnps<-data.frame(t(elthlarvalPeeliSnps))

larvalPeeliSnps<-data.frame(t(larvalPeeliSnps))

#tidy up temporary files
rm("eltwlarvalPeeliSnps1","eltwlarvalPeeliSnps2","twthlarvalPeeliSnps1","twthlarvalPeeliSnps2","elthlarvalPeeliSnps1","elthlarvalPeeliSnps2")
