#Use dplyr to "Group by" and get means etc for the groups.

require(dplyr)
tmp <- group_by(larv, LarvalRecords_NestID)
nestMeansDoY <- summarise(tmp, 
                   count = n(), 
                   mdoy = mean(hatchedDoY, na.rm = TRUE),
                   sddoy= sd(hatchedDoY, na.rm = TRUE),
                   mDist= mean(Distance.to.Angle.Crossing..m., na.rm = TRUE),
                   mAge=mean(estimatedAge,na.rm = TRUE),
                   nestDist=mDist - ((mAge-7)*1200/2))
nestMeansDoY

#nestMeansDoY<-data.frame(nestMeansDoY$LarvalRecords_NestID,nestMeansDoY$count,nestMeansDoY$mdoy,nestMeansDoY$sddoy)
nestMeansDoY<-data.frame(nestMeansDoY)
nestMeansDoY<-filter(nestMeansDoY,!LarvalRecords_NestID =="")

#rename(nestMeansDoY, nest=nestMeansDoY.LarvalRecords_NestID, count=nestMeansDoY.count,meanDOY=nestMeansDoY.mdoy,sdDOY=nestMeansDoY.sddoy, mDist=nestMeansDoY.mDist, sdDist=nestMeansDoY.sdDist)

write.csv(nestMeansDoY,"nests.csv")

