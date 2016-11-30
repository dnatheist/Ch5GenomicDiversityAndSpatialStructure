## @knitr halfSibsHatchDiff

# routine to plot and identify patterns in half siblings and their hatch day differences within years to see if we can infer polyandry, polygyny etc.

#Create appropriate subset.
halfSibs<-output$relatedness
halfSibs$ind1.id<-as.integer(halfSibs$ind1.id)
halfSibs$ind2.id<-as.integer(halfSibs$ind2.id)
halfSibs<-inner_join(halfSibs, larv, by =c("ind1.id"="LarvalRecords_LarvaID"))
halfSibs<-select(halfSibs,ind1.id, ind2.id, wang,hatchDoY) #is this correct age to use?
halfSibs<-inner_join(halfSibs, larv, by =c("ind2.id"="LarvalRecords_LarvaID"))
halfSibs<-select(halfSibs,ind1.id, ind2.id, wang, hatchDoY.x, hatchDoY.y,YearOnly)
halfSibs$hatchDiff<-halfSibs$hatchDoY.x-halfSibs$hatchDoY.y

#actually filter half siblings
halfSibs<-filter(halfSibs, wang<0.4, wang>0.15)

#plot 
hist(halfSibs$hatchDiff)
hist(abs(halfSibs$hatchDiff), 12)
#plot(halfSibs$wang, halfSibs$hatchDiff)


