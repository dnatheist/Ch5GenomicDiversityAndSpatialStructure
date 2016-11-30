


mdf <- melt(newdata)
ggplot(mdf) +
        geom_density(aes(x = relvalues, color = Relationship))

#To extract intersect points off the density plot
intersect(newdata$Full, newdata$Half) #http://stackoverflow.com/questions/21212352/find-two-densities-point-of-intersection-in-r/21213177#21213177

#Compare estimators with each other 

plot(output$relatedness[,5:11]) # the wang estimator seems to give best correlation with the ML estimators.

plot(output$inbreeding[,2:3]) # corelation between inbreeding estimators. Should include the ML methods but does not output them for some reason despite the doco saying it should.

hist(output$inbreeding$LH)
hist(output$inbreeding$LR)

#remove to free memory for simulation
rm(eltwlarvalPeeliSnps)
rm(goodDArTsnps)
rm(larv)
rm(larvalPeeliSnps)
rm(qslAllLarvaInfo)
rm(qslAllLarvaInfoApr2016)
rm(qslMPeeliiForRelated)
rm(Report.DMac15.1861)
rm(thlarvalPeeliSnps)
rm(twthlarvalPeeliSnps)
