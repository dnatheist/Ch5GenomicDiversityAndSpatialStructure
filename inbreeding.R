#Inbreeding 
require(ggplot2)
load("outData/2011-2013 Larval snps coancestoryOutput")
inbreds<-data.frame(output$inbreeding)
inb<-merge(inbreds,qslMPeeliiForRelated, by.x="ind.id", by.y = "LarvaID")

inbredPlot<-qplot(inb$mumPulldown,inb$LH, xlab="",ylab="Inbreeding Coefficient")
inbredPlot +
        theme(  axis.text.x = element_text(angle = 90, hjust = 1),
                panel.background = element_rect(fill = NA, colour="grey")
                )
                
        

        