#For creating Gephi files. Need to have fields in right order too.

require(dplyr)
test<-cbind(output$relatedness,output$delta7)
View(test)
names(test)
test2<-test[c(-12,-13,-14,-15)]
test2<-merge(test2,output$inbreeding,by.x = "ind1.id", by.y = "ind.id")
names(test2)
head(test2)
test2<-rename(test2, delta7 = trioml.1)
test2<-rename(test2, source = ind1.id)
test2<-rename(test2, target = ind2.id)
test2<-rename(test2, inbredCoef = L3)
test2$LH <- NULL
test2$LR <- NULL
test2$L2 <- NULL
test2$from<-as.numeric(test2$from)
test2<-test2[order(test2$from),]
test2$type<-"undirected"

names(test2)
head(test2)

plot(test2$trioml,test2$delta7)
plot(test2$trioml,test2$inbredCoef)

plot(test2$inbredCoef,test2$delta7)

edges<-test2
write.csv(edges,file = "relatedEdgesForGephi.csv")

nodes<-test2
nodes$Label<-nodes$source
write.csv(nodes,file = "relatedNodesForGephi.csv")
