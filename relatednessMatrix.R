#This converts the fish pairs (dyads) into a matrix making it easier to identify family groups. It exports a csv file allowing colouring in excel

mdata<-output$relatedness[,c(2,3,6)]

require(tidyr)
data_wide <- spread(mdata, ind2.id, wang)
dwm<-as.matrix(data_wide)
rm(data_wide)
write.table(dwm,"relProbs.csv", sep=",", col.names = TRUE, row.names = FALSE ,quote = FALSE)

#groupMatrix
groups<-output$relatedness[,c(2,3,4)]
data_wide <- spread(groups, ind2.id, group)
groupsMat<-as.matrix(data_wide)
rm(data_wide)
write.table(groupsMat,"relGroups.csv", sep=",", col.names = TRUE, row.names = FALSE ,quote = FALSE)

#Perhaps this is all done (groups - all individuals in relation to fish number 185)
#groups<-output$relatedness[,c(3,4)]
#groups<-groups[1:285,]
