## @knitr createNet

require(igraph)

relData<-output$relatedness[,c(2,3,estimator)]

#to see how the data is distributed
#hist(relData[,3], main="Histogram of estimator") 

#Makes links DF
lrelData<-relData
colnames(lrelData)[1] <- "from"
colnames(lrelData)[2] <- "to"
colnames(lrelData)[3] <- "weight"
lrelData$type<-"probRel"

relDataNoRows<-nrow(relData) #get number of rows and store for later use.
#create new DF for nodes
nrelData <- data.frame(matrix(ncol = 4, nrow = relDataNoRows))
colnames(nrelData) <- c("id", "name", "type", "label")

#Populate the DF with the 'to' nodes from original list
nrelData[,1]<-relData[,2]
nrelData[,2]<-relData[,2]

#add the first entry from the 'from' list (otherwise it wont get in as it is not
# in the 'to list')
newrow = c(relData[1,1],relData[1,1],NA,NA)
nrelData = rbind(nrelData,newrow)

# nrow(nrelData); length(unique(nrelData$id))
# nrow(lrelData); nrow(unique(lrelData[,c("from", "to")]))

nrelData<-data.frame(unique(nrelData[ , 1:4]))

#Apply types "adult" or "larvae" to fish based on id number.
tst<-ifelse(grepl("^A", nrelData$id),nrelData$type<-"adult",nrelData$type<-"larvae")
nrelData$type<-tst
rm(tst)

# Add data from larv to relatedNodes so that they are available for filtering on year etc.
larvNrel<-merge(nrelData,larv, by.x="id", by.y="LarvalRecords_LarvaID")
larvNrel<-larvNrel[,c(1:4,66)]

#Have a look at two data frames used to make net
# head(nrelData)
# head(larvNrel)

#Create the network data
net<-graph_from_data_frame(d=lrelData, vertices=larvNrel, directed=FALSE)

#fileName=paste(c("./outData/",main,"net"),sep="", collapse=" ")
fileName=paste("./outData/",main," net",sep="", collapse=" ")
save(net,file = fileName)# save net files for next chunk on.

