require(igraph)

relData<-output$relatedness[,c(2,3,11)]
hist(relData$dyadml) #to see how the data is distributed

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

#Populate the DF with the nodes from original list
nrelData[,1]<-relData[,2]
nrelData[,2]<-relData[,2]

#add the first entry (otherwise it wont get in)
#newrow = c("185","185",NA,NA)
newrow = c(relData[1,1],relData[1,1],NA,NA)
nrelData = rbind(nrelData,newrow)

#have a look
head(nrelData)
head(lrelData)
nrow(nrelData); length(unique(nrelData$id))
nrow(lrelData); nrow(unique(lrelData[,c("from", "to")]))

nrelData<-data.frame(unique(nrelData[ , 1:4]))


#Apply types "adult" or "larvae" to fish based on id number.
tst<-ifelse(grepl("^A", nrelData$id),nrelData$type<-"adult",nrelData$type<-"larvae")
nrelData$type<-tst
rm(tst)

# Add data from larv to relatedNodes so that they are available for filtering on year etc.
larvNrel<-merge(nrelData,larv, by.x="id", by.y="LarvalRecords_LarvaID")
larvNrel<-larvNrel[,c(1:4,66)]
larvNrel$label<-NULL
# #Filter out unrelateds based on weights from related. Actually I need to delete them out - making a weight =0 does not prevent the edges appearing on the network
# tst<-ifelse(lrelData$weight<=0.7,0,lrelData$weight)
# lrelData$weight<-tst
# rm(tst)
#Actually can do this with delete_edges (see below)


# Collapse multiple links of the same type between the same two nodes
# by summing their weights, using aggregate() by "from", "to", & "type":
# (we don't use "simplify()" here so as not to collapse different link types)
# lrelData <- aggregate(lrelData[,3], lrelData[,-3], sum)
# lrelData <- lrelData[order(lrelData$from, lrelData$to),]
# colnames(lrelData)[4] <- "weight"
# rownames(lrelData) <- NULL
# hist(lrelData$weight)
#this is not required I think.

net<-graph_from_data_frame(d=lrelData, vertices=nrelData, directed=FALSE)

# Generate colors for adults and larvae:
# colrs <- c("tomato", "gold")
# V(net)$color <- colrs[V(net)$type]

# Set node size based on audience size:
#V(net)$size <- V(net)$audience.size*0.7

# The labels are currently node IDs.
# Setting them to NA will render no labels:
#V(net)$label.color <- "black"
#V(net)$label <- NA

# Set edge width based on weight:
#E(net)$width <- E(net)$weight*3

#change arrow size and edge color:
# E(net)$arrow.size <- .2
# E(net)$edge.color <- "gray80"

#make more sparse.for full sibs
net.FS <- delete_edges(net, E(net)[weight<0.4])
l <- layout_with_fr(net.FS)

#to colour larvae
plot(net.FS, edge.arrow.size=0, edge.curved=0, vertex.size=5,vertex.color=c("dark red", "slategrey")[(V(net.FS)$type=="larvae")+1], vertex.frame.color="#555555",vertex.label=V(net)$name, vertex.label.color="black",vertex.label.cex=.7,layout=l, main = main)

#to colour by years (need larvNrel above)
V(net.FS)$color=V(net.FS)$YearOnly #assign the "YearOnly" attribute as the vertex color
V(net.FS)$color=gsub("2011","indianred",V(net.FS)$color) #2011 will be red
V(net.FS)$color=gsub("2012","lightgoldenrod1",V(net.FS)$color) #2012 will be blue
V(net.FS)$color=gsub("2013","lightgreen",V(net.FS)$color) #2013 will be blue

plot(net.FS, edge.arrow.size=0, edge.curved=0, vertex.size=5, vertex.color=V(net.FS)$color,vertex.frame.color="#555555",vertex.label=V(net)$name, vertex.label.color="black",vertex.label.cex=.7,layout=l, main = main)

# #make more sparse.for half sibs
# net.HS <- delete_edges(net, E(net)[weight>0.3|weight<0.2])
# l <- layout_with_fr(net.FS)
# 
# plot(net.HS, edge.arrow.size=0, edge.curved=0, vertex.size=5,vertex.color=c("dark red", "slategrey")[(V(net.HS)$type=="larvae")+1], vertex.frame.color="#555555",vertex.label=V(net)$name, vertex.label.color="black",vertex.label.cex=.7,layout=l, main = "Half Siblings")
# 
# #make more sparse.for unrelateds
# net.UR <- delete_edges(net, E(net)[weight>0.15])
# l <- layout_with_fr(net.FS)
# 
# plot(net.HS, edge.arrow.size=0, edge.curved=0, vertex.size=5,vertex.color=c("dark red", "slategrey")[(V(net.HS)$type=="larvae")+1], vertex.frame.color="#555555",vertex.label=V(net)$name, vertex.label.color="black",vertex.label.cex=.7,layout=l, main = "Unrelated")

net.copy <- delete_edges(net, which(E(net)$weight <0.2))
net.copy <- delete.vertices(net.copy,which(degree(net.copy)<1))
l <- layout_with_fr(net.FS)#layout_as_tree(net.copy)
plot(net.copy, edge.arrow.size=0, edge.curved=0, vertex.size=15,vertex.color=c("dark red", "slategrey")[(V(net.copy)$type=="larvae")+1], vertex.frame.color="#555555",vertex.label=V(net.copy)$name, vertex.label.color="black",vertex.label.cex=.7, main = "Full Siblings", frame=T, vertex.label.degree=pi/2)

#edge labels and colours example below:
net.HS <- delete_edges(net, E(net)[weight>0.4])
net.HS <- delete_edges(net.HS, E(net.HS)[weight<0.22])

E(net.HS)[ weight > .32 ]$color <- "red"
E(net.HS)[ weight < .32 ]$color <- "yellow"
E(net.HS)[ weight < .25 ]$color <- "green"

l <- layout_with_fr(net.HS)

plot(net.HS, edge.arrow.size=0, edge.curved=0.2, vertex.size=5,vertex.color=c("orange", "cyan")[(V(net.HS)$type=="larvae")+1], vertex.frame.color="#555555",vertex.label=V(net)$name, vertex.label.color="black",vertex.label.cex=.7,layout=l, main = main, edge.label = E(net.HS)$weight,edge.label.cex=.5)
