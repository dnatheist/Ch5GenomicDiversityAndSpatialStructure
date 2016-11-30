#Prepare data
prelDataA<-cbind(larvNrel[,c(1,13)],rep("1",nrow(larvNrel)),rep("mother",nrow(larvNrel)))

colnames(prelDataA)[1]<-"from"
colnames(prelDataA)[2]<-"to"
colnames(prelDataA)[3]<-"weight"
colnames(prelDataA)[4]<-"type"

prelDataB<-cbind(larvNrel[,c(1,14)],rep("1",nrow(larvNrel)),rep("father",nrow(larvNrel)))
colnames(prelDataB)[1]<-"from"
colnames(prelDataB)[2]<-"to"
colnames(prelDataB)[3]<-"weight"
colnames(prelDataB)[4]<-"type"
prelData<-rbind(prelDataA,prelDataB)
prelData[prelData==""] <- NA
rm(prelDataB);rm(prelDataA)
tmp<-prelData[complete.cases(prelData),]

#Create Vertices Data Frame
prelVert <- data.frame(matrix(ncol = 4, nrow = nrow(tmp)))
colnames(prelVert) <- c("id", "name", "type", "label")

prelVert$id<-tmp$to
prelVert$name<-tmp$to
prelVert$type<-tmp$type
prelVert<-rbind(prelVert, nrelData)
prelVert<-unique(prelVert)


prelData<-prelData[complete.cases(prelData),]

prelData <- subset(prelData, !from == 177) #these line needed to remove a few potential contaminants
prelData <- subset(prelData, !from == 314) #although I doubt there was actually contam as the otherhs are the same.
prelData <- subset(prelData, !from == 366)


require(igraph)
net.parents<-graph_from_data_frame(d=prelData, vertices=prelVert, directed=FALSE)

l <- layout_with_fr(net.parents)

#Colour vertices :parents and larvae
V(net.parents)$color=V(net.parents)$type #assign the "type" attribute as the vertex color then assign a colour based on that type.
V(net.parents)$color=gsub("mother","lightpink",V(net.parents)$color)
V(net.parents)$color=gsub("father","lightblue",V(net.parents)$color)
V(net.parents)$color=gsub("larvae","lightgreen",V(net.parents)$color)
V(net.parents)$color=gsub("adult","orange",V(net.parents)$color)

plot(net.parents, edge.arrow.size=0, edge.curved=0, vertex.size=5, vertex.color=V(net.parents)$color,vertex.frame.color="#555555",vertex.label=V(net.parents)$name, vertex.label.color="black",vertex.label.cex=.7,layout=l, main = "Full Siblings With Parents")

##Now to add in sib weights so we can include half sibs and exclude full sib edges.
prelData$weight<-as.numeric(as.character(prelData$weight))
allData<-rbind(prelData, lrelData)
allData<-allData[complete.cases(allData),]
net.all<-graph_from_data_frame(d=allData, vertices=prelVert, directed=FALSE)

#set parmeters for edges to delete (so they dont show and clutter the graph)
net.all <- delete_edges(net.all, E(net.all)[weight>0.379])
net.all <- delete_edges(net.all, E(net.all)[weight<0.2])

#to colour half sibling edges according to weight 
E(net.all)[ weight > .34 ]$color <- "darkgreen"
E(net.all)[ weight < .3 ]$color <- "green"
E(net.all)[ weight < .24 ]$color <- "yellow"

#to colour parents and larvae
V(net.all)$color=V(net.all)$type #assign the "type" attribute as the vertex color
V(net.all)$color=gsub("mother","lightpink",V(net.all)$color) #mums will be pink
V(net.all)$color=gsub("father","lightblue",V(net.all)$color) #dads will be blue
V(net.all)$color=gsub("larvae","lightgreen",V(net.all)$color) #larvae will be green
V(net.all)$color=gsub("adult","orange",V(net.all)$color) #adult will be orange

plot(net.all, edge.arrow.size=0, edge.curved=0, vertex.size=5, vertex.color=V(net.all)$color,vertex.frame.color="#555555",vertex.label=V(net.all)$name, vertex.label.color="black",vertex.label.cex=.7,layout=l, main = "Half Siblings With Parents")

#Now another look coloured by Years
net.all<-graph_from_data_frame(d=allData, vertices=prelVert, directed=FALSE)

#set parmeters for edges to delete (so they dont show and clutter the graph)
net.all <- delete_edges(net.all, E(net.all)[weight>0.379])
net.all <- delete_edges(net.all, E(net.all)[weight<0.2])

#to colour half sibling edges according to weight 
E(net.all)[ weight > .34 ]$color <- "darkgreen"
E(net.all)[ weight < .3 ]$color <- "green"
E(net.all)[ weight < .24 ]$color <- "yellow"

#to colour vertices by yearsparents and larvae
V(net.all)$color=gsub("mother","lightpink",V(net.all)$color) #mums will be pink
V(net.all)$color=gsub("father","lightblue",V(net.all)$color) #dads will be blue
V(net.all)[name<134]$color<-"red"#2011 will be red
V(net.all)[name>133]$color<-"yellow"#2012 will be yellow
V(net.all)[name>184]$color<-"lightgreen" #2013 will be green


plot(net.all, edge.arrow.size=0, edge.curved=0, vertex.size=5, vertex.color=V(net.all)$color,vertex.frame.color="#555555",vertex.label=V(net.all)$name, vertex.label.color="black",vertex.label.cex=.7,layout=l, main = "Half Siblings With Some Parents Coloured by Year")