#function below run first addnew
#pretty rough as at 8 July 16.

require(dplyr)
links<-output$relatedness[,2:5]
links$group<-NULL

lookup<-data.frame(lookupVariable="ind1.id",larv$LarvalRecords_LarvaID,newVariable="ind1.id",larv$LarvalRecords_NestID, source="")

lookup<-rename(lookup, lookupValue=larv.LarvalRecords_LarvaID, newValue=larv.LarvalRecords_NestID)
write.csv(lookup,"lookup.csv")

allowedVars <- c("ind1.id")
tmp1Pass<-addNewData("lookup.csv", links, allowedVars)

lookup<-data.frame(lookupVariable="ind2.id",larv$LarvalRecords_LarvaID,newVariable="ind2.id",larv$LarvalRecords_NestID, source="")

lookup<-rename(lookup, lookupValue=larv.LarvalRecords_LarvaID, newValue=larv.LarvalRecords_NestID)
write.csv(lookup,"lookup.csv")

allowedVars <- c("ind2.id")
tmp2Pass<-addNewData("lookup.csv", tmp1Pass, allowedVars)

#Make Net part
require(igraph)
#to see how the data is distributed
#hist(tmp2Pass[,3], main="Histogram of estimator") 

#Makes links DF
colnames(tmp2Pass)[1] <- "from"
colnames(tmp2Pass)[2] <- "to"
colnames(tmp2Pass)[3] <- "weight"
tmp2Pass$type<-"probRel"

tmp2PassNoRows<-nrow(tmp2Pass) #get number of rows and store for later use.
#create new DF for nodes
tmp2PassData <- data.frame(matrix(ncol = 4, nrow = tmp2PassNoRows))
vertData<-tmp2PassData 
colnames(vertData) <- c("id", "name", "type", "label")

#Populate the DF with the 'to' nodes from original list
vertData[,1]<-tmp2Pass[,2]
vertData[,2]<-tmp2Pass[,2]

#add the first entry from the 'from' list (otherwise it wont get in as it is not
# in the 'to list')
newrow = c(tmp2Pass[1,1],tmp2Pass[1,1],NA,NA)
vertData = rbind(vertData,newrow)
vertData$type<-"nest"
vertData$label<-vertData$name
# nrow(vertData); length(unique(vertData$id))
# nrow(lrelData); nrow(unique(lrelData[,c("from", "to")]))

vertData<-data.frame(unique(vertData[ , 1:4]))


# # Add data from larv to relatedNodes so that they are available for filtering on year etc.
# larvNrel<-merge(vertData,larv, by.x="id", by.y="LarvalRecords_LarvaID")
# larvNrel<-larvNrel[,c(1:4,66)]

#Have a look at two data frames used to make net
# head(vertData)
# head(larvNrel)

#Create the network data
net<-graph_from_data_frame(d=tmp2Pass, vertices=vertData, directed=FALSE)
l <- layout_with_fr(net)

net <- delete_edges(net, E(net)[weight<0.4])

E(net)[ weight > .34 ]$color <- "darkgreen"
E(net)[ weight < .3 ]$color <- "green"
E(net)[ weight < .24 ]$color <- "yellow"

net<-simplify(net, remove.loops = TRUE)

plot(net, edge.arrow.size=0, edge.curved=0, vertex.size=5, vertex.color=V(net)$color,vertex.frame.color="#555555",vertex.label=V(net)$name, vertex.label.color="black",vertex.label.cex=.7,layout=l, main = "Nests Relatedness")






##' Modifies 'data' by adding new values supplied in newDataFileName
##'
##' newDataFileName is expected to have columns 
##' c(lookupVariable,lookupValue,newVariable,newValue,source)
##' 
##' Within the column 'newVariable', replace values that
##' match 'lookupValue' within column 'lookupVariable' with the value
##' newValue'.  If 'lookupVariable' is NA, then replace *all* elements
##' of 'newVariable' with the value 'newValue'.
##'
##' Note that lookupVariable can be the same as newVariable.
##'
##' @param newDataFileName name of lookup table
##' @param data existing data.frame
##' @param allowedVars vector of permissible variable names for newVariable
##' @return modified data.frame
addNewData <- function(newDataFileName, data, allowedVars){
        
        import <- readNewData(newDataFileName, allowedVars)
        
        if( !is.null(import)){    
                for(i in seq_len(nrow(import))){  #Make replacements
                        col.to <- import$newVariable[i] 
                        col.from <- import$lookupVariable[i]
                        if(is.na(col.from)){ # apply to whole column
                                data[col.to] <- import$newValue[i]
                        } else { # apply to subset
                                rows <- data[[col.from]] == import$lookupValue[i]
                                data[rows,col.to] <- import$newValue[i]
                        }
                }   
        }      
        data
}

##' Utility function to read/process newDataFileName for addNewData
##' 
##' @param newDataFileName name of lookup table
##' @param allowedVars vector of permissible variable names for newVariable
##' @return data.frame with columns c(lookupVariable,lookupValue,newVariable,newValue,source)
readNewData <- function(newDataFileName, allowedVars){
        
        if( file.exists(newDataFileName)){
                import <- read.csv(newDataFileName, header=TRUE, stringsAsFactors=FALSE,
                                   strip.white=TRUE)
                if( nrow(import)> 0 ){
                        
                        #Check columns names for import are right
                        expectedColumns<- c("lookupVariable","lookupValue","newVariable","newValue")
                        nameIsOK <-  expectedColumns %in% names(import)
                        if(any(!nameIsOK))
                                stop("Incorrect name in lookup table for ",
                                     newDataFileName, "--> ", paste(expectedColumns[!nameIsOK],
                                                                    collapse=", "))
                        
                        #Check values of newVariable are in list of allowed variables
                        import$lookupVariable[import$lookupVariable == ""] <- NA
                        nameIsOK <- import$newVariable %in% allowedVars
                        if(any(!nameIsOK))
                                stop("Incorrect name(s) in newVariable column of ",
                                     newDataFileName, "--> ", paste(import$newVariable[!nameIsOK],
                                                                    collapse=", "))
                } else {
                        import <- NULL
                }
        } else {
                import <- NULL
        }
        import
}