#A trial to run outbreds seperately - perhaps to be incorporated into munge and pipeline.

library(dplyr)
require(related)
load("outData/PICorderJune2016/2011-2013 Larval snps coancestoryOutput")
outbreds<-output$inbreeding
obs<-filter(outbreds,LH< (-0.4))
nobs<-filter(outbreds,LH>= (-0.4))
larvalPeeliSnpsTmp<-data.frame(t(larvalPeeliSnps))

larvalPeeliSnpsOutbreds<-larvalPeeliSnpsTmp[which(larvalPeeliSnpsTmp$X1%in% obs$ind.id),]
larvalPeeliSnpsOutbredsNOT<-larvalPeeliSnpsTmp[which(larvalPeeliSnpsTmp$X1%in% nobs$ind.id),]

larvalPeeliSnpsOutbreds<-data.frame(t(larvalPeeliSnpsOutbreds))
larvalPeeliSnpsOutbredsNOT<-data.frame(t(larvalPeeliSnpsOutbredsNOT))



df<-larvalPeeliSnpsOutbreds; main<-"2011-2013 Larval snps outbreds only"
df<-df[-c(1:4),]
pic<-df[2:nrow(df),14] #check pic in order
plot(as.character(pic),main="Polymorphic Information Content (PIC)", ylab="PIC")# to show PIC is in order (ish) ideally cut off would be 1300 or so.
rm(pic)

#But choose only first 800 so PC does not crash
df<-df[-c(802:nrow(df)),]

#For each row create a unique alleleID (I wanted a shorther unique allele name - there is about 12000 of them).
require(digest)
for (i in 1:length(df$X..2))   {     
        df$unique[i]<-digest(df$X..2[i],algo="xxhash32",seed=1)
}

df<-cbind(df$unique, df)
names(df)[names(df)=="df$unique"] <- "allele"
df$unique<-NULL
df[1:10,c(1,4)]

df2<-df
#remove unwanted (all) meta columns
df2[2:18]<-list(NULL)
#transpose
df2<-setNames(data.frame(t(df2[,-1])), df2[,1]) #and make colnames from new alleleID
rownames(df2)<-df2[,1]    # sensible row names
fishID<-df2[,1] #keep names for later

df2[,1]<-NULL
colnames(df2)<-paste("X",colnames(df2),sep="") #Can't have numbers for colnames in R 
df2<-as.data.frame(apply(df2,2,function(x)gsub('\\s+', '',x)))#remove white space etc. # Maybe this should be right up front. And in other scripts.

# Make numeric as currently factor and related requires integers
indx <- sapply(df2, is.factor)
df2[indx] <- lapply(df2[indx], function(x) as.numeric(as.character(x)))
rm(indx)
# replace each value 1 with the column number (thus making unique Allele ID number)
x<-df2
library(data.table)
setDT(x)
for(j in seq_along(x)){
        set(x, i= which(x[[j]]==1 & !is.na(x[[j]])), 
            j=j, value= j)
}

d1<-data.frame(x)
rm(x);rm(i);rm(j)
#############################
#If value 0 in second column of pair change it to the integer in the first column of the pair.# If value 0 in first column of pair change it to the integer in the second column of the pair
d1.lst <- lapply(seq(1, ncol(d1), by=2), function(x){return(d1[, x:(x+1)])})
d1.lst.fil <- lapply(d1.lst, function(x){
        x[,1][!is.na(x[,1]) & x[,1]=="0"] <- x[,2][!is.na(x[,1]) & x[,1]=="0"]
        x[,2][!is.na(x[,2]) & x[,2]=="0"] <- x[,1][!is.na(x[,2]) & x[,2]=="0"]
        return(x)
})
df2 <- do.call(cbind, d1.lst.fil)
rm(d1);rm(d1.lst);rm(d1.lst.fil)
#############################

#add column with fish IDs before writing
df2<-cbind(fishID,df2)
rownames(df2)<-df2[,1]    # sensible row names again
df2[,1]<-NULL

# Write out and then read in. 
#why? Because as noted on page 7 of https://frasierlab.files.wordpress.com/2015/03/tutorial1.pdf that with this method we get all of the other associated data, such as number of loci, frequencies, etc.) 
#Also smaller sample just so 'related::sim' does not crash. Ok with 600 aleles - crashed with more doing (100 sims)
#write.table(df2[,1:800],"outRelated.txt", sep=" ", col.names = FALSE, row.names = TRUE ,quote = FALSE) #(df2[,1:600], to limit number for testing)
fileName=paste0("outData/",main," outRelatedOutbreds.txt",sep=" ")
write.table(df2[,1:800],file=fileName, sep=" ", col.names = FALSE, row.names = TRUE ,quote = FALSE) #(df2[,1:600], to limit number for testing)

#From here down is the 'related' tutorial code essentially. 
require(related)
input <- readgenotypedata(fileName)
fileName=paste0("outData/",main," famSimInputOutbreds",sep=" ")
save(input,file = fileName)

# See COANCESTRY paper itself (Wang 2011). Use method=2 for 95%CI

#Likelihood Estimators
#output <- coancestry (input$gdata, wang=1)
output <- coancestry (input$gdata, trioml=1, wang=1, lynchli=1, lynchrd=1, ritland=1, quellergt=1, dyadml=1)
#or to include the ML methods and therefore allow for inbreeding:
#output <- coancestry (input$gdata, allow.inbreeding=TRUE, wang=1, dyadml = 1, trioml = 1)


fileName=paste("./outData/",main," coancestoryOutputOutbreds",sep="", collapse=" ")
save(output,file = fileName)# save net files for next chunk on.
#save(output,file = paste(c("outData/",main,"coancestoryOutput"), collapse=" "))# save these files for on.


# fileName=paste0("outData/",main," coancestoryOutbreds.csv",sep=" ")
# write.table(output$relatedness, fileName, sep=",", col.names = TRUE, row.names = TRUE) #capture related data in files,for the use of.

#Simulation code direct from 'related' tutorial paper
fileName=paste0("./outData/",main," famSimInputOutbreds",sep=" ")
load(file = fileName)
sim <- familysim(input$freqs,100)
outputSim <- coancestry ( sim , wang=1)
##output2 <- coancestry ( sim , trioml=1) # takes 22 hours so make an export for Bernd
simrel <- cleanuprvals (outputSim$relatedness , 100)
fileName=paste0("outData/",main," outputSimOutbreds",sep=" ")
save(simrel,file = fileName)
relvalues <- simrel [, 6] #5 for trioml, 6 for wang, 7 for lynchli,8 lynchrd, 9 ritland,10quellergt, 11 dyadml
label1 <- rep ("PO", 100)
label2 <- rep (" Full ", 100)
label3 <- rep (" Half ", 100)
label4 <- rep (" Unrelated ", 100)
labels <- c( label1 , label2 , label3 , label4 )


#Plot as a density plot
Relationship <- labels
newdata <- as.data.frame ( cbind ( Relationship , relvalues ) )
newdata$relvalues <- as.numeric ( as.character ( newdata$relvalues ))
qplot ( relvalues , ..density.. , data = newdata , geom ="density", colour = as.factor ( Relationship ) , xlab =" Relatedness Value ", ylab ="Density",main=main)