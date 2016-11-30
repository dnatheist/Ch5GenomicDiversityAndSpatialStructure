#Gather files and generate metaData
fishNos<-data.frame(t(larvalPeeliSnps[5,18:ncol(larvalPeeliSnps)])) #collect numbers
fishNos<-as.integer(as.character(fishNos$X5)) 
larvMeta <- larv[larv$LarvalRecords_LarvaID %in% fishNos, ]#obtain subset of meta data
larvMeta[1]<-NULL #remove label column up front. Not sure if this neccessary
larvMeta<-larvMeta[,c(1,62,65,96,97)] #choose metadata to include
colnames(larvMeta)[1]<-"id"
colnames(larvMeta)[3]<-"pop"
colnames(larvMeta)[4]<-"lat"
colnames(larvMeta)[5]<-"lon"
larvMeta<-larvMeta[match(fishNos,larvMeta$id),] #order meta data file (or creating gl will break)
#Save MetaData file to read in
write.csv(larvMeta,"otherData/larvMeta.csv", row.names=FALSE) #create metadata file

#Create genlight object #64 bit R only
gl<-gl.read.dart(datafile="otherData/larvalPeeliSnps.csv", topskip=5, nmetavar=17, nas="-", ind.metafile="otherData/larvMeta.csv")

# setdiff(larvMeta$id,fishNos)
# setdiff(fishNos,larvMeta$id)

### the following will be best in seperate file after restartign rstudio with R 32 bit
input.df <- gl2related(gl)
input<-readgenotypedata(input.df)#this adds a few counts and stuff even though we have already read in the file

require(related) #32 bit R only.
output <- coancestry (input$gdata, wang=1, dyadml = 1, trioml = 1)
