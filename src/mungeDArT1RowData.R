#This file is to load 1 row data and get an output that has row and col names sensible so it can be used for 'public' data and iterating correlation through each snp. 

# Having trouble with * column name and extracting unwanted fist at the second. 18/2/2016

library(dplyr)
##Load DArT 1 row snps data frame
df<-read.csv("./data/Report-DMac15-1861SNP1Row.csv")

df$X.<-sub("\\|.*", "", df$X.)

#Create Unique identifier
df<-transform(df,id=paste((rownames(df)),df$X.,sep=""))

rownames(df)<-df$id
df<-df %>% select(id,everything())
#df$id<-NULL

df<-setNames(data.frame(t(df[,-1])), df[,1])# Just transpose this way larvae become observations and snps become the variables.

## Create a unique identifier as DArT did not import comments with larval fish ID
df$X1<-paste(df[,2],df[,3],df[,4],sep="")

vec<-df$X1# List of DArT 'larval identifiers' for use recoding names with recoderFunc which translates the DaRT identifiers into something more meaningful.

recoderFunc <- function(data, oldvalue, newvalue) {
        
        # convert any factors to characters
        
        if (is.factor(data))     data     <- as.character(data)
        if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
        if (is.factor(newvalue)) newvalue <- as.character(newvalue)
        
        # create the return vector
        
        newvec <- data
        
        # put recoded values into the correct position in the return vector
        
        for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
        
        newvec
        
}

vec<-recoderFunc(vec, allIdentifiers$identifier, allIdentifiers$Comments)

df[,1]<-vec # Load more meaningful names (Labels) in place of Larval ID numbers
rm(vec)# Tidy up and now new df data frame now has correct fish labels 

#These two lines are in because long names break NewHybrids.
df<- as.data.frame(sapply(df,gsub,pattern="corin-trial1",replacement="Bend"))
df<- as.data.frame(sapply(df,gsub,pattern="narranderah-trial1",replacement="Narra"))


df<-as.data.frame(t(df))
write.csv(df,".\\otherData\\oneRowDArTsnps.csv", row.names=FALSE)

###OK to here

# ##Now to delete data not required for distance matrix creation.
# df<-df[-c(1:18),] # 17 rows of DArT calculations and allele sequences
# df<-df[,-c(2:4)] # 4 columns of DArT descriptors
# rownames(df) <- df[,1] #Make rownames meaningful
# 
# colnames(df)[1] <- "newname1"
# colnames(df)[2] <- "newname2"

##Create List with Trial Fish Samples
trialFish<-c(subset(df, grepl("trial",newname2),newname1))
trialFish<-trialFish$newname1
#trialFish<-c(trialFish,"221","222")

##Remove degraded and unwanted trial samples from data frame
degradedFish<-c("243","248","249","251","252","253")
df<-df[-which(rownames(df)%in% degradedFish),] #leaving 285
df<-df[-which(rownames(df)%in% trialFish),] #leaving 277 all good SNPs (MC and TC)

##Create list with Adult Fish identifiers
adultFish<-c(subset(df, grepl("A", newname1),newname1))
adultFish<-adultFish$X1

df$newname1 <- NULL #Remove as now redundant.
df$newname2 <- NULL #Remove as now redundant.
df$X1 <- NULL #Remove as now redundant.

######################################
# From here to a new r script
##################################

dfm<-merge(df,larv, by.x = "row.names", by.y = "LarvalID")  #The problem with this is that it drops all the adults that were in df. oh well. 250 instead of 277 (merely because they are not in larv.)

##junk
library(corrplot)
tst<-dfm[,2:6365]
indx <- sapply(tst, is.factor)
tst[indx] <- lapply(tst[indx], function(x) as.numeric(as.character(x)))

m<-cor(tst, use="pairwise.complete.obs")


corrplot(m, type="upper", order="hclust")

library(Hmisc)
rcorr(dfm$"66661459",dfm$Distance.to.Angle.Crossing..m.)

tst$`606651772`


n<-cor(dfm[,2:6365], use="pairwise.complete.obs")

library(dplyr)
n<-m[m>0.39]
n<-n[n<1]

