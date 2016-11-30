## @knitr related

# This script allows running of 'related' package from a two row DArT file. This is an R package that allows users to estimate pairwise relatedness for pairs of individuals based on codominant molecular markers (microsatellites, SNPs, etc.), and also has simulation capabilities for comparing the performance of different estimators and for testing the resolution of a data set. It can also test for relatedness patterns within groups. Relatedness can be estimated using any of seven different methods, and can incorporate inbreeding and genotyping errors. The underlying code for implementing these relatedness estimators is Jinliang Wang’s Fortran code for his COANCESTRY program (Wang 2011).
# Alan Couch, Institute for Applied Ecology, University of Canberra

####################
# Note that one needs to run 32 bit version of R
# (set this in Global options under Tools menu)
# because 'related' relies on some Fortran code that does not run under 64 bit R.
####################

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
fileName=paste0("outData/",main," outRelated.txt",sep=" ")
write.table(df2[,1:800],file=fileName, sep=" ", col.names = FALSE, row.names = TRUE ,quote = FALSE) #(df2[,1:600], to limit number for testing)

#From here down is the 'related' tutorial code essentially. 
require(related)
input <- readgenotypedata(fileName)
fileName=paste0("outData/",main," famSimInput",sep=" ")
save(input,file = fileName)

# See COANCESTRY paper itself (Wang 2011). Use method=2 for 95%CI

#Likelihood Estimators
#output <- coancestry (input$gdata, wang=1)
output <- coancestry (input$gdata, allow.inbreeding=TRUE, trioml=1, wang=1, lynchli=1, lynchrd=1, ritland=1, quellergt=1, dyadml=1)
#or to include the ML methods and therefore allow for inbreeding:
#output <- coancestry (input$gdata, allow.inbreeding=TRUE, wang=1, dyadml = 1, trioml = 1)


fileName=paste("./outData/",main," coancestoryOutput",sep="", collapse=" ")
save(output,file = fileName)# save net files for next chunk on.
#save(output,file = paste(c("outData/",main,"coancestoryOutput"), collapse=" "))# save these files for on.


fileName=paste0("outData/",main," coancestory.csv",sep=" ")
write.table(output$relatedness, fileName, sep=",", col.names = TRUE, row.names = TRUE) #capture related data in files,for the use of.

