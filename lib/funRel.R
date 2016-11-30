funRelated<-function (df){
#For each row create a unique alleleID (maybe redundant but I wanted a shorther unique allele name (about 12000 of them)).
require(digest)
for (i in 1:length(df$X..2))   {     
        df$unique[i]<-digest(df$X..2[i],algo="xxhash32",seed=1)
}

df<-cbind(df$unique, df)
names(df)[names(df)=="df$unique"] <- "allele"
df$unique<-NULL
df[1:10,1:10]

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

#Write out and then read in. why? Because as noted on page 7 of https://frasierlab.files.wordpress.com/2015/03/tutorial1.pdf that with this method we get all of the other associated data, such as number of loci, frequencies, etc.) 
#smaller sample just so 'related::sim' does not crash. Ok with 600 aleles - crashed with  (100 sims)
write.table(df2[,1:600],"outRelated.txt", sep=" ", col.names = FALSE, row.names = TRUE ,quote = FALSE) #(df2[,1:600], to limit number for testing)


#df[, sample(ncol(df), 60)] is to randomly select 60 cols. Might want to do this rather than the first 60 but would need to do it in pairs of columns?

#From here down is the 'related' tutorial code essentially. 
require(related)
input <- readgenotypedata ("outRelated.txt")

input$nloci #An integer containing the number of loci used.
input$nalleles #A series of integers specifying the number of alleles at each locus.
input$ninds #An integer containing the number of individuals in the genotype file.

# See COANCESTRY paper itself (Wang 2011). Use method=2 for 95%CI

#Likelihood Estimators
output <- coancestry (input$gdata, allow.inbreeding=TRUE, dyadml = 2)

head(output$relatedness)# A data frame containing all pairwise estimates of relatedness. Thiswill always have 11 columns: (1) an integer for the pair number; (2) the ID for individual #1; (3) the ID for individual #2; (4) the group assignment (see section 4.5); (5 - 11) for the 7 relatedness estimators—contain values of 0 for estimators not chosen.
#head(output$delta7) # A data frame that contains the ∆7 estimates for the relatedness estimators that use it (trioml, wang, lynchrd, dyadml). This data frame contains one row for each pairwise comparison, and 8 columns: (1) an integer for the pair number; (2) the ID for individual #1; (3) the ID for individuals #2; (4) the group assignment (See section 4.5); and (5 - 8) estimates of ∆7 for the 4 relevant estimators (trioml, wang, lynchrd, dyadml), with values of 0 for estimators not chosen.
#head(output$delta8) # A data frame that contains the ∆8 estimates for the relatedness estimators that use it (trioml, wang, lynchrd, dyadml). This data frame contains one row for each pairwise comparison, and 8 columns: (1) an integer for the pair number; (2) the ID for individual #1; (3) the ID for individuals #2; (4) the group assignment (See section 4.5); and (5 - 8) estimates of ∆7 for the 4 relevant estimators (trioml, wang, lynchrd, dyadml), with values of 0 for estimators not chosen.
head(output$inbreeding) # A data frame that contains the inbreeding estimates for each individual, as used in the relatedness estimators. Only four of the relatedness estimators can account for inbreeding: dyadml, lynchrd, ritland, trioml. This data frame contains one row for each individual, and 5 columns: (1) individual ID; (2-5) inbreeding estimates for the four relatedness estimators. Estimators not used will have a zero (0) in the corresponding column.


# Simulation code direct from 'related' tutorial paper
sim <- familysim(input$freqs,100)
output2 <- coancestry (sim, wang =1)
simrel <- cleanuprvals ( output2$relatedness , 100)
relvalues <- simrel [, 6]
label1 <- rep ("PO", 100)
label2 <- rep (" Full ", 100)
label3 <- rep (" Half ", 100)
label4 <- rep (" Unrelated ", 100)
labels <- c( label1 , label2 , label3 , label4 )

#Plot the Relatedness
qplot ( as.factor ( labels ) , relvalues , geom ="boxplot", ylab =" Relatedness Values ", xlab ="Relatedness ")

#Plot as a density plot
Relationship <- labels
newdata <- as.data.frame ( cbind ( Relationship , relvalues ) )
newdata$relvalues <- as.numeric ( as.character ( newdata$relvalues ))
qplot ( relvalues , ..density.. , data = newdata , geom ="density", colour = as.factor ( Relationship ) , xlab =" Relatedness Value ", ylab ="Density")

#From the simulations of 'known' relationships we can choose cut off values for identifying relationships and colouring them in the matrix. In this case FS or PO above 0.4; HS 0.1-0.4; unrelated fish have a relatedness value below 0.15.

## Compare Estimators

require(related)
compareestimators (input , 100)
}
