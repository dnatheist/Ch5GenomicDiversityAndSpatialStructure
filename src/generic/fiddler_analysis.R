# fiddler_analysis.R

# Terry Bertozzi
# 16.iii.2015

# This script takes the output (SNPs_Selected_recoded.csv) of the PERL script dart2minal.pl and 
# performs data checks and corrections before running a DAPC analysis and producing the PCA plot
# and other statistical values and graphs presented in the publication:
#
# Donnellan, S.C., Foster, R., Junge, C., Huveneers, C., Kilian, A. and Bertozzi, T. (2015) 
# Fiddling with the proof: the Magpie Fiddler Ray is a colour pattern variant of the common 
# Southern Fiddler Ray (Rhinobatidae: Trygonorrhina). _Zootaxa_ (submitted)
#
# This script requires the following packages to be installed: adegenet, devtools, ggplot2.
# Package version details can be found in README.md


# ===== get, check and tidy the data =====

# load the addNewData function
library(devtools, quietly=TRUE)
source_gist("https://gist.github.com/dfalster/5589956")

# get the data
data.full <- read.csv("SNPs_Selected_recoded.csv",stringsAsFactors=TRUE)

# add grouping variables
allowedvars <- c("pop", "state")
data.full <- addNewData("pop_data.csv", data.full, allowedvars)

# reorder the columns
data.full <- data.full[,c(1,4150,4151,2:4149)]
data.full$pop <- as.factor(data.full$pop)
data.full$state <- as.factor(data.full$state)

# remove the two technical replicate magpie samples
data <- data.full[! data.full$sample %in% c("abtc127320_1","abtc127320_S2"),]

# calculate the minor allele frquency
x <- colSums(data[,4:4151], na.rm=TRUE)
y <- apply(data[,4:4151], 2, FUN=function(x)length(which(!is.na(x))))
minAlleleFreq <- x/(2*y)

# check if any loci are coded the wrong way (i.e. by major allele freq)
sum(minAlleleFreq > 0.5) #386 in this case

# recode the loci as minor alleles
for(n in 1:4148){
        if (minAlleleFreq[n] > 0.5){
                recode <- gsub(0,9,data[,n + 3])
                recode <- gsub(2,0,recode)
                data[,n + 3]  <- as.integer(gsub(9, 2, recode))
        }else{
                next #last loop doesn't seem to work properly unless this is here
        }
}

# re-calculate the minor allele frquency
x <- colSums(data[,4:4151], na.rm=TRUE)
y <- apply(data[,4:4151], 2, FUN=function(x)length(which(!is.na(x))))
minAlleleFreq <- x/(2*y)

# plot the minor allele frequencies
hist(minAlleleFreq, breaks =40,main=NULL,xlab="Minor allele frequency")

# how many loci have a frequency > 0.1?
sum(minAlleleFreq >= 0.1)

# proportion of missing loci by individual
locixind <- apply(data[,4:4151], 1, FUN=function(x)length(which(!is.na(x)))/4148*100)
missing <- data.frame("sample"=data$sample,"pcentmissloci"=(100-locixind),stringsAsFactors = FALSE)

# number of missing individuals by locus
indxloci <- apply(data[,4:4151],2, FUN=function(x)length(which(!is.na(x))))
table(indxloci)

# save the cleaned up data
write.csv(data,"fiddler_snp.csv",row.names=FALSE)


# ===== DAPC analysis =====
library(adegenet)

# get the data into a genlight object
data.gl <- new("genlight", gen=data[,4:4151], ind.names=data$sample, ploidy=2, pop=data$pop,other=as.list(data$state))

# find clusters
grps <- find.clusters(data.gl,max.n.clust=20) # first interactive run
#grps <- find.clusters(data.gl,n.pca=25,n.clust=2) # subsequent runs

# plot the assignment of pops to inferred groups
table.value(table(pop(data.gl), grps$grp), col.lab=paste("inf",1:2))


dapc1 <- dapc(data.gl,grps$grp) # first interactive run
#dapc1 <- dapc(data.gl,n.pca=25,n.da=1)  # subsequent runs

# plot desities of individuals on the discriminant function 
scatter(dapc1)

# PCA plot
library(ggplot2)                                                   
df<-data.frame(PC1=dapc1$tab[, 1], PC2=dapc1$tab[, 2], sample=data$sample, state=data$state, pop=data$pop)
df2 <- subset(df, pop == "magpie")
df1 <- subset(df, pop == "southern" | pop == "eastern" )
col_palette <- c("#AC9393", "#E69F00", "#56B4E9", "#483737") #0072B2

ggplot(df1, aes(x=PC1, y=PC2, group = state)) +
        geom_point(aes(colour = state, shape = pop, fill=state), size=3, alpha=1) + # symbol colour depends on "State"
        scale_shape_manual(values = c(24,8,21))   +  # Manually select the symbols (eastern,magpie,southern)
        scale_colour_manual(values=col_palette) + # and colours
        scale_fill_manual(values=col_palette) + #and fill
        geom_point(data=df2,aes(shape =pop),colour="black",fill="black", size=3, alpha=1) + # add the hybrid and magpie layer
        theme(legend.position="none", panel.background=element_rect(fill=NA),panel.border=element_rect(colour="black",fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +#remove the legend and format plot area
        theme(axis.ticks=element_line(colour="black"), axis.text=element_text(colour="black")) #format axis ticks and labels


# ===== Fixed locus and heterozygosity evaluation =====
# remove the putative hybrids
data.sub <- data[! data$sample %in% c("abtc89717","abtc84053"),]

# rename the magpie fiddlers
data.sub$pop[9] <- "southern"
data.sub$pop[14] <- "southern"

col.names <- colnames(data.sub)

fixed <- 0

locus.num <- integer()
locus.name <- character()
for (i in 4:4151){#4151
        if (!any(data.sub[,i] %in% 1,na.rm=TRUE)){ #check if there are any hets
                d <-mean(data.sub[which(data.sub$pop=='southern'),i], na.rm=TRUE)
                f <-mean(data.sub[which(data.sub$pop=='eastern'),i], na.rm=TRUE)
                if ((d==0 & f==2) | (d==2 & f==0)){
                        fixed <- fixed + 1
                        locus.num <- c(locus.num,i)
                        locus.name <- c(locus.name,col.names[i])
                }
        }
}
df <- data.frame(locus.num,locus.name)

het.count<-0
for (i in 4:4151){
        if (any(data.sub[,i] %in% 1,na.rm=TRUE)){ #check if there are any hets
                het.count<-het.count+1
        }    
}    

# heterozygosity per sample
het.sample <- apply(data[,4:4151],1, FUN=function(x)length(which(x==1)))
data.frame(data$sample,"het.loci" = het.sample,"percent.het.loci" = het.sample/4148*100)

# heterozygozity per locus
het.locus <- apply(data.sub[,4:4151],2, FUN=function(x)length(which(x==1)))    
hist(het.locus, breaks=33, xlab="Heterozygotes / locus")


# ===== contribution of loci to PCA =====
vc <- dapc1$var.contr
vc.sorted <-sort(vc, decreasing=TRUE)
cumm <- data.frame("proportion"=vc.sorted,"cumm.propn"=cumsum(vc.sorted))
cumm$num <- 1:4148
plot(cumm$num, cumm$cumm.propn, xlab="locus count", ylab="cummulative proportion of variance ", cex=0.5)

# ===== top PCA loadings for PC1 =====
# get PC1 loadings and the locus names
pc1.load <- data.frame("loadings"=dapc1$pca.loadings[1:4148,1])
pc1.load$locus <-names(data[4:4151])

# convert negative loadings to positive
pc1.load["loadings"]<-lapply(pc1.load["loadings"], function(x) abs(x))

# sort by loadings
pc1.load <- pc1.load[with(pc1.load, order(-loadings,locus)),]

# grab the first 200 loci
loci.200 <- pc1.load[1:200,"locus"]

# subset "data"
data1 <- data[,loci.200]


# ===== conversion of data for NewHybrids =====

# copy the locus data
data1 <-data[,4:4151]

# recode the loci as HewYhbrids lumped format
for(n in 1:ncol(data1)){
        data1[,n]<- gsub(2,22,data1[,n])
        data1[,n] <- gsub(1,12,data1[,n])
        data1[,n] <- gsub(0,11,data1[,n])
        data1[,n][is.na(data1[,n])] <- 0
}

# make numeric sample IDs 
sampleID<-seq_along(data$sample)

# write the NewHybrids header to MAC classic file
sink("fiddler_snp_newhybrid_all.dat")
cat(sprintf("NumIndivs %d\r", length(sampleID)))
cat(sprintf("NumLoci %d\r", length(names(data1))))
cat("Digits 1\r")
cat("Format Lumped\r")
cat("\r")
cat("LocusNames ")
cat(sprintf("% s", names(data1)))
cat("\r\r")
sink()


# add sample IDs to dataframe and reorder
data1$sampleID <-sampleID
data1 <-data1[,c(ncol(data1),1:(ncol(data1)-1))]

# write the genotpyes to file 
write.table(data1,"fiddler_snp_newhybrid_all.dat",append=TRUE, eol="\r", quote=FALSE, sep=" ",col.names=FALSE, row.names=FALSE)