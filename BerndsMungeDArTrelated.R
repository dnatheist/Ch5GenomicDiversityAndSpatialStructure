#Bernds Munge for DArT to 'related'. Not working yet. One problem is that adegent for 32 bit not installed.

source("src/generic/dart.r")
source("src/generic/read.dart.r")

#library(devtools)
#dev_mode(on=T)
#dev_mode(on=F)

#Install non- cran version (the older version of adgenet has a serious flaw and is not updated yet.)
#install_github("thibautjombart/adegenet")
#install_github("green-striped-gecko/PopGenReport")


all.dart <- read.dart("OtherData/allDArTsnps.csv", topskip = 5)
gl2 <- dart2genlight(all.dart, covfilename = "OtherData/qslDartCovariatesAll.csv")



gi <- genlight2genind(gl2)

glibrary(pegas)

library(pegas)

geno <-data.frame(as.loci(gi))

genotype <- matrix(-1, nrow=nrow(geno), ncol=2*ncol(geno))

for (i in 1:ncol(geno))
{
        
        for (ii in 1:nrow(geno))
        {
                
                if (is.na(geno[ii,i]))   { genotype[ii,(i-1)*2+1] <- 0;   genotype[ii,(i-1)*2+2] <- 0}   else {
                        if (geno[ii,i] == "A/A") { genotype[ii,(i-1)*2+1] <- 1;  genotype[ii,(i-1)*2+2] <- 1}
                        if (geno[ii,i] == "A/B") { genotype[ii,(i-1)*2+1] <- 1;  genotype[ii,(i-1)*2+2] <- 2}
                        if (geno[ii,i] == "B/B") { genotype[ii,(i-1)*2+1] <- 2;  genotype[ii,(i-1)*2+2] <- 2}
                        
                }
        }
        
        
}


genotype <- cbind(ind=rownames(geno), as.data.frame(genotype))

genotype$ind <- as.character(genotype$ind)

require(related)

input <- readgenotypedata(genotype)
output <- coancestry(input$gdata, lynchrd=2, quellergt=2, wang=2)
