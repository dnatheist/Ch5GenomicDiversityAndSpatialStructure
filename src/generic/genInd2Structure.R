# Function to export to STRUCTURE format from genind object.
# genind objects are created in the R package adegenet.  The function below is an R function.
# Lindsay V. Clark, 14 March 2013

# obj: genind object
# file: file name to write
# pops: whether to include population info in the file
# Function is flexible with regards to ploidy, although genotypes are
# considered to be unambiguous.
# Missing data must be recorded as NA in obj@tab.

# example use: 
# data(nancycats)
# genind2structure(nancycats, file="nancy_structure.txt", pops=TRUE)
genind2structure(gi.dart, file="structure.txt", pops=FALSE)

genind2structure <- function(obj, file="", pops=FALSE){
        if(!"genind" %in% class(obj)){
                warning("Function was designed for genind objects.")
        }
        
        # get the ploidy of the dataset
        pl <- obj@ploidy
        # get the number of individuals
        S <- length(obj@ind.names)
        # column of individual names to write; set up data.frame
        tab <- data.frame(ind=rep(obj@ind.names, each=pl))
        # column of pop ids to write
        if(pops){
                popnums <- 1:length(obj@pop.names)
                names(popnums) <- obj@pop.names
                popcol <- rep(popnums[as.character(pop(obj))], each=pl)
                tab <- cbind(tab, data.frame(pop=popcol))
        }
        
        # begin going through loci
        gen <- obj@tab * pl # genotypes with allele counts
        loci <- names(obj@loc.names) # use the simplified names here, convert later
        for(L in loci){
                thesegen <- gen[,grep(L, dimnames(gen)[[2]])] # genotypes by locus
                al <- 1:dim(thesegen)[2] # numbered alleles
                loccol <- c() # column to add to table
                for(s in 1:S){
                        if(!all(!is.na(thesegen[s,]))){
                                loccol <- c(loccol, rep(-9, pl))
                        } else {
                                loccol <- c(loccol, rep(al, times=thesegen[s,]))
                        }
                }
                if(length(loccol) != S*pl) stop("Missing data must be recorded as NA.")
                tab <- cbind(tab, data.frame(x=loccol))
                names(tab)[length(names(tab))] <- obj@loc.names[L]
        }
        
        # export table
        write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
}
