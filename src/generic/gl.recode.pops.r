# An R function to recode populations in a DArT genlight SNP file
#
# This program recodes population assignments and/or deletes populations from a DaRT genlight SNP file
# The program also eliminates monomoprhic loci and loci with all missing data that arises when populations are deleted.
#
# INPUT
# x -- name of an input genlight object containing the DArT data [Required]
# pop.recode -- name of a recode table file, in csv, to recode populations (including "Delete") [Default NULL]
#
# OUTPUT
# a genlight file with the recoded and reduced data -- returned
#
# USEAGE
# gl <- gl.recode.pops(x=gl, pop.recode="recode_table.csv")
#
# Author: Arthur Georges
# Last update: 31-Jan-16
#
# Dependencies

library(plyr)
library(adegenet)
library(tidyr)
library(reshape2)

gl.recode.pops <- function(x, pop.recode){
  
#  x=gl
#  pop.recode="pop_recode_table_0.csv"

# RECODE POPULATIONS
  cat("Reassigning entities to populations as per ", pop.recode, "\n")
  recode.table <- read.csv(pop.recode, stringsAsFactors=FALSE, header=FALSE);
# Apply the recode to the populations
  pop.list <- as.character(pop(x));
  ntr <- length(recode.table[,1])
  for (i in 1:nInd(x)) {
    for (j in 1:ntr) {
      if (pop.list[i]==recode.table[j,1]) {pop.list[i] <- recode.table[j,2]}
    }
  }
  pop(x) <- pop.list
  
# Remove rows flagged for deletion
  cat("Removing entities flagged for deletion in ", pop.recode, "\n")
  x <- x[!x$pop=="Delete"]
  
# REPORT A SUMMARY
    cat ("\n\nSummary of Genlight object\n\n")
    cat ("------------------------------\n")
    cat ("\nNo. of entities (individuals) =", nInd(x), "[count can be accessed using nInd(x)]\n")
    cat("  aggregated into ", nPop(x), "populations [count can be accessed using nPop(x)]\n")
    cat ("No. of attributes (loci) =", nLoc(x), "[count can be accessed using nLoc(x)]\n")
    cat ("Data can be accessed using as.matrix(x)\n\n")
    cat("Names of entities (individuals) can be accessed with indNames(x)\n")
    cat("Names of attributes (loci) can be accessed with locNames(x)\n")
    cat("Names of populations can be accessed with popNames(x)\n\n")
    cat("A list of other associated data can be obtained with names(x\b$other)\n")
    cat ('  which can then be accessed with statements of the form x@metrics["clid"]')

    return <- x

}

