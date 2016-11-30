gl.filter.monomorphs <- function (x) {

# An R function to remove monomorphic loci from an agegenet genlight object (gl)
#
# When a population is deleted (using gl.recode.r), loci for remaining populations may be monomorphic, or may comprise all na values.
# This program is used in conjunction with gl.recode. 
#
# INPUT
# x -- name of a genlight object containing SNP data from DArT [Required]
# OUTPUT
# a genlight file with the polymorphic data -- returned
#
# USEAGE
# gl <- gl.remove.monomorphs(gl)
#
# Author: Arthur Georges
# Updated 31-Jan-16; 25-Feb-2016
#
  cat("Identifying monomorphic loci\nGo for a coffee .....\n")
# Create a vector to hold test results
  a <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {a[i] <- NA}
# Set up the progress counter
  pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  getTxtProgressBar(pb)
# Identify polymorphic, monomorphic and 'all na' loci
  # TRUE if monomorphic
  xmat <-as.matrix(x)    
  for (i in (1:nLoc(x))) {
    if(any(!is.na(xmat))) {
      a[i] <- (all(!xmat[,i],rm.na=TRUE) || all(!xmat[,i],rm.na=TRUE))
    }
    setTxtProgressBar(pb, i/nLoc(x))  
  }
# Count the number of monomorphic loci (TRUE), polymorphic loci (FALSE) and loci with no scores (all.na) 
  counts <- count(a)
  cat("\nPolymorphic loci:", counts[1,2], "\nMonomorphic loci:", counts[2,2], "\nLoci with no scores (all NA):" , counts[3,2] ,"\n")
#Treat all na loci as monomorphic
  # TRUE if monomorphic or all na
  a[is.na(a)] <- TRUE 
# Write the polymorphic loci to a new genlight object
  cat("Deleting monomorphic loci and loci with no scores\n")
  x <- x[,(a==FALSE)]
  x@other$metrics <- x@other$metrics[(a==FALSE),]  

return <- x

}

