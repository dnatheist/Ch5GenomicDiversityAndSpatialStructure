gl.filter.RepAvg <- function(x, threshold) {
  
  # Filter loci on reproducibility (RepAvg)
  
  # A function to remove loci from a genlight object based on the estimate
  # of reproducibility (=RepAvg).
  #
  # Input: A genlight object containing SNP genotypes.
  #
  # Returns a genlight object with loci with a RepAvg less
  # than the specified threshold deleted.
  #
  # Useage: result <- gl.filter.RepAvg(gl, 0.99)
  #
  # Author: Arthur Georges
  # Last update: 25-Feb-2016
  # Dependencies: library(adegenet)
  #
  # BEGIN FUNCTION
  #
    n0 <- nLoc(x)
    cat("Initial no. of loci =", n0, "\n")  
  # Remove SNP loci with RepAvg >= threshold
    x2 <- x[, x@other$metrics["RepAvg"]>=threshold]
  # Remove the corresponding records from the loci metadata
    x2@other$metrics <- x@other$metrics[x@other$metrics["RepAvg"]>=threshold,]  
    cat ("  no. of loci deleted =", (n0-nLoc(x2)), "\nLoci retained =", nLoc(x2))
    
  return(x2)
  
}
