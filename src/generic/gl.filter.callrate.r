gl.filter.CallRate <- function(x, threshold=0.95) {
  
  # Filter loci on call rate
  
  # A function to remove loci from a genlight object based on the rate
  # of missing values (=call rate).
  #
  # Input: A genlight object containing SNP genotypes.
  #        A threshold fraction of acceptable non-missing values. Default threshold=0.95.
  #
  # Returns a genlight object 
  #
  # Useage: result <- gl.filter.callrate(gl, 0.95)
  #
  # Author: Arthur Georges
  # Dependencies: library(adegenet)
  # Updated: 27-Feb-2016
  #
  # Dependencies
  library(adegenet)
  
  # BEGIN FUNCTION

    # Determine starting number of loci
    n0 <- nLoc(x)
    cat("Initial no. of loci =", n0, "\n")
    # Remove loci with NA count < 1-threshold 
    x2 <- x[ ,glNA(x,alleleAsUnit=FALSE)<((1-threshold)*nInd(x))]
    x2@other$metrics <- x@other$metrics[glNA(x,alleleAsUnit=FALSE)<((1-threshold)*nInd(x)),]    
    cat ("  no. of loci deleted =", (n0-nLoc(x2)), "\nLoci retained =", nLoc(x2),"\n")
    
    return(x2)

}

