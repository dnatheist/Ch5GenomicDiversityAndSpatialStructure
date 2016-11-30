gl.report.CallRate <- function(x) {
  
  # Report loci on call rate
  
  # A function to report call rate of loci in a genlight object
  # of missing values (=call rate).
  #
  # Input: A genlight object containing SNP genotypes.
  #
  # Output: Breadown of loci based on call rate (to screen).
  #
  # Returns the mean callrate for the genlight object.
  #
  # Useage: result <- gl.report.CallRate(gl)
  #
  # Author: Arthur Georges
  # Dependencies: library(adegenet)
  #
  # BEGIN FUNCTION
  #
  cat("\nNo. of loci =", nLoc(x), "\n\n")
  s <- sum(glNA(x,alleleAsUnit=FALSE)==0)
  cat("    Zero missing values =",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(glNA(x,alleleAsUnit=FALSE)<=((1-0.99)*nInd(x)))
  cat("  < 1% missing values =",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(glNA(x,alleleAsUnit=FALSE)<=((1-0.975)*nInd(x)))
  cat("  < 2.5% missing values =",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(glNA(x,alleleAsUnit=FALSE)<=((1-0.95)*nInd(x)))
  cat("  < 5% missing values =",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(glNA(x,alleleAsUnit=FALSE)<=((1-0.90)*nInd(x)))
  cat("  < 10% missing values =",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(glNA(x,alleleAsUnit=FALSE)<=((1-0.80)*nInd(x)))
  cat("  < 20% missing values =",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(glNA(x,alleleAsUnit=FALSE)<=((1-0.75)*nInd(x)))
  cat("  < 25% missing values =",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(glNA(x,alleleAsUnit=FALSE)<=((1-0.5)*nInd(x)))
  cat("  < 50% missing values =",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n\n")

  r <- round(mean((100-(glNA(x,alleleAsUnit=TRUE)*100/nInd(x)))), digits=2)
  
  return(r)
  
}


