gl.report.RepAvg <- function(x) {
  
  # A function to report loci from a genlight object based on the DArT
  # estimate of reproducibility in calling the snp (=RepAvg).
  #
  # Input: A genlight object.
  #
  # Output: A genlight object with loci with a missing value rate greater
  # than the specified threshold deleted.
  #
  # Returns the mean RepAVG for the genlight object.
  #
  # Useage: result <- gl.report.RepAvg(gl)
  #
  # Author: Arthur Georges
  # Dependencies: library(adegenet)
  #
  # BEGIN FUNCTION
  #
  cat("\nNo. of loci =", nLoc(x), "\n\n")
  s <- sum(x@other$metrics[,which(names(x@other$metrics)=="RepAvg")]==1)
  cat("    No. loci with RepAvg = 1.00",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(x@other$metrics[,which(names(x@other$metrics)=="RepAvg")]>=0.995)
  cat("    No. loci with RepAvg >= 0.995",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(x@other$metrics[,which(names(x@other$metrics)=="RepAvg")]>=0.99)
  cat("    No. loci with RepAvg >= 0.99",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(x@other$metrics[,which(names(x@other$metrics)=="RepAvg")]>=0.985)
  cat("    No. loci with RepAvg >= 0.985",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(x@other$metrics[,which(names(x@other$metrics)=="RepAvg")]>=0.98)
  cat("    No. loci with RepAvg >= 0.98",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(x@other$metrics[,which(names(x@other$metrics)=="RepAvg")]>=0.975)
  cat("    No. loci with RepAvg >= 0.975",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")
  s <- sum(x@other$metrics[,which(names(x@other$metrics)=="RepAvg")]>=0.97)
  cat("    No. loci with RepAvg >= 0.97",s,"[",round((s*100/nLoc(x)),digits=1),"%] loci\n")

  r <- round(mean(x@other$metrics[,which(names(x@other$metrics)=="RepAvg")]), digits=2)
  
  return(r)
  
}


