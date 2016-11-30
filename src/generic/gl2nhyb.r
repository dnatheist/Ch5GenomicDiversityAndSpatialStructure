gl2nhyb <- function(gl, outfile="nhyb.out", p0=NULL, p1=NULL, threshold=0) {
  
  # A function to remove loci from a genlight object where nominated parental 
  # populations share alleles and to create a datafile suitable to input to newHybrids.
  #
  # This function compares two sets of parental populations to identify loci
  # that exhibit a fixed difference, and retains them in an output gl object.
  # A fixed difference occurs when a SNP allele is present in all individuals 
  # of one population and absent in the other. There is provision for setting
  # a level of tollerance, e.g. threshold = 0.05 which considers alleles present
  # at 5% in one population and 95% in the other to be a fixed difference.
  # 
  # Input: Genlight object containing SNP genotypes.
  #        List of "parental" populations for species 1.
  #        List of "parental" populations for species 2.
  #        Threshold value for identifying fixed differences.
  #
  # Output: A file in NewHybrids format
  #
  # Returns The reduced genlight object
  #
  # Dependencies: library(MASS) library(reshape2) library(plyr)
  #
  # Useage: m <- gl2nhyb(gl, c("Pop1", "Pop4"), c("Pop7", "Pop9"), threshold=0)
  #
  # Author: Arthur Georges
  #
  # BEGIN FUNCTION
  #
  # load support functions
  source("src/generic/is.fixed.r")

  library(reshape2)
  library(plyr)
  library(MASS)

  # EXTRACT THE SNP DATA
    cat("Extracting the SNP data\n")
    gl.mat <- as.matrix(gl)
    ind.names <- rownames(gl.mat)
    fixed.loci <- NA
    length(fixed.loci) <- 0

  # PROCESS AS FOLLOWS IF BOTH PARENTAL POPULATIONS ARE SPECIFIED
    if (!is.null(p0) & !is.null(p1)) {
      cat("  Both parental populations have been specified \n")
      # Replace those names with Z0 for Parental Population 0, and z1 for Parental Population 1
      rownames(gl.mat) <- replace(rownames(gl.mat), is.element(rownames(gl.mat), p0), "z0")
      rownames(gl.mat) <- replace(rownames(gl.mat), is.element(rownames(gl.mat), p1), "z1")
      # Create a vector containing the flags for the parental population
      par.names <- rownames(gl.mat)
      par.names <- replace(par.names, (par.names != "z0" & par.names != "z1"), " ")
      # Discard non-parental populations
      gl.mat <- gl.mat[(rownames(gl.mat)=="z0" | rownames(gl.mat)=="z1"),]
      # REFORMAT THE DATA BY BREAKINGDOWN AGAINST POPULATION AND LOCUS AND CALCULATE ALLELE FREQUENCIES
      # Disaggregate the columns into one column against Population and Locus
      gl.mat.long <- melt(gl.mat, na.rm=FALSE)
      colnames(gl.mat.long) <- c("popn", "locus", "snp")
      # Calculate sums and counts broken down by population and locus
      gl.mat.sum<-ddply(gl.mat.long, 
                        c("popn","locus"), 
                        summarize, 
                        sums=sum(as.numeric(as.character(snp)),na.rm=TRUE), 
                        count=sum(!is.na(snp)), 
                        missing=sum(is.na(snp))
      )
      names(gl.mat.sum)<-c("popn","locus","sum","nobs","nmissing")
      rm(gl.mat.long)
      # Calculate some new variables, and in particular, percentage frequencies
      attach(gl.mat.sum)
      gl.mat.sum$frequency<-sum*100/(2*nobs)
      gl.mat.sum$frequency[is.nan(gl.mat.sum$frequency)]<-NA
      gl.mat.sum$n<-nobs+nmissing
      # Sort the data on locus and population in preparation for analysis
      gl.mat.sum<-gl.mat.sum[order(locus, popn),]
      detach(gl.mat.sum)
      # CONVERT MARGINAL FREQUENCIES TO EXTREMES
      gl.mat.sum$frequency <- replace(gl.mat.sum$frequency, gl.mat.sum$frequency<=threshold, 0)
      gl.mat.sum$frequency <- replace(gl.mat.sum$frequency, (gl.mat.sum$frequency>=100-threshold), 100)
      gl.mat.sum$frequency <- replace(gl.mat.sum$frequency, (gl.mat.sum$frequency>threshold) & (gl.mat.sum$frequency<100-threshold), NA)
      gl2 <- gl.mat.sum[order(gl.mat.sum$locus, gl.mat.sum$popn),]
      rm(gl.mat.sum)
      # IDENTIFY LOCI WITH FIXED DIFFERENCES BETWEEN P0 AND P1
      cat("  Identifying loci with fixed difference between parental stocks\n")
      # Establish a vector to hold the loci
      npops<-2
      nloci<-nlevels(gl2$locus)
      nloci
      fixed.loci <- NA
      length(fixed.loci) <- nloci*2
      # Cycle through the data to identify the fixed loci
      for (i in seq(1, nloci*2, 2)) {                          
        if (as.character(gl2$locus[i]) != as.character(gl2$locus[i+1])) {
          cat("Error: Loci do not agree on is.fixed comparison\n")
        }
        if (!is.na(is.fixed(gl2$frequency[i],gl2$frequency[i+1]))) {
          if (is.fixed(gl2$frequency[i],gl2$frequency[i+1])) {
            fixed.loci[i] <- as.character(gl2$locus[i])
          }
        }
      }
      # Remove the NAs
      fixed.loci <- fixed.loci[!is.na(fixed.loci)]
      # If no loci remain, set the data matrix to be the original matrix
      if (length(fixed.loci) == 0) {
        gl2 <- as.matrix(gl)
        cat("  No fixed differences between parental populations \n   Retaining all loci\n")
        flag <- "bothparnonefixed"
      } else {  
      # Set the data matrix to contain only the fixed differences
        gl.reduced <- gl[, (locNames(gl) %in% fixed.loci)]
        gl2 <- as.matrix(gl.reduced)
        cat("  Only loci with fixed differences between parental populations have been retained\n")
        flag <- "bothpar"
      }
    }
    
    # PROCESS AS FOLLOWS IF ONLY ONE PARENTAL POPULATION IS SPECIFIED
    if ((!is.null(p0) & is.null(p1)) || (is.null(p0) & !is.null(p1))) {
      gl2 <- as.matrix(gl)
      cat("  Only one parental population specified, retaining all loci\n")
      flag <- "onepar"
    }
    
    # PROCESS AS FOLLOWS IF NO PARENTAL POPULATION IS SPECIFIED
    if (is.null(p0) & is.null(p1)) {
      gl2 <- as.matrix(gl)
      cat("  No parental population specified, retaining all loci\n")
      flag <- "nopar"
    }
    
    # Specify the parental names
    # Replace those names with Z0 for Parental Population 0, and z1 for Parental Population 1
    gl.mat <- as.matrix(gl)
    rownames(gl.mat) <- replace(rownames(gl.mat), is.element(rownames(gl.mat), p0), "z0")
    rownames(gl.mat) <- replace(rownames(gl.mat), is.element(rownames(gl.mat), p1), "z1")
    # Create a vector containing the flags for the parental population
    par.names <- rownames(gl.mat)
    par.names <- replace(par.names, (par.names != "z0" & par.names != "z1"), " ")
    # CREATE THE NEWHYBRIDS INPUT FILE
    # Convert to NewHrbrids lumped format
    # Recode values
    #  0 to 11
    #  1 to 12
    #  2 to 22
    #  NA to 0
    cat("Converting data to NewHybrids format\n")
    gl2[gl2 == 2] <- 22
    gl2[gl2 == 1] <- 12
    gl2[gl2 == 0] <- 11
    gl2[is.na(gl2)] <- 0
    n.loci <- ncol(gl2)
    # Create sequential row number
    rownum <- seq(1:nrow(gl2))
    # Bind to the matrix
    gl2 <- data.frame(gl2)
    
    #NOTE: NewHybrids does not seem to work with the addition of specimen names.
    #if (flag=="bothpar" || flag == "onepar") {
    #  gl2 <- cbind(rownum, " n ", ind.names, par.names, gl2)
    #  metarows <- 4
    #  cat("  Adding sequential number, individual names, and flagging parental stock\n")
    #}
    #if (flag=="nopar") {
    #  gl2 <- cbind(rownum, " n ", ind.names, gl2)
    #  metarows <- 3
    #  cat("  Adding sequential number and individual names, parental stock not identified\n")
    #}
    
    if (flag=="bothpar" || flag == "onepar") {
      gl2 <- cbind(rownum, par.names, gl2)
      metarows <- 2
      cat("  Adding sequential number and flagging parental stock\n")
    }
    if (flag=="nopar") {
      gl2 <- cbind(rownum, gl2)
      metarows <- 1
      cat("  Adding sequential number, parental stock not identified\n")
    }
    
    # Clear row and column names
    rownames(gl2) <- NULL
    colnames(gl2) <- NULL
    # Output only the first 200 loci, to accommodate limits in newhybrids
    if(n.loci > 200) {
      cat("   Warning: NewHybrids may have difficulty if the number of loci exceed 200\n")
      cat("     Retaining first 200 loci only\n")
      gl2 <- gl2[,1:(metarows+200)]
      n.loci <- 200
    }
    
    # Output data file
    cat("Writing the NewHybrids input file", outfile, "\n")
    sink(outfile)
    cat(c("NumIndivs ", nrow(gl2), "\n"))
    cat(c("NumLoci ", n.loci, " \n"))
    cat(c("Digits 1\nFormat Lumped \n"))
    write.matrix(gl2[,1:(ncol(gl2))], sep=" ")
    sink()
    
    return(gl)
    
}  
