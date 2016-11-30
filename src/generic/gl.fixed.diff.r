gl.fixed.diff <- function(x, threshold=0) {

  # A function to generate a matrix of fixed differences from a genelight object
  
  # A function take SNP data grouped into populations in a genlight object 
  # and generate a matrix of fixed differences between populations taken pairwise.
  #
  # Input: A genlight object containing SNP genotypes.
  #        A threshold value for tollerance in when a difference is regarded as fixed. threhsold = 0.05 means that difference [96,4] will
  #        be regarded as a fixed difference.
  #
  # Returns a matrix of fixed differences.
  #
  # Useage: m <- gl.fixed.diff(gl, threshold=0.01)
  #
  # Author: Arthur Georges
  # Updated 1-Feb-16
  #  x=gl
  #  threshold=0
  #
  # BEGIN FUNCTION
  #
  # EXTRACT THE DATA
    gl.mat <- as.matrix(x)
  # Assign population names to the individuals, discarding AA numbers
    rownames(gl.mat) <- pop(gl)

  # REFORMAT THE DATA BY BREAKINGDOWN AGAINST POPULATION AND LOCUS AND CALCULATE ALLELE FREQUENCIES
  # Disaggregate the columns into one column against Population and Locus
    gl.mat.long <- melt(gl.mat, na.rm=FALSE)
    colnames(gl.mat.long) <- c("popn", "locus", "snp")
    rm(gl.mat)
  # Calculate sums and counts broken down by population and locus
    cat("Tallying allele frequencies, this may take some time\n")
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
    gl.mat.sum$frequency<-sum*100/(2*nobs); 
    gl.mat.sum$frequency[is.nan(gl.mat.sum$frequency)]<-NA
    gl.mat.sum$n<-nobs+nmissing
  # Sort the data on locus and population in preparation for analysis
    gl.mat.sum<-gl.mat.sum[order(locus, popn),]
    detach(gl.mat.sum)
    
  # CONVERT MARGINAL FREQUENCIES TO EXTREMES
    gl.mat.sum$frequency <- replace(gl.mat.sum$frequency, gl.mat.sum$frequency<=threshold, 0)
    gl.mat.sum$frequency <- replace(gl.mat.sum$frequency, (gl.mat.sum$frequency>=100-threshold), 100)
    gl.mat.sum$frequency <- replace(gl.mat.sum$frequency, (gl.mat.sum$frequency>threshold) & (gl.mat.sum$frequency<100-threshold), NA)
    gl.mat.sum <- gl.mat.sum[order(gl.mat.sum$locus, gl.mat.sum$popn),]

  # GENERATE A MATRIX OF PAIRWISE FIXED DIFFERENCES
  # Establish an array to hold the fixed differences
    npops<-nlevels(gl.mat.sum$popn)
    nloci<-nlevels(gl.mat.sum$locus)
    fixed <- array(-1, c(npops, npops))
    # Read in sample size thresholds (threshold table)
    thold <- read.csv("threshold.csv", stringsAsFactors=FALSE, header=FALSE); class(threshold)
    tmp <- thold
    for (i in 1:10) {
      for (j in i:10) {
        tmp[j,i] <- thold[i,j]
      }
    }
    thold <- round(tmp*nloci/42000, digits=1) ########## 42,000 Value from simulations ##############
    rm(tmp)
    # Set up the progress counter
    pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
    getTxtProgressBar(pb)
    #cat("/nThreshold values used for significance (0 = significant)\n", thold, "\n\n")
    # Cycle through the data to sum the fixed differences into a square matrix
    for (i in 1:nloci) {                           # For each locus
      countj<-0
      for (j in ((i-1)*(npops)+1):(i*(npops))){    # For each population against that locus
        # bring the threshold index value within bounds 1:10
        n1 <- max(c(1,gl.mat.sum$nobs[j])) # Zero value captured by NA trap below
        n1 <- min(c(10,n1))
        countj<-countj+1
        countk<-0
        for (k in ((i-1)*(npops)+1):(i*(npops))) { # For each population, compare pairwise
          # bring the threshold index value within bounds
          n2 <- max(c(1,gl.mat.sum$nobs[k]))
          n2 <- min(c(10,n2))
          countk<-countk+1 
          # If the two populations to be compared have sufficient sample size (cf threshold)
          if (thold[n1,n2] <= 1) {
            # Compare and if not missing, increment
            cf <- is.fixed(gl.mat.sum$frequency[j],gl.mat.sum$frequency[k])
            if (!is.na(cf)) {
              if (fixed[countj,countk] == -1) {
                fixed[countj,countk] <- cf 
              } else {
                fixed[countj,countk] <- fixed[countj,countk] + cf 
              }
            }
          } 
        }
      }
      setTxtProgressBar(pb, i/nloci)
    }
  # Convert missing values to NA
    fixed[fixed == -1] <- NA
  # Tidy up adding row and column names
    rownames(fixed)<-levels(gl.mat.sum$popn)
    colnames(fixed)<-levels(gl.mat.sum$popn)
  # Return the matrix
    fixed <- fixed[order(rownames(fixed)), order(colnames(fixed))]

}
