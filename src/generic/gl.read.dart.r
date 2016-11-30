# An R function to take the input from DArT and convert to agegenet genlight format (gl)
#
# DaRT provide the data as a matrix of entities (individual turtles) across the top and
# attributes (SNP loci) down the side in a format that is unique to DArT. This program
# reads the data in to adegenet format (genlight) for consistency with 
# other programming activity. The script may require modification as DArT modify their
# data formats from time to time.
#
# INPUT
# datafile -- name of an input file in csv obtained from DArT IN 2 ROW FORMAT [Required]
# topskip -- number of header lines to skip which does NOT include the line containing entity (individual) names [Required]
# nmetadata -- number of leading columns containing locus metadata (RepAvg is usually the last column) [Required]
# nas -- a string containing the symbol used for missing data (absences arising from restriction site polymorphism) [Default '-']
# ind.metafile -- name of the csv file containing the metadata associated with each entity (individual) [Default NULL]
#
# OUTPUT
# An Rdata file with the data in genlight format (on disk)
# a genlight file with the data -- returned
#
# USEAGE
# gl <- gl.read.dart(datafile="SNP_DFwt15-1908_scores_2Row.csv", topskip=6, nmetavar=16, nas="-", ind.metafile="metadata.csv" )
#
# Author: Arthur Georges
# Date of last update: 31-Jan-16
#
# Dependencies

library(plyr)
library(adegenet)
library(tidyr)
library(reshape2)

gl.read.dart <- function(datafile, topskip, nmetavar, nas="-", ind.metafile=NULL)
{

# INPUT THE DATA TO PRELIMINARY STORAGE

  cat("Reading data from file:", datafile,"\n")
  cat("  This may take some time, please wait!\n")
  x <- read.csv(datafile, na.strings=nas, skip = topskip, check.names=FALSE)
  cat("The following locus metadata was identified: ", names(x[1:nmetavar]),"\n")
# Error checks
  if(any(names(x) == "AlleleID")) {
    cat("  includes key variable AllelID\n")
  } else {
      cat("Fatal Error: Dataset does not include key variable AlleleID!\n"); stop()
  }
  if(any(names(x) == "SNP")) {
    cat("  includes key variable SNP\n")
  } else {
    cat("Fatal Error: Dataset does not include key variable SNP!\n"); stop()
  }
  if(any(names(x) == "SnpPosition")) {
    cat("  includes key variable SnpPosition\n")
  } else {
    cat("Fatal Error: Dataset does not include key variable SnpPosition!\n"); stop()
  }
  if(any(names(x) == "RepAvg")) {
    cat("  includes key variable RepAvg\n")
  } else {
    cat("Warning: Dataset does not include variable RepAvg which may limit your filtering options in later analyses!\n")
  }  
  
  #length(levels(pop(x)))
  #if (length(levels(pop(x)))) != length(ind.metadata[,1]) {cat"Population list in the metadatafile does not match the population list in the raw datafile\n"); stop()}
  
# Extract names of the entities (individuals)
  ind.names <- colnames(x)[(nmetavar+1):ncol(x)]
  cat("Data identified for ",ncol(x)-nmetavar, "individuals, ", nrow(x)/2, "loci")
# More error checks
  if (length(ind.names)!= length(unique(ind.names))) {
    cat("Fatal Error: Individual names are not unique!\n"); stop()
  }
# Extract the SNP data  
  snpdata <- x[, (nmetavar+1):ncol(x)]
# Extract the standard metadata for loci
 locus.metadata <- x[, 1:nmetavar]
# More error checks
  if(max(snpdata,na.rm=TRUE)!=1 || min(snpdata,na.rm=TRUE)!=0) {
    cat("Fatal Error: SNP data must be 0 or 1!\n"); stop()
  }
# Calculate number of entities (individuals) and attributes (loci)
  nind <- ncol(snpdata)
  nloci <- nrow(locus.metadata)/2

# CONVERT TO GENLIGHT FORMAT

cat("\nStarting conversion to genlight object ....\n")
cat("Please note conversion of bigger data sets will take some time!\n")
pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
getTxtProgressBar(pb)
# 2 Row format -- Consider every second line only
  seqby2 = seq(2,nrow(snpdata),2)
# Extract the SNP position
  pos <- locus.metadata$SnpPosition[seqby2]
# Extract the SNP transitions (metavariable SNP)
  state <- as.character(locus.metadata$SNP)[seqby2]
  state <- substr(state,nchar(state)-2,nchar(state))
  state <-  sub(">","/", state)
  uid <- as.character(locus.metadata$AlleleID)[seqby2]
  a.list <- strsplit(uid, "F")
  uid <- sapply(a.list, "[", 1)
# Create locus names
  locname <- paste(uid, state, sep="")
# Initialize the data matrix
  x <- matrix(NA, nrow=nloci, ncol=nind)
# For each individual, convert the 2 line data to required one line format (0, homozygous reference; 1, heterozygous; 2 homozygous mutant)
  for (i in 1:nind) {
    isnp = paste(snpdata[seqby2-1,i],snpdata[seqby2,i], sep="/")
    g <- isnp
    g <- gsub("0/1",2,g)
    g <- gsub("1/0",0,g)
    g <- gsub("1/1",1,g)
    g <- gsub("NA/NA",NA,g)
    x[,i] <- as.numeric(g)
    setTxtProgressBar(pb, i/nind)
  }
# Create the genlight object
  gl <- new("genlight", gen=t(x), ploidy=2, ind.names=colnames(snpdata), loc.names=locname ,loc.all=state, position=pos, parallel=F)

  close(pb)

# Add in the standard metadata
  gl@other$metrics <- locus.metadata[seqby2,]

#
# Add in extra metadata -- population assignments
if (!is.null(ind.metafile)) {
  cat("Adding population assignments and other additional individual metadata from file :", ind.metafile,"\n")
  ind.metadata <- read.csv(ind.metafile, header=T, stringsAsFactors=T)
  # Remove leading and trailing spaces that could lead to a spurious mismatch
  ind.metadata$id <- gsub("^\\s+|\\s+$", "", ind.metadata$id)
  # Check that the number of individuals in the metafile is the same as in the dataset
#  if(length(ind.metadata[1])!=nInd(gl)) {
#    cat("Fatal Error: Number of individuals in metadata file does not equal number of individuals in the dataset\n"); stop()
#  }
  # Check for an entry for every individual
  id.col = match( "id", names(ind.metadata))
  if (is.na(id.col)) {
    cat ("Fatal Error: No id column present!\n") ;stop()
    } else {
    if (sum(ind.metadata[,id.col] == names(snpdata))== nind ) {
      cat ("Ids of individual metadata file match!\n") 
    }else {
        cat("Fatal Error: Ids in files ",datafile,"and ",ind.metafile," do not match!\n\n");stop()
    }
  }
  pop.col = match( "pop", names(ind.metadata))
  # Check for population assignment
  if (is.na(pop.col)) {
    cat ("Warning: No pop column present\n") 
  } else {
    pop(gl) <- as.factor(ind.metadata[,pop.col])
    cat("Populations assigned to individuals\n")
  }
  # Check for latitude and longitude data  
  lat.col = match( "lat", names(ind.metadata))
  lon.col = match( "lon", names(ind.metadata))
  if (is.na(lat.col)) {
   cat ("Warning: No lat column present\n")
  }
  if (is.na(lon.col)) {
    cat ("Warning: No lon column present\n") 
  }
  if (!is.na(lat.col) & !is.na(lon.col))  {
    gl@other$latlong <- ind.metadata[,c(lat.col, lon.col)]
    cat("Added latlon data\n" )
  }
  # Check for other metadata  
  known.col <- names( ind.metadata) %in% c("id","pop", "lat", "lon")
  # known.col <- ifelse(is.na(known.col), , known.col)
  other.col <- names(ind.metadata)[!known.col]
  if (length(other.col>0) )
  {
    gl@other$covariates<-ind.metadata[,other.col]
    cat("Added ",other.col," to the other$covariates slot\n")
  }
}
# Report
  cat("Genlight object created")

  return <- gl

}
