#' Import DarT data to R
#'
#' @param filename path to file (csv file only currently)
#' @param nas a character specifying NAs (default is "-")
#' @param topskip a number specifying the number of rows to be skipped (something between 4-6 normally, has no default currently)
#' @param stdmetrics a vector of column headings that are extracted. AlleleID and CloneID are compulsory, the rest are needed for filtering
#' @param addmetrics add additional headers/columns by name
#' @param lastmetric specifies the last non genetic column (Default is "RepAvg"). Be sure to check if that is true, otherwise the number of individuals will not match. You can also specify the last column by a number.
#' @return a list of length 4. #individuals, #snps, #non genetic metrics, #genetic data (still two line format, rows=snps, columns=individuals)
#' @examples
#' dart1 <- read.dart("dart.csv", topskip=6)

read.dart <- function(filename, nas = "-", topskip, stdmetrics =c("AlleleID","CloneID", "SNP","SnpPosition","RepAvg","CallRate", "AvgCountRef", "AvgCountSnp", "FreqHomRef", "FreqHomSnp", "FreqHets","OneRatioSnp"), addmetrics=NULL, lastmetric ="RepAvg")
{
  snpraw <- read.csv(filename,   na.strings=nas, skip = topskip, check.names=FALSE)
  
  
  if (is.character(lastmetric))
  {
    lmet <- which(lastmetric==names(snpraw))
    if (length(lmet)==0)  cat ("Could not determine data columns based on RepAvg!\n")
  } else lmet  <- lastmetric
  
  ind.names <- colnames(snpraw)[(lmet+1):ncol(snpraw) ]
  
  if (length(ind.names)!= length(unique(ind.names))) {cat("Individual names are not unique. You need to change them!\n"); stop()}
  
  datas <- snpraw[, (lmet+1):ncol(snpraw)]
  
  stdmetricscols <- which(  names(snpraw)   %in% stdmetrics )
  
  if (length(stdmetricscols) != length(stdmetrics))
  { cat(paste("\nCould not find all standard metrics.\n",stdmetrics[!(stdmetrics %in% names(snpraw)   )]
              ," is missing.\n Carefully check the spelling of your headers!\n"))
    stop()
  }
  
  if (!is.null(addmetrics)) 
  {
    addmetricscols <- which(  names(snpraw)   %in% addmetrics )
    if (length(addmetricscols) != length(addmetrics))
    { cat(paste("\nCould not find all additional metrics.\n",addmetrics[!(addmetrics %in% names(snpraw)   )]
                ," is missing.\n Carefully check the spelling of your headers! or set addmetrics to NULL\n"))
      stop()
    }
    stdmetricscols <- c(stdmetricscols, addmetricscols)
  } 
  cat ("Added the following cometrics:\n")
  cat (paste(paste(names(snpraw)[stdmetricscols], collapse=" "),".\n"))
  covmetrics <-  snpraw[,stdmetricscols]
  
  #####Various checks (first are there two rows per allele?
  
  covmetrics$CloneID = as.character(covmetrics$CloneID)
  
  #check that there are two lines per locus...
  covmetrics = tidyr::separate(covmetrics, CloneID, into  = c("clid","clrest"),sep = "\\|", extra="merge")
  
  covmetrics$AlleleID = as.character(covmetrics$AlleleID)
  
  #check that there are two lines per locus...
  covmetrics = tidyr::separate(covmetrics, AlleleID, into  = c("allid","alrest"),sep = "\\|", extra="merge")
  
  
  covmetrics$uid <- paste(covmetrics$allid, 1:nrow(covmetrics),sep="-")
  ### there should be only twos (and maybe fours)
  tt <- table(table(covmetrics$clid) )
  cat(paste("Number of rows per AlleleID. Should be only multiples of 2s:\n "))
  cat(names(tt))
  cat("\n")
  nind <- ncol(datas)
  nsnp <- nrow(covmetrics)/2
  
  cat(paste("Recognised:", nind, "indivduals and",nsnp," SNPs in ", filename,"\n"))
  
  out <- list(nind=nind, nsnp=nsnp, covmetrics= covmetrics, gendata =datas)
  
  out
}

