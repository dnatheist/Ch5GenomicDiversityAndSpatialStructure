#### Load DarT files
#necessary libraries

library(plyr)
library(adegenet)
library(tidyr)


#' Convert DarT to genlight
#' 
#' @param dart a dart object created via \code{\link{read.dart}}
#' @param covfilename optional file in csv format with covariates for each individual (see details for explanation)
#' @param probar show progress bar
#' @return a genlight object is returned. Including all available slots are filled. loc.names, ind.names, pop, lat, lon (if provided via the covariate file)
#' @details the covariate file needs to have very specific headings. First an heading called id. Here the ids have to match the ids in the dart object \code{colnames(dart[[4]])}. The following column headings are optional. pop: specifies the population membership of each individual. lat and lon specify spatial coordinates (perferable in decimal degrees WGS1984 format). Additional columns with individual covariates can be imported (e.g. age, gender).
#' @examples
#' \dontrun{
#' dgl <- dart2genlight(dart, "covariates.csv")
#' }
dart2genlight <- function(dart, covfilename=NULL, probar = TRUE)
 {
cat("Start conversion....\n")
cat("Please note conversion of bigger data sets will take some time!\n")
cat("Once finished, we recommend to save the object using save(object, file=\"object.rdata\")\n")

if(probar) pb <- txtProgressBar(min=0, max=1, style=3, initial=NA)

 
 
#### out contains the dart data
nind <- dart[["nind"]]
nsnp <- dart[["nsnp"]]
sraw <- dart[["covmetrics"]]

if (sum(c("SNP", "SnpPosition") %in% names(sraw))!=2)
  {
  cat("Could not find SNP or SnpPosition in Dart file. Check you headers!!!")
  stop()
  }
sdata <- dart[["gendata"]]
#every second line only....
esl = seq(2,nrow(sdata),2)

pos <- sraw$SnpPosition[esl]
alleles <- as.character(sraw$SNP)[esl]
a1 <- substr(alleles,nchar(alleles)-2,nchar(alleles))
a2 <-  sub(">","/", a1)
locname <- paste(sraw$uid[esl], a2,sep="-")
geninddata <- matrix(NA, nrow=nsnp, ncol=nind)
for (i in 1:nind)
{
 isnp = paste(sdata[esl-1,i],sdata[esl,i], sep="/")
 g <- isnp
 g <- gsub("0/1",2,g)
 g <- gsub("1/0",0,g)
 g <- gsub("1/1",1,g)
 g <- gsub("NA/NA",NA,g)
 geninddata[,i] <- as.numeric(g)
if (probar) setTxtProgressBar(pb, i/nind)
}
gout <- new("genlight", gen=t(geninddata), ploidy=2, ind.names=colnames(sdata), loc.names=locname ,loc.all=a2, position=pos, parallel=F)

if (probar) close(pb)
gout@other$metrics <- sraw[esl,]

####
#additional covariates and long lat to the data file are stored in other

if (!is.null(covfilename))
{
cat(paste("Try to add covariate file:", covfilename,".\n"))
###### population and individual file to link AAnumbers to populations...
ind.cov <- read.csv(covfilename, header=T, stringsAsFactors=T)
# is there an entry for every individual

id.col = match( "id", names(ind.cov))

if (is.na(id.col)) {cat ("There is no id column\n") ;stop()} else {

#reorder
if (length(ind.cov[,id.col]) !=length(names(sdata)))  {cat ("Ids of covariate file does not match the number of ids in the genetic file\n") ;stop()} 

if (length( match(names(sdata), ind.cov[,id.col])) ==nind ) cat ("Ids of covariate file are matching!\n") else {cat("Ids are not matching!!!!\n");stop()}
}

ord <- match(names(sdata), ind.cov[,id.col])

 pop.col = match( "pop", names(ind.cov))

 if (is.na(pop.col)) {cat ("Please note:there is no pop column\n") }  else {
    pop(gout) <- as.factor(ind.cov[ord,pop.col])
    cat("Added pop factor.\n")
    }
 
 lat.col = match( "lat", names(ind.cov))
 lon.col = match( "lon", names(ind.cov))

  if (is.na(lat.col)) {cat ("Please note:there is no lat column\n") }
  if (is.na(lon.col)) {cat ("Please note:there is no lon column\n") }
  if (!is.na(lat.col) & !is.na(lon.col))
    {
    gout@other$latlong <- ind.cov[ord,c(lat.col, lon.col)]
    cat("Added latlon data.\n" )
    }

 known.col <- names( ind.cov) %in% c("id","pop", "lat", "lon")
# known.col <- ifelse(is.na(known.col), , known.col)
 other.col <- names(ind.cov)[!known.col]
 if (length(other.col>0) )
    {
    gout@other$covariates<-ind.cov[ord,other.col]
    cat(paste("Added ",other.col," to the other$covariates slot.\n"))
    }
}
gout
}

#' Converts a genlight object to genind object
#' 
#' @param snp a genind object
#' @param probar switch off progress bar
#' @return a genind object, with all slots filled.
#' @details this function uses a faster version of df2genind (from the adgegenet package)

genlight2genind <- function(snp, probar=TRUE)
{

cat("Start conversion....\n")
ptm <- proc.time()[3]
cat("Please note conversion of bigger data sets will take some time!\n" )
cat("Once finished, we recommend to save the object using >save(object, file=\"object.rdata\")\n")
#convert to genind....
x <- as.matrix(snp[,])
if (probar) pb <- txtProgressBar(min=0, max=1, style=3, initial=NA)

for (i in 1:nrow(x))
  {
  for (ii in 1:ncol(x))
    {

    inp <- x[i,ii]
    if (!is.na(inp))
      {
      if (inp==0) x[i,ii] <- "A/A" else if (inp==1) x[i,ii] <- "A/B" else if (inp==2) x[i,ii] <- "B/B"
      }
    }
if (probar)   setTxtProgressBar(pb, i/nrow(x))
  }
  
cat("\nMatrix converted.. Prepare genind object...\n")

gen<-df2genind(x[,], sep="/", ncode=1, ind.names=snp@ind.names, pop = snp@pop, ploidy=2)#, probar=probar)
gen@other <- snp@other

cat(paste("Finished! Took", round(proc.time()[3]-ptm),"seconds.\n") )
gen
}



################################################################################ 
# faster version of df2genind (no longer necessary once adegent2 is out)
################################################################################
# 
# df2genindb<-
# function (X, sep = NULL, ncode = NULL, ind.names = NULL, loc.names = NULL,
#     pop = NULL, missing = NA, ploidy = 2, type = c("codom", "PA"), probar=probar)
# {
# if (probar)   pb <- txtProgressBar(min=0, max=2, style=3, initial=NA)
#     if (is.data.frame(X))
#         X <- as.matrix(X)
#     if (!inherits(X, "matrix"))
#         stop("X is not a matrix")
#     res <- list()
#     type <- match.arg(type)
#     n <- nrow(X)
#     nloc <- ncol(X)
#     ploidy <- as.integer(ploidy)
#     if (ploidy < 1L)
#         stop("ploidy cannot be less than 1")
#     if (is.null(ind.names)) {
#         ind.names <- rownames(X)
#     }
#     if (is.null(loc.names)) {
#         loc.names <- colnames(X)
#     }
#     if (!is.null(pop)) {
#         if (length(pop) != n)
#             stop("length of factor pop differs from nrow(X)")
#         pop <- as.factor(pop)
#     }
#     if (toupper(type) == "PA") {
#         mode(X) <- "numeric"
#         rownames(X) <- ind.names
#         colnames(X) <- loc.names
#         temp <- apply(X, 2, function(c) all(is.na(c)))
#         if (any(temp)) {
#             X <- X[, !temp]
#             warning("entirely non-type marker(s) deleted")
#         }
#         temp <- apply(X, 1, function(r) all(is.na(r)))
#         if (any(temp)) {
#             X <- X[!temp, , drop = FALSE]
#             pop <- pop[!temp]
#             warning("entirely non-type individual(s) deleted")
#         }
#         temp <- apply(X, 2, function(loc) length(unique(loc[!is.na(loc)])) ==
#             1)
#         if (any(temp)) {
#             X <- X[, !temp, drop = FALSE]
#             warning("non-polymorphic marker(s) deleted")
#         }
#         prevcall <- match.call()
#         res <- genind(tab = X, pop = pop, prevcall = prevcall,
#             ploidy = ploidy, type = "PA")
#         return(res)
#     }
#     mode(X) <- "character"
#     if (is.null(sep)) {
#         if (!is.null(ncode)) {
#             temp <- nchar(X[!is.na(X)])
#             if (ncode < max(temp))
#                 stop("some character strings exceed the provided ncode.")
#         }
#         if (is.null(ncode)) {
#             temp <- nchar(X[!is.na(X)])
#             ncode <- max(temp)
#         }
#         if ((ncode%%ploidy) > 0)
#             stop(paste(ploidy, "alleles cannot be coded by a total of",
#                 ncode, "characters", sep = " "))
#     }
#     tempX <- X
#     if (!is.null(sep))
#         tempX <- gsub(sep, "", X)
#     tempX <- gsub("^0*$", NA, tempX)
#     tempX <- gsub("(NA)+", NA, tempX)
#     temp <- apply(tempX, 2, function(c) all(is.na(c)))
#     if (any(temp)) {
#         X <- X[, !temp, drop = FALSE]
#         tempX <- tempX[, !temp, drop = FALSE]
#         loc.names <- loc.names[!temp]
#         nloc <- ncol(X)
#         warning("entirely non-type marker(s) deleted")
#     }
#     temp <- apply(tempX, 1, function(r) all(is.na(r)))
#     if (any(temp)) {
#         X <- X[!temp, , drop = FALSE]
#         tempX <- tempX[!temp, , drop = FALSE]
#         pop <- pop[!temp]
#         ind.names <- ind.names[!temp]
#         n <- nrow(X)
#         warning("entirely non-type individual(s) deleted")
#     }
#     n <- nrow(X)
#     X[is.na(tempX)] <- NA
#     fillWithZero <- function(M, targetN) {
#         naIdx <- is.na(M)
#         keepCheck <- any(nchar(M) < targetN)
#         while (keepCheck) {
#             mat0 <- matrix("", ncol = ncol(M), nrow = nrow(M))
#             mat0[nchar(M) < targetN] <- "0"
#             M <- matrix(paste(mat0, M, sep = ""), nrow = nrow(mat0))
#             keepCheck <- any(nchar(M) < targetN)
#         }
#         M[naIdx] <- NA
#         return(M)
#     }
#     if (is.null(sep) | ploidy == as.integer(1)) {
#         X <- fillWithZero(X, targetN = ncode)
#         splitX <- list()
#         for (i in 1:ploidy) {
#             splitX[[i]] <- substr(X, 1, ncode/ploidy)
#             X <- sub(paste("^.{", ncode/ploidy, "}", sep = ""),
#                 "", X)
#         }
#     }
#     if (!is.null(sep)) {
#         if (ploidy > 1) {
# #            temp <- t(as.matrix(as.data.frame(strsplit(X, sep))))
# 
#             temp <- do.call(rbind, strsplit(X, sep))
#             splitX <- list()
#             for (i in 1:ncol(temp)) {
#                 splitX[[i]] <- matrix(temp[, i], nrow = n)
#             }
#         }
#         else {
#             splitX <- list()
#             splitX[[1]] <- X
#         }
#         temp <- unlist(splitX)
#         temp <- temp[!is.na(temp)]
#         ncode <- max(nchar(temp)) * ploidy
#         splitX <- lapply(splitX, function(Y) fillWithZero(Y,
#             targetN = ncode/ploidy))
#     }
#     loc.all <- list()
#     for (i in 1:nloc) {
#         temp <- unlist(lapply(splitX, function(e) e[, i]))
#         loc.all[[i]] <- sort(unique(temp[!is.na(temp)]))
#     }
#     names(loc.all) <- loc.names
#     temp <- lapply(1:nloc, function(i) matrix(0, nrow = n, ncol = length(loc.all[[i]]),
#         dimnames = list(NULL, loc.all[[i]])))
#     names(temp) <- loc.names
#     findall <- function(cha, loc.all) {
#         if (is.na(cha))
#             return(NULL)
#         return(which(cha == loc.all))
#     }
#     for (k in 1:ploidy) {
#         for (i in 1:n) {
#             for (j in 1:nloc) {
#                 allIdx <- findall(splitX[[k]][i, j], loc.all[[j]])
#                 temp[[j]][i, allIdx] <- temp[[j]][i, allIdx] +
#                   1
#                 if (is.null(allIdx)) {
#                   temp[[j]][i, ] <- NA
#                 }
#             }
# if (probar)        setTxtProgressBar(pb, (i+(k-1)*n)/(n))
#         }
#     }
#     cat("\nAlmost done...\n")
#     flush.console()
#     nall <- unlist(lapply(temp, ncol))
#     loc.rep <- rep(names(nall), nall)
#     col.lab <- paste(loc.rep, unlist(loc.all, use.names = FALSE),
#         sep = ".")
#     mat <- matrix(unlist(temp), nrow = nrow(temp[[1]]))
#     mat <- mat/ploidy
#     colnames(mat) <- col.lab
#     rownames(mat) <- ind.names
#     if (!is.na(missing)) {
#         if (missing == 0) {
#             mat[is.na(mat)] <- 0
#         }
#         if (toupper(missing) == "MEAN") {
#             moy <- apply(mat, 2, function(c) mean(c, na.rm = TRUE))
#             for (j in 1:ncol(mat)) {
#                 mat[, j][is.na(mat[, j])] <- moy[j]
#             }
#         }
#     }
#     prevcall <- match.call()
#     res <- genind(tab = mat, pop = pop, prevcall = prevcall,
#         ploidy = ploidy, type = type)
#     return(res)
# }


##########
### Hs returns NA, if only NA at one loci (which is not good for SNPs)
######### no longer needed?!
# 
#  Hsb <-
# function (x, truenames = TRUE) 
# {
#   if (is.genind(x)) {
#     x <- genind2genpop(x, quiet = TRUE)
#   }
#   if (!is.genpop(x)) 
#     stop("x is not a valid genpop object")
#   if (x@type == "PA") 
#     stop("not implemented for presence/absence markers")
#   x.byloc <- seploc(x, truenames = truenames)
#   lX <- lapply(x.byloc, function(e) makefreq(e, quiet = TRUE, 
#                                              truenames = truenames)$tab)
#   lres <- lapply(lX, function(X) 1 - apply(X^2, 1, sum))
#   res <- apply(as.matrix(data.frame(lres)), 1, mean,na.rm=T)
#   return(res)
# }

############################################ 
### pairwise LD function across all loci
############################################

LDallp <- function(geni, name=NULL, save=TRUE,  nchunks=2, ncores=1, chunkname=NULL)
{
  library(doParallel)
  library(adegenet)
  library(data.table)
  cat(paste("Start to calculate LD for all pairs of loci...\n"))
  cat(paste("Using", ncores,"cores in", nchunks," chunks.\n"))
  cat("Depending on the number of loci this may take a while...\n")
  cat("nchunks specifies the number of steps in the progress bar and the number of intermediate saves, but slows the computation a bit. nchunks = 1 is fastest.\n")
  cat(paste("Seperate all",length(indNames(geni)),"loci...\n"))
  flush.console()
  
  #convert into list of 
  slg <- seploc(geni)
    for (i in 1:length(slg)) slg[[i]] <- slg[[i]]@tab
  
  
  cat(paste("Generate all possible pairs:", length(slg)*(length(slg)-1)/2,"...\n"))
  flush.console()
  allp <- combn(length(slg),2)
  resnames <- c("loc1" ,"loc2","D", "Dprime", "r", "R2", "n", "X2", "p")
  lddone <- NULL
  chunknr <- 0  #to make sure old chunks are not overridden
  if (!is.null(chunkname)) 
  {
    cat("You specified results from a previous runs ...\n")
    cat(paste("Loooking for LD_chunks_", chunkname, "files.\n"))
    chunkfiles <- list.files(pattern=paste0("LD_chunks_",chunkname))
    if (length(chunkfiles>0))
      {
      cat(paste("Found", length(chunkfiles), "file(s).\n"))
      for (i in 1:length(chunkfiles))
        {
        load(paste0("LD_chunks_", chunkname,"_",i,".rdata"))
        lddone[[i]] <- ldc
        }
      chunknr<-length(chunkfiles)
      lddone <- rbindlist(lddone)
      
      cat(paste("Found", nrow(lddone),"pairs...\n"))
      setnames(lddone,resnames)
      done <- nrow(lddone)
      if (done==ncol(allp)) {cat("Already everyting is calculated. If you want to recalculate please delete al LD_chunk files or specify a different chunkname.\n Aborting function...\n");return(lddone)}
      allp <- allp[,-c(1:done)]  
      cat(paste("...only", ncol(allp), "pairs left to be done...\n"))
      } else cat(paste("No chunkfiles with LD_chunks_",chunkname,"_x.rdata found. \nTherefore I restart to calculate all pairs.\n"))
  }
  cat(paste("Calculate LD for all pairs...\n"))
  flush.console()
  n<- ncol(allp)
  ptm <- proc.time()[3]
  #ld functions here
  
  LD.fast <- function(g1,g2,...)
  {
    prop.A <- colMeans(g1, na.rm=T)/2
    if(length(prop.A)!=2) return(list(NA,NA,NA,NA,NA,NA,NA,NA))
    names(prop.A) <- c("A","B")
    prop.B <- colMeans(g2, na.rm=T)/2
    if(length(prop.B)!=2) return(list(NA,NA,NA,NA,NA,NA,NA,NA))
    names(prop.B) <- c("A","B")
    
    major.A <- names(prop.A)[which.max(prop.A)]
    major.B <- names(prop.B)[which.max(prop.B)]
    pA <- max(prop.A, na.rm=TRUE)
    pB <- max(prop.B, na.rm=TRUE)
    pa <- 1-pA
    pb <- 1-pB
    
    Dmin <- max(-pA*pB, -pa*pb)
    pmin <- pA*pB + Dmin;
    
    Dmax <- min(pA*pb, pB*pa);
    pmax <- pA*pB + Dmax;
    
    mja <- which.max(colSums(g1, na.rm=T))
    mjb <- which.max(colSums(g2, na.rm=T))
    counts <- table(g1[,mja], g2[,mjb], useNA="no")
    
    #counts <- counts[c("1","2"),c("1","2")]
    
    #     counts <- table(
    #       allele.count(g1, major.A),
    #       allele.count(g2, major.B) )
    #rownames(counts)<-1:2
    #colnames(counts)<-1:2
    
    n3x3 <- matrix(0, nrow=3, ncol=3)
    colnames(n3x3) <- rownames(n3x3) <- 0:2
    
    # ensure the matrix is 3x3, with highest frequency values in upper left
    for(i in rownames(counts))
      for(j in colnames(counts))
        n3x3[3-as.numeric(i),3-as.numeric(j)] <- counts[i,j]
    
    
    loglik <- function(pAB,...)
    {
      (2*n3x3[1,1]+n3x3[1,2]+n3x3[2,1])*log(pAB) +
        (2*n3x3[1,3]+n3x3[1,2]+n3x3[2,3])*log(pA-pAB) +
        (2*n3x3[3,1]+n3x3[2,1]+n3x3[3,2])*log(pB-pAB) +
        (2*n3x3[3,3]+n3x3[3,2]+n3x3[2,3])*log(1-pA-pB+pAB) + 
        n3x3[2,2]*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB))
    }
    
    # SAS code uses:
    #
    #s <- seq(pmin+0.0001,pmax-0.0001,by=0.0001)
    #lldmx <- loglik(s)
    #maxi <- which.max(lldmx)
    #pAB <- s[maxi]
    
    # but this should be faster:
    solution <- optimize(
      loglik,
      lower=pmin+.Machine$double.eps,
      upper=pmax-.Machine$double.eps,
      maximum=TRUE
    )
    pAB <- solution$maximum
    
    estD <- pAB - pA*pB
    if (estD>0)  estDp <- estD / Dmax else    estDp <- estD / Dmin
    
    n <-  sum(n3x3)
    
    corr <- estD / sqrt( pA * pB * pa * pb )
    
    dchi <- (2*n*estD^2)/(pA * pa * pB* pb)
    dpval <- 1 - pchisq(dchi,1)
    
    retval <- list(
      call=match.call(),
      "D"=estD,
      "D'"=estDp,
      "r" = corr,
      "R^2" = corr^2,
      "n"=n,
      "X^2"=dchi,
      "P-value"=dpval
    )
    class(retval) <- "LD"
    retval
  }
  
  #split into nchunks steps for the progress bar...
  runs <- 1:n
  if (nchunks>length(runs)) nchunks <- length(runs)
  chunks <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  if (nchunks<2 ) splitruns <- list( runs) else splitruns <- chunks(runs, nchunks)
  ldchunks <- list()
  cl<-makeCluster(ncores) #adjust the number of cores of your computer!!!!
  registerDoParallel(cl)
  pb <- txtProgressBar(min=0, max=nchunks, style=3, initial=NA)
  for (i in 1:nchunks)
  {
    iter <- splitruns[[i]]
    if (ncores <2) it2 <- list(iter) else it2 <- chunks(iter, ncores)
    ll <- foreach (ip=1:length(it2), .combine=rbind) %dopar% {
        res <- matrix(NA, ncol=9, nrow=length(it2[[ip]]))
        for (ii in 1:length(it2[[ip]]))
          {
          l1 <- allp[1,it2[[ip]][ii]]
          l2 <- allp[2,it2[[ip]][ii]]
          s1 <- slg[[l1]]
          s2 <- slg[[l2]]
          if (  (length(colSums(s1, na.rm=T)) +  length(colSums(s2, na.rm=T)) ) ==4)
          {
          r <- LD.fast(s1,s2)
          res[ii,] <-c(l1,l2, do.call(cbind, r[2:8]))
          } else res[ii,] <- c(l1,l2, rep(NA,7))
          }
        res
        }
  ldchunks[[i]] <-as.data.frame(ll)
  setTxtProgressBar(pb, i)
  ldc <- ldchunks[[i]]
  save(ldc, file=paste0("LD_chunks_",chunkname,"_",i+chunknr,".rdata"))
  }
stopCluster(cl)
LDres2 <- rbindlist(ldchunks)
#LDres2 <- t(do.call(cbind, ldchunks))
setnames(LDres2,resnames)
if (!is.null(lddone)) LDres2 <- rbindlist(list(lddone, LDres2))
#colnames(LDres2)<- resnames

cat(paste("\n# Simulations:", n,". Took", round(proc.time()[3]-ptm),"seconds.\n"))
if (save) 
{
  if (!is.null(name)) 
  { 
    nobj <- name
    filename <- paste0("LD_",name,".rdata")
  } else
  {
    fx <- deparse(substitute(geni))
    #cat(paste('name:',fx,"\n"))
    fx <- gsub("\\[|\\]|[.]|[:]|[,]|[&]|[%]","_",fx) 
    fx <- gsub("__","_", fx)
    fx <- gsub(" ","", fx)
    fx <- gsub("_$","",fx)
    if (length(fx)>0) 
    {
      nobj <- paste0("LD_",fx)
      filename=paste0(fx,".rdata") 
    } else
    {
      filename="LDallp.rdata"
      nobj <- "LDallp"
    }
  }
  cat(paste0("\n Results are saved as object ", nobj," under ", filename,".\n"))
  (cat(paste("Once you have checked you can delete your LD_chunks_",chunkname,"files.\n")))
  assign(nobj, LDres2)
  save(list=nobj, file=filename)
}
LDres2
}


#####
#filter function to make filtering a bit easier
#####
filter.dart<- function(...)
{
  library(knitr, quietly = T)
  filnames <-   deparse(substitute(list(...)))
  filnames <-  gsub("[()]","",filnames)
  filnames <-  gsub(" ","",filnames)
  filnames <- gsub("list","", filnames)
  filters <- eval(eval(quote(substitute(list(...)))))
  nfilters <- length(filters)
  #filters <- gsub("list( | ")
  names(filters) <- unlist(strsplit( filnames,","))
  lenfil <- unlist(lapply(filters, length))
  tabfilters <- sapply(filters, table)
  tabfilters <- do.call(rbind, tabfilters)
  tabfilters <- cbind(tabfilters, sum=lenfil, freq=round(tabfilters[,2]/lenfil,3))
  print(kable(tabfilters))
  nsnp <- min(lenfil)
  dffil <- do.call(rbind, filters)
  index.all <- colSums(dffil)==nfilters
  index.comb <- rbind(index.all,dffil)
  index.combout <- index.comb
  index.combout[1,] <- index.combout[1,]*2
  rs <-rowSums(index.comb)
  image(1:ncol(index.comb), 1:nrow(index.comb),t(index.combout), col=c("white", "red", "orange"), axes=F, xlab="loci", ylab="",main="Filters for each loci")
  axis(2, at= 1:nrow(index.comb),label= gsub("index.","",row.names(index.comb)), las=2)
  text(nsnp/2,1:nrow(index.comb), rs)
  axis(1)#, at=c(1,round(nsnp/2), nsnp), label=c(1,round(nsnp/2), nsnp))
  lines(c(1,nsnp),c(rep(nrow(index.comb),2))-0.5)
  box()
  #filters
  data.frame(t(index.comb))
}


#########
### HWEperpop
########
HWEperpop <- function(geni, pvalue=0.05, plot=TRUE)
{
  library(SNPassoc)   #package for LD and Hardy Weinberg...
  library(pegas)
  sep.troch <- seppop(geni)  #creates a list of all populations
  #define a function
  hwepop <- function(pop)
  {
    xxx <- as.loci(pop)[,-1]
    sup <-  setupSNP(data.frame(xxx), 1:ncol(xxx)) 
    sum.sup <- summary(sup)
    index.HWE <- sum.sup$HWE < pvalue 
    index.HWE <- ifelse(is.na(index.HWE), FALSE, index.HWE)      #what to do with NA?????
    index.HWE
  }
  indices.HWE.pop <- do.call(rbind,lapply(sep.troch, hwepop)) #let's do that for every population
  if (plot) barplot( colSums(indices.HWE.pop,na.rm=T))
  indices.HWE.pop
}


#### 
## outflank function to identify loci under selection
###
#### loci under selection
#using outflank
#install outflank; see readme at github
#install qvalue from bioconductor
# source("http://bioconductor.org/biocLite.R")
# biocLite("qvalue") 
# install_github("Whitlock/OutFLANK")
outflank <- function(genl, plot=TRUE)
{
  library(OutFLANK)
  # outflank requires (of course) a different format 
  # missing value is 9!!! tempted to rewrite their model to be able to use genlight directly....
  snpmat <- as.matrix(genl)#(matrix(NA, nrow=nind, ncol=nsnp)
  snpmat <- replace(snpmat, is.na(snpmat), 9) 
  mdfm <- MakeDiploidFSTMat(SNPmat = snpmat, list(colnames(snpmat)), list(as.character(genl@pop)))
  #run outflank
  outf <- OutFLANK(mdfm, LeftTrimFraction = 0.05, RightTrimFraction = 0.05, Hmin = 0.1, NumberOfSamples = length(levels(genl@pop)), qthreshold = 0.05)
  if (plot) OutFLANKResultsPlotter(outf)
  #the filter if we believe in Whitlock and why wouldn't we?
  index.outflank <- !(outf$results$OutlierFlag) ## 6650 inliers and 188 outliers 
  index.outflank
  #sum(!index.outflank) #number of outliers 
}





