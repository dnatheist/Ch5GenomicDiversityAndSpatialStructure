# ================================================== #
# BCB420 / JTB2020                                   #
# March 2014                                         #
# Clustering                                         #
#                                                    #
# Boris Steipe <boris.steipe@utoronto.ca>            #
# ================================================== #
# http://steipe.biochemistry.utoronto.ca/abc/index.php/R_Gene_expression_clustering
# This is an R script for the exploration of clustering
# methods, especially on gene expression data. You will
# get the most out of it if you execute the commands one
# by one, and try to understand completely what they
# mean. All functions have help pages that can be accessed
# through the R console, and if you can't figure
# something out, or if you would like to change a particular
# command or function and don't know how, you are 
# encouraged to contact me for guidance and so I can
# update and improve this script. Be curious, play with 
# this and enjoy!


# ==================================================
# Data
# ==================================================
# 
# Let's find some cell-cycle data in GEO, for clustering.
# The goal is to identify coregulated genes, but we don't
# know what their response to time in the cell-cycle will
# be. Going up? Going down?

# The following segement of code is slightly adapted from
# performing a standard GEO2R analyis on the NCBI website for 
# "Cell cycle expression profiles in HeLa cells" (GSE26922)
# see: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26922
# The dataset contains triplicate measurements for t0 (blocked) and 
# t= 2,4,6,8 and 12h post block-release.

# First, we need to install some analysis packages from bioconductor
# The following commands do for the bioconductor system what 
# install.packages() does for CRAN.

source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")
biocLite("GEOquery")
biocLite("limma")

# Then we load the libraries....
library(Biobase)
library(GEOquery)
library(limma)

# Then load series and platform data from GEO ...
gset <- getGEO("GSE26922", GSEMatrix =TRUE)

if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Check what we have

head(gset)

# The code below is pretty much verbatim GEO2R ...

# Make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Group names for all samples
sml <- c("G0","G0","G0","G1","G1","G1",
         "G2","G2","G2","G3","G3","G3",
         "G4","G4","G4","G5","G5","G5");

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# Set up the data and proceed with analysis
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G5-G0, G1-G0, G2-G1, G3-G2, G4-G3, G5-G4, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

# load NCBI platform annotation
gpl <- annotation(gset)
platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))

# replace original platform annotation
tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gset), "ID"))]
tT <- merge(tT, ncbifd, by="ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))

#### so far, the GEO2R code ....
# It has returned to us the 250 top-differentially expressed
# genes across the groups we have defined.

# ==================================================
# Exploring the data
# ==================================================
# 
# Let's familiarize ourselves a bit with the structure 
# of the data.

head(tT)

# The top gene has the ID 8117594: what are the original values?

exprs(gset)["8117594",]

# Note how we use a string constant to get a data row from he table
# of expression values. We can also use a vector - the expression below 
# returns the data rows for the top three differentially expressed genes:
exprs(gset)[c("8117594","7900167", "8151871"),]

# ==================================================
# Processing the data for cluster analysis
# ==================================================
# 
# For cluster analysis, it's useful to make a table from
# these data that contains only numbers, and just a single
# value (mean) for the biological replicates.

gSym <- tT$Gene.symbol

dat <- c()
for (i in 1:nrow(tT)) {
        v <- c()
        v  <- c(v,  mean(exprs(gset)[tT$ID[i], 1:3]))  # t = 0
        v  <- c(v,  mean(exprs(gset)[tT$ID[i], 4:6]))  # t = 2 hours
        v  <- c(v,  mean(exprs(gset)[tT$ID[i], 7:9]))  # etc...
        v  <- c(v,  mean(exprs(gset)[tT$ID[i], 10:12]))
        v  <- c(v,  mean(exprs(gset)[tT$ID[i], 13:15]))
        v  <- c(v,  mean(exprs(gset)[tT$ID[i], 16:18]))
        dat <- rbind(dat, v)
}
colnames(dat) <- c("t0", "t2", "t4", "t6", "t8", "t12")

# We could use the IDs as rownames, like so ...
#    rownames(dat) <- tT$ID
# ... or the gene symbols, since the IDs don't really
# tell us anything useful. But: are the gene symbols unique?
# If they are not unique, we'll have all sorts of trouble later
# on when we select by rowname...
# R has the function duplicated() to find repeated values
# in a vector...
as.character(gSym[duplicated(gSym)])

# Ha! There are eleven symbols that reappear. Some of them  are
# formatted like "FAM72A///FAM72D///FAM72B" which may mean that
# a spot on the microarray doesn't distinguish between three
# isoforms ... and some are simply the empty string "".
# Since duplicated() gives us a convenient logical vector to 
# identify them, we can simply remove them. This is good enough
# for our clustering exercise, for "real" work we should go back
# to the platform information, find out why there are duplicated
# gene symbols, and address this issue.
dat <- dat[!duplicated(gSym), ]
rownames(dat) <- gSym[!duplicated(gSym)]

# This completes the creation of our expression dataset for clustering.

# You could store the data in a local file ...
write.csv(dat, file="/my/R/folder/GSE26922.dat")

# and then read it back in like so...
dat <- read.csv(file="/my/R/folder/GSE26922.dat",
                row.names = 1,
                header=TRUE)


# ==================================================
# First explorations: Heatmap
# ==================================================
# 
# Heatmaps are a staple of gene expression analysis.
# You can tweak many of the parameters, but for a first look
# will just heatmap the data with default parameters.

# This is a standard view that can be applied to all manners
# of multidimensional data, not just genes.
heatmap(dat)

# Just for illustration and readability let's map only
# every fifth gene
heatmap(dat[seq(1, nrow(dat), by=5), ])

# Study the heatmap, and consider what it tells you.
# For example, there seem to be genes that are low at t4 and t6
# but high at t0, t2 ...
set1 <- c("TPX2", "CCNA2", "AURKA", "CEP55", "CCNB1")
# ... and there are genes for which the inverse is true:
set2 <- c("MAB21L3", "CCNE1", "TCF19///TCF19", "ZBTB14") 

# We can use a "parallel coordinates" plot - matplot()
# to look at the actual expression levels. Note that 
# matplot expects the values column-wise ordered, thus
# we have to transpose - t() - the data!
matplot(t(dat[set1,]), 
        type="l", lwd=2, col="skyblue", lty=1, 
        ylim=c(8,14), xlab="time", ylab="log expression value")

# Then we can use lines() to superimpose the genes for set2.
# No transpose here :-)
for (i in 1:length(set2)) {
        lines(dat[set2[i], ], type="l", lwd=2, col="firebrick")
}

# Indeed, these genes - visibly different in the heatmap
# are mutualy similar in their expression profiles and different
# from each other.

# ==================================================
# Hierarchical clustering
# ==================================================
# 
# Hierarchical clustering is probably the most basic technique.
# The dendrograms on the rows and columns of the heatmap
# were created by hierarchical clustering.

# For hierarchical clustering, first we need to produce
# a distance table. There are many ways to define distances
# let's just go with the default: "Euclidian distance".
distDat <-dist(dat)

# Then we use the clustering distance matrix to produce a 
# dendrogram in which the most similar genes are connected, and then
# similar genes or connected groups are added. There are 
# several ways to define "most-similar", lets just go with the
# default for now: "complete linkage" hierarchical clustering
hc <- hclust(distDat)

plot(hc)

# Not bad. But do note that both distance as well as clustering
# method matter, and there is not really a "best" way that
# works for all data. You'll need to explore: what you are looking for
# is a distance metric that gives the clearest block structure.

df <- function(x) dist(x, method="euclidian")
heatmap(dat[seq(1, nrow(dat), by=4), ], labRow="", labCol="", distfun = df)

df <- function(x) dist(x, method="canberra")
heatmap(dat[seq(1, nrow(dat), by=4), ], labRow="", labCol="", distfun = df)

df <- function(x) dist(x, method="maximum")
heatmap(dat[seq(1, nrow(dat), by=4), ], labRow="", labCol="", distfun = df)

df <- function(x) dist(x, method="minkowski")
heatmap(dat[seq(1, nrow(dat), by=4), ], labRow="", labCol="", distfun = df)

# You are not confined to the default distance functions, it
# is quite straightforward to define your own, for example
# using correlation properties. Here is a distance function
# defined as 1- abs(pearson correlation)...

df <- function(x) as.dist(1 - abs(cor(t(x))))
heatmap(dat[seq(1, nrow(dat), by=4), ], labRow="", labCol="", distfun = df)

# Back to our original dendrogram ...
plot(hc)

# To get clusters from a dendrogram, we need to "cut" it at some
# level. The tree then falls apart into sub-trees and each of these
# is one "cluster"...

# Draw rectangles at different cut-levels, to give the desired number
# of clusters.
rect.hclust(hc,k=2)
rect.hclust(hc,k=5)
rect.hclust(hc,k=10)
rect.hclust(hc,k=20)
rect.hclust(hc,k=50)

# Now retrieve the actual indices and use them to generate
# parallel coordinate plots.

class <-cutree(hc, k = 20)

# Explain the output...
class

# The table() function allows us to count the number of
# occurences in each class ...
table(class)
sort(table(class))

# Let's plot the four largest classes (in parallel, into the same window)
# Look at this carefully. See how the selection statement on class
# generates a logical vector: TRUE in all rows for which the statement is true,
# and how this is used to select the rows of dat that we want to plot ...

oPar <- par(mfrow=c(2,2))
matplot(t(dat[class==2,]),type="l", xlab="time",ylab="log expression value")
matplot(t(dat[class==4,]),type="l", xlab="time",ylab="log expression value")
matplot(t(dat[class==11,]),type="l", xlab="time",ylab="log expression value")
matplot(t(dat[class==6,]),type="l", xlab="time",ylab="log expression value")
par(oPar)


# As an alternative, try Wards- linkage clustering (and read up on the 
# options: single-, complete- and average-linkage clustering)
hc.ward <-hclust(distDat, method = "ward", members=NULL)

plot(hc.ward)

# draw rectangles
rect.hclust(hc.ward,k=9)

# This looks reasonable ...
# Now retrieve the actual indices and use them to generate
# paralell coordinate plots.

class.ward<-cutree(hc.ward, k = 9)
sort(table(class.ward))

# get some nice colors
install.packages("RColorBrewer")
library(RColorBrewer)
# what spectra are there in the package .. ?
display.brewer.all()

niceCols <- brewer.pal(9, "Spectral")

oPar <- par(mfrow=c(3,3))
for (i in 1:9) {
        matplot(t(dat[class.ward == i,]),
                type="l", col=niceCols[i], 
                xlab="time",ylab="log expression value")
}
par(oPar)

# While this may be aesthetically somewhat satisfactory, it is clear that
# the clusters are not homogenous as we might need them for biological
# interpretation. This is a general problem with clustering methods that
# fix the number of cluster centres either directly as in Kmeans (see
# below), or indirectly by cutting trees at a fixed level. It is also
# a problem with the data, where differences in absolute values might 
# override separation into clusters that might better be defined in terms
# of relative values.

# Here is a package that adresses the dynamic range problem.
# Read about it here: http://cran.r-project.org/web/packages/dynamicTreeCut/dynamicTreeCut.pdf
install.packages("dynamicTreeCut")
library(dynamicTreeCut)

class.dynamic <- cutreeDynamic(dendro = hc.ward, distM = as.matrix(distDat), cutHeight=100)

niceCols <- brewer.pal(8, "Spectral")

oPar <- par(mfrow=c(3,3))
for (i in 1:8) {
        matplot(t(dat[class.dynamic == i,]),
                type="l", 
                col=niceCols[i],
                xlab="time",
                ylab="log expression value")
}
par(oPar)

# One thing our clustering algorithms do is to pull apart profiles
# that have similar shape, but different absolute levels. This is
# because we have not normalized our data. Let's thus try
# clustering merely based on profile shape, i.e.
# relative expression levels, by scaling all rows between zero
# and one.

datNorm <- t(apply(dat, 1, function(x)(x-min(x))/(max(x)-min(x))))
distDatNorm <-dist(datNorm)

hc.Norm <-hclust(distDatNorm)

class.dynamic <- cutreeDynamic(dendro = hc.Norm, distM = as.matrix(distDatNorm), cutHeight=15)

niceCols <- brewer.pal(6, "Spectral")

oPar <- par(mfrow=c(3,2))
for (i in 1:6) {
        matplot(t(datNorm[class.dynamic == i,]),
                type="l", 
                col=niceCols[i],
                xlab="time",
                ylab="log expression value")
}
par(oPar)

# With hierarchical clustering, this is probably as good
# as we can get - the clusters are of reasonable size -
# but from a biological point of view one would argue 
# that several of them are not really different.


# ==================================================
# Partitioning clustering
# ==================================================

# === K-means ======================================

# K-means clusters by assigning elements to a fixed
# number of cluster centres, so that similarity
# within a cluster is maximized.

?kmeans

k <- 4
cl<-kmeans(dat, k)

niceCols <- brewer.pal(k, "Spectral")

plot(dat[,"t0"],dat[,"t6"],col=niceCols[cl$cluster])
points(cl$centers, col = niceCols[1:k], pch = 8, cex=2)

# But: be aware ...
# ... K-means does not guarantee a globally optimal solution,
# merely a locally converged one.

# === K-medoids ======================================

# load library "cluster" for K-medoid partitioning
library(cluster)
set.seed(112358)
k <- 4
cl<-pam(dat, 4)
plot(dat[,"t0"],dat[,"t6"], col=niceCols[cl$cluster])
plot(cl) # shows boundary and silhouette plots

# ==================================================
# Affinity propagation clustering
# ==================================================
# Based on B. J. Frey and D. Dueck. Clustering by
# passing messages between data points. 
# Science, 315(5814):972–976, 2007

install.packages("apcluster")
library(apcluster)

apRes <- apcluster(negDistMat(r=2), dat)
apRes

heatmap(apRes)

# try this on the normalized data
apRes <- apcluster(negDistMat(r=2), datNorm)
heatmap(apRes)

# The clear and pronounced block structure shows that this
# is a successful clustering...

length(apRes)
cutree(apRes)

oPar <- par(mfrow=c(3,2))
matplot(t(datNorm[unlist(apRes[2]),]),type="l",xlab="time",ylab="log expression value")
matplot(t(datNorm[unlist(apRes[3]),]),type="l",xlab="time",ylab="log expression value")
matplot(t(datNorm[unlist(apRes[8]),]),type="l",xlab="time",ylab="log expression value")
matplot(t(datNorm[unlist(apRes[14]),]),type="l",xlab="time",ylab="log expression value")
matplot(t(datNorm[unlist(apRes[15]),]),type="l",xlab="time",ylab="log expression value")
matplot(t(datNorm[unlist(apRes[9]),]),type="l",xlab="time",ylab="log expression value")
par(oPar)



# ==================================================
# t-stochastic Neighbour embedding
# ==================================================

# tsne - pioneered by Geoff Hinton - is an embedding
# procedure that guarantees that points in a high-
# dimensional feature space that are close together
# remain close together when projected into a low-
# dimensional space. It can give astoundingly good
# results that help figuring out the internal structure
# of a dataset.

install.packages("tsne")
library(tsne)

# The clustering method uses a "callback" to  execute
# a plotting routine as it improves its embedding
# through many cycles.

# define a plotting routine
plotProgress <- function(x){ 
        plot(x, type='n');
        #	text(x, labels = rownames(dat), cex=0.5)
        points(x, pch=21, col="#6677FF", bg="firebrick")
}

set.seed(112358)
tsneDat <- tsne(dat, epoch_callback = plotProgress)

# presumably the outliers here are due to the non-normalized
# data ....

set.seed(112358)
tsneDatNorm <- tsne(datNorm, epoch_callback = plotProgress)


# I've run this many times to find a
# seed that gives a low error, and run until the iteration
# stabilizes, it's worthwhile experimenting a bit ...

# ... this will run a few minutes  :-)
set.seed(270745)
tsneRef <- tsne(datNorm,
                epoch = 1000, max_iter=8000,
                epoch_callback = plotProgress, perplexity=10)


# We see that the normalized data is better behaved. And we see
# intriguing structure in the data. What does this mean though?
# How do the results compare to what we have seen previously?
# Lets plot this plot with the color values we got from our
# affinity-propagation (AP) clustering.

# ==================================================
# Digression: color scales
# ==================================================

# Our AP clustering had identified 16 clusters. The color brewer palette
# "Spectrum" supports only a max. of 11 colors - and there is 
# some wisdom in not overloading your plot with colors. But lets define
# a color palette for more colors anyway. You could use one of the 
# inbuilt palettes ...

n<-20
pie(rep(1,n), col=rainbow(n), main="rainbow()")
pie(rep(1,n), col=heat.colors(n), main="heat.colors()")
pie(rep(1,n), col=terrain.colors(n), main="terrain.colors()")
pie(rep(1,n), col=topo.colors(n), main="topo.colors()")
pie(rep(1,n), col=cm.colors(n), main="cm.colors()")

# ... or we could do something a bit more fancy: here is
# code that will generate a spectrum loosely based on the Munsell
# color wheel of perceptually equidistant colors (with some
# tweaks to increase the spread in the blue-green while
# keeping within the display gamut).

eqSpect <- colorRampPalette(
        c("#f2003c", "#F66900", "#F19100", "#F1B100",
          "#EFD300", "#f0ea00", "#CBDB00", "#9DD501",
          "#5ED108", "#00AF63", "#00A78E", "#00a3ac",
          "#0093af", "#0082b2", "#006ebf", "#4F37C2",
          "#8000D3", "#A001BF", "#BA00A4", "#D0007F"),
        space="rgb",
        interpolate="spline")

# Such a perceptually tuned spectrum is quite a bit more pleasing
# than one that is computed from extrapolating between rgb values.
# Adjacent colors are easier to distinguish, in particular hues that
# are close to the primary reds, greens and blues...

oPar <- par(mfrow=c(2,1))
pie(rep(1,n), col=rainbow(n), main="rainbow()")
pie(rep(1,n), col=eqSpect(n), main="eqSpect()")             
par(oPar)

# ==================================================
# Coloring tsne results
# ==================================================

# Let's use our eqSpect() colors to plot our tsn-embedded points,
# and color them according to their AP class:

plot(tsneRef[,1], tsneRef[,2], type='n');  # set up the frame

n <- length(apRes)
for (i in 1:n) {
        points(tsneRef[unlist(apRes[i]),1],
               tsneRef[unlist(apRes[i]),2],
               pch=21,
               col="#444444",
               bg= eqSpect(n)[i])
}

# Nice. We see that there is quite a bit of overlap between
# the clusters that we "see" in the embedding, and those that
# AP has defined. There are also clusters that look mixed up
# and should probably be joined.


# ==================================================
# Selecting from tsne results
# ==================================================

# Time for something a little bit more advanced. Looking
# at this plot I would like to know how the visual clusters
# are structured internally, i.e.  I would like to identify
# some points, find out what what their row-numbers are and
# plot them separately.

# The simplest approach is just to plot the row-numbers using
# the text() function.

plot(tsneRef[,1], tsneRef[,2], type='n');
text(tsneRef[,1], tsneRef[,2],
     labels=as.character(1:nrow(tsneRef)));

# ... well lets make the labels a bit smaller, for less overlap

plot(tsneRef[,1], tsneRef[,2], type='n');
text(tsneRef[,1], tsneRef[,2],
     labels=as.character(1:nrow(tsneRef)),
     cex=0.4);

# ... and color the labels according to the 16 AP clusters:
#
n <- length(apRes)
colV <- rep("", nrow(dat))
for (i in 1:n) {
        colV[as.integer(unlist(apRes[i]))] <- eqSpect(n)[i]
}
plot(tsneRef[,1], tsneRef[,2], type='n');
text(tsneRef[,1], tsneRef[,2],
     labels=as.character(1:length(tsneRef[,1])),
     cex=0.5, col=colV);

# We could now use this information and get a matplot() of the
# clusters by hand ... like so:
# First record  a set of datapoints we are interested in -
# from the lower left corner of the plot:
set1 <- c(5, 13, 44, 12, 80, 64, 16, 10, 48, 102, 31, 109, 42, 106, 147, 57, 85)
set2 <- c(189, 54, 35, 202, 119, 20, 50, 130, 161, 167, 8, 77, 236, 1)
set3 <- c(156, 187, 27, 223, 132, 224, 209, 73)


# Now open a separate window for the paralell coordinates plot
# so we can compare ...
w1 <- dev.cur()   # remember which one is our main device
dev.new()         # open a second graphics window
w2  <- dev.cur()  # store its ID

dev.set(w2)
matplot(t(datNorm[set1,]), 
        type="l", lwd=2, lty=1, col=eqSpect(16)[4], 
        xlab="time",ylab="log expression value")
for (i in 1:length(set2)) {
        lines(datNorm[set2[i], ], type="l", lwd=1, lty=1, col=eqSpect(16)[3])
}
for (i in 1:length(set3)) {
        lines(datNorm[set3[i], ], type="l", lwd=1, lty=1, col=eqSpect(16)[15])
}
dev.set(w1)

# Take some time to study this. It shows quite well how the internal
# structure within a set of similarly expressed genes is separated into
# clusters by AP, and how tsne gives us a good indication of that structure.

# ==================================================
# Interactive graphics with tsne results
# ==================================================

# But typing all those numbers by hand is tedious! Wouldn't it be nice if we
# could select them interactively? 

# Yes.

# R provides two functions for interactive work with 2D plots: 
# identify(), and locate().

# identify() returns the row number of a picked point.
# In order to pick points, we have to plot them
plot(tsneRef[,1], tsneRef[,2], pch=19, cex=0.8, col=colV); 

# Try this ...
identify(tsneRef)
# ... and press <esc> to stop (or a different mouse button)

# We can write a function to show us what we picked, display the
# actual profiles in a matplot() in our second window, and finally
# return the row numbers, so we can use them (e.g. to retrieve
# the gene symbols).

pickGenes <- function(x) {
        # picks genes from our tsne output plot.
        all <- c()
        dev.set(w2)
        matplot(t(datNorm[set1,]), type="n",
                xlab="time",ylab="log expression value")
        dev.set(w1)
        while(TRUE) {
                pick <- identify(x, n=1, plot=FALSE)
                if (!length(pick)) break
                points(x[pick,1], x[pick,2], pch=19, cex=0.8, col="#FFFFFFAA")
                print(pick)
                all <- c(all, pick)
                dev.set(w2)
                lines(datNorm[pick, ], type="l", lwd=1, lty=1, col=colV[pick])
                dev.set(w1)
        }
        return(all)
}

# Plot again ...
plot(tsneRef[,1], tsneRef[,2], pch=19, cex=0.8, col=colV); 

# ... and start picking. Press <esc> to stop.
myPicked <- pickGenes(tsneRef)
rownames(dat[myPicked,])

# ==================================================
# tsne 3D embedding
# ==================================================

# tsne can embed in 2- 3- or higher dimensions. There are
# inbuilt ways to work with 3D data in R - try ...
demo(persp)
# ... but for "real" work with 3D data, the standard
# is the rgl package. This may require you to update
# a few libraries on your machine ...

install.packages("rgl")
library(rgl)
# rgl has its own demo function that you could explore
# demo(rgl)

# here we simply plot coloured points in 3D by calling
# the plot3d() function as our callback for tsne() ...

plotTsne3D <- function(x){ 
        plot3d(x[,1], x[,2], x[,3], type='p', col=colV)
        #	text3d(x[,1], x[,2], x[,3], text=names(colV), col=colV, cex=1.4)
}

# ... and setting k=3 as the number of dimensions for tsne to embed
# into.

set.seed(112358)
tsne3D <- tsne(datNorm, epoch = 100, k=3, max_iter=3000, epoch_callback = plotTsne3D, perplexity=10)

# That's all.

# [End]