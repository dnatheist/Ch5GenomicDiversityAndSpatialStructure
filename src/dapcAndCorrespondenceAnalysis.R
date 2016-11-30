## @knitr dapcCorroStructure

#Looking for DAPC, fixed differences, and Correspondence Analysis between populations.

require(reshape2)
require(dart)
require(pegas)
all.dart <- read.dart("otherData/larvalPeeliSnps.csv", topskip = 5)
gl.dart <- dart2genlight(all.dart, covfilename = "otherData/qslLarvalPeeliiMetaForPCA.csv") # this glObject is then suitable for sharing via GitHub or such for open data as required by PeerJ etc. Need better file but.

# cant work gl.fixed.diff(gl.dart, t=0)

#PCA
pca1 <- glPca(gl.dart, nf = 3 ,parallel = FALSE)
pca1

scatter(pca1, posi="topleft")
title("PCA of Murray Cod data\n axes 1-2")

myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)

require(ape)
tre <- nj(dist(as.matrix(gl)))
tre

plot(tre, typ="fan", cex=0.7)
title("NJ tree of the Upper Murrumbidgee Murray Cod data")

plot(tre, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=myCol, cex=4)
title("NJ tree of the Upper Murrumbidgee Murray Cod data")

#Discriminant Component Analysis 
dapc2 <- dapc(gl.dart, n.pca=10, n.da=1, parallel=FALSE)
scatter(dapc2,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE, txt.leg=paste("group", 1:5), col=c("red","blue","green","orange", "pink"))
compoplot(dapc2, col=c("red","blue","green","orange", "pink"),lab="", txt.leg=paste("group", 1:5), ncol=5)




#Analysis on genInd object (individuals)
#Hardy Weinberg
source("genlight2genInd.R")
gI<-genlight2genind(gl.dart)
hwe<-hw.test(gI, B=0)
hwedf<-data.frame(hwe)
hist(hwedf$Pr.chi.2...)


#There are 
sum(is.na(gI$tab)) #missing values in the data
forPCOA<-scaleGen(gI, NA.method="mean") #scale the data
forPCOA[1:5,1:5]

pca1 <- dudi.pca(forPCOA,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
s.label(pca1$li)
title("PCA of cod dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
s.class(pca1$li, pop(gI))

s.class(pca1$li,pop(gI),xax=1,yax=3,sub="PCA 1-3",csub=2)
title("PCA of cod dataset\naxes 1-3")
add.scatter.eig(pca1$eig[1:20],nf=3,xax=1,yax=3)

col <- funky(15)
s.class(pca1$li, pop(gI),xax=1,yax=3, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE)

colorplot(pca1$li[c(1,3)], pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
title("PCA of cod dataset\naxes 1-3")
abline(v=0,h=0,col="grey", lty=2)

#Correspondance Analysis (genepop object)
obj <- genind2genpop(gI)
ca1 <- dudi.coa(tab(obj),scannf=FALSE,nf=3)
barplot(ca1$eig,main="Correspondence Analysis eigenvalues",
        col=heat.colors(length(ca1$eig)))

s.label(ca1$li, sub="Correspondance Analysis",csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=1,yax=2,posi="topleft")

#Isolation By Distance
# Dgen <- dist.genpop(obj,method=2)
# Dgeo <- dist(gI$other$latlong) #need xys for sites in gI for this to work.
# ibd <- mantel.randtest(Dgen,Dgeo)
# ibd





