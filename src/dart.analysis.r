source("src/generic/dart.r")
source("src/generic/read.dart.r")

#library(devtools)
#Install non-cran version (the older version of adgenet has a serious flaw and is not updated yet.)
#install_github("thibautjombart/adegenet")
#install_github("green-striped-gecko/PopGenReport")


all.dart <- read.dart("OtherData/goodDArTsnps.csv", topskip = 5)
gl.dart <- dart2genlight(all.dart, covfilename = "OtherData/qslDartCovariatesGoodSpecies.csv")

#======================================================
#FST

library(StAMPP)
system.time(snpfst <-stamppFst(gl.dart,nboots=1, percent=95, nclusters=8 ))


#line 58 and 59 in munge code file causes clustering and pcato break
#===========================================================
# K means clustering
system.time(grp <- find.clusters(gl.dart, n.pca=100, choose=TRUE, stat="BIC", parallel=F))
plot(grp$Kstat, type="o", xlab="number of clusters (K)",
     ylab="BIC",main="find.clusters on a genlight object\n(chose n groups)")
#========================
#DAPC see: http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
dapc1<-dapc(gl.dart, pop=NULL, n.pca=NULL, n.da=NULL,
scale=FALSE, var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
pca.select=c("nbEig", "percVar"), perc.pca=NULL, glPca=NULL,parallel=FALSE)
scatter(dapc1)
myCol <- c("darkblue","purple","green","orange","red","blue")
scatter(dapc1, posi.da="bottomright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft")
#or
scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:6))
#or
scatter(dapc1, ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=myCol, solid=.4, cex=3, clab=0,
        mstree=TRUE, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, txt.leg=paste("Cluster",1:6))
par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=myCol)

myInset <- function(){
        temp <- dapc1$pca.eig
        temp <- 100* cumsum(temp)/sum(temp)
        plot(temp, col=rep(c("black","lightgrey"),
                           c(dapc1$n.pca,1000)), ylim=c(0,100),
             xlab="PCA axis", ylab="Cumulated variance (%)",
             cex=1, pch=20, type="h", lwd=2)
}
add.scatter(myInset(), posi="bottomright")#,inset=c(-0.03,-0.01), ratio=.28, bg=transp("white"))

#Assess which alleles responsible for variation
set.seed(4)
contrib <- loadingplot(dapc1$var.contr, thresh=0.00048, axis=2,lab.jitter=1)

#=============================================================================
#PCA plot

###out to a win metafile
gl <- gl.dart
system.time(pca1 <- glPca(gl[,], parallel=F, nf=3))
win.metafile("Maccullochella.wmf", width = 17, height=17)
plot(pca1$scores, pch=16, col=pop(gl), asp=1)
text(pca1$scores, labels=gl@ind.names, col=as.numeric(pop(gl)), cex=0.8)
#points(pca1$scores, pch=16, col=rainbow(274)
dev.off()

### Output ggplot a bit nicer
require(ggplot2)
gg <- data.frame(pca1$scores,cluster=pop(gl))
gg$labels  <- rownames(gg)
win.metafile("maccullochella2.wmf", width = 17, height=17)
ggplot(gg, aes(x=PC1, y=PC2))+   geom_text(aes(label=labels, color=factor(cluster)),hjust=-0.2) +     geom_hline(yintercept=0,linetype=2) +     geom_vline(xintercept=0,linetype=2) +    scale_color_discrete(name="Cluster") + ggtitle("Principal Component SNPS ")+ theme(plot.title = element_text(lineheight=.8, face="bold", size=20)) + annotate("text", x = 3, y = 0.7, label = "Murray Cod") + annotate("text", x = 70, y = -0.3, label = "Trout Cod") + annotate("text", x = 18, y = 0.3, label = "F1 Backcross") + annotate("text", x = 35, y = 0.1, label = "F1 Hybrid") +xlim(-5,80)
dev.off()

###Also this which is less nice pcoa

system.time(pca1 <- glPca(gl.dart  , parallel=F, nf=3))
s.class(pca1$scores, pop(gl.dart), col=rainbow(nlevels(pop(gl.dart))))

scatter(pca1, ratio=.2)
barplot(pca1$eig, main="eigenvalues", col=heat.colors(length(pca1$eig)))

s.class(pca1$scores, pop(gl.dart), col=colors()[c(131,134)])
add.scatter.eig(pca1$eig,2,1,2)

###or 3d
library(rgl)
plot3d(pca1$scores[,1:3],col="red")

#====================================================================================
#To convert a genlight to genind (needed to then convert to Pritchards 'structure' input format)
gi.dart <- genlight2genind(gl.dart)

#% !Rnw weave = Sweave
library(PopGenReport)
pgr <- popgenreport(gi.dart, mk.complete=TRUE )
