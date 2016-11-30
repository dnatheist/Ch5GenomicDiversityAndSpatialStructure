# make a larval M. peelii PCA plot to seperate the outbreds

require(dart)
all.dart <- read.dart("otherData/larvalPeeliSnps.csv", topskip = 5)
gl.dart <- dart2genlight(all.dart, covfilename = "otherData/qslLarvalPeeliiMetaForPCA.csv") # this glObject is then suitable for sharing via GitHub or such for open data as required by PeerJ etc. Need better file but.

#PCA plot
gl <- gl.dart
system.time(pca1 <- glPca(gl[,], parallel=F, nf=3))
require(ggplot2)
gg <- data.frame(pca1$scores, cluster=gl$other$covariates$raceCladeName)   #cluster=pop(gl)) or cluster=gl$other$covariates$YearOnly)
gg$labels  <- rownames(gg)

ggplot(gg, aes(x=PC1, y=PC2)) +
        geom_text(aes(label=labels, color=factor(cluster)),hjust=-0.2, size=8) +
        geom_hline(yintercept=0,linetype=2) +
        geom_vline(xintercept=0,linetype=2) +
        scale_color_discrete(name="Cluster") + 
        ggtitle("Principal Components - Maccullochella peelii ") +
        theme(plot.title = element_text(lineheight=.8, face="bold", size=29)) +
        theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),legend.text = element_text(size = 20),legend.title = element_text(size = 20))