## @knitr PCAtSNEStructure

#First we do PCA and tSNE using peelii larvae and site
#Vanilla PCA
require(adegenet)
require(ggplot2)

system.time(pca1 <- glPca(gl[,], parallel=F, nf=3))
gg <- data.frame(pca1$scores,cluster=pop(gl))
gg$labels  <- rownames(gg)

ggplot(gg, aes(x=PC1, y=PC2))+
        geom_text(aes(label=labels, color=factor(cluster)),hjust=-0.2, size=4) +
        geom_hline(yintercept=0,linetype=2) +
        geom_vline(xintercept=0,linetype=2) +
        scale_color_discrete(name="Site Name") +
        ggtitle("Principal Components - Peelii larvae and sites")+
        theme(plot.title = element_text(lineheight=.8, face="bold", size=29))+
        theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),legend.text = element_text(size = 20),legend.title = element_text(size = 20)) 

#tsne 2d
require(tsne)
colors<- rainbow(length(unique(gg$cluster)))
names(colors) <- unique(gg$cluster)
#ecb <-function(x,y){ plot(x,t='n',main="Larval M. peelii by Sites"); text(x,labels= row.names(gg),col=colors[gg$cluster])} #to make series
#tsne_gg <- tsne(gg[,1:3], max_iter=10000,k=2, epoch_callback = ecb, perplexity=5,epoch=1000)
tsne_gg <- tsne(gg[,1:3], max_iter=10000,k=2, perplexity=5)
plot(tsne_gg,t='n',main="Larval M. peelii by Sites"); text(tsne_gg,labels= row.names(gg),col=colors[gg$cluster])

#Second we do PCA and tSNE with the same larvae coloured by year 
gg<- data.frame(pca1$scores,cluster=as.character(gl@other$covariates))
gg$labels  <- rownames(gg)
ggplot(gg, aes(x=PC1, y=PC2))+
        geom_text(aes(label=labels, color=factor(cluster)),hjust=-0.2, size=4) +
        # geom_hline(yintercept=0,linetype=2) +
        # geom_vline(xintercept=0,linetype=2) +
        scale_color_discrete(name="Year") +
        ggtitle("Principal Components - Peelii larvae and Years")+
        theme(plot.title = element_text(lineheight=.8, face="bold", size=29))+
        theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),legend.text = element_text(size = 20),legend.title = element_text(size = 20))

#tsne 2d
require(tsne)
colors<- rainbow(length(unique(gg$cluster)))
names(colors) <- unique(gg$cluster)
#ecb <-function(x,y){ plot(x,t='n',main="Larval M. peelii by Years"); text(x,labels= row.names(gg),col=colors[gg$cluster])} #to make series 
#tsne_gg <- tsne(gg[,1:3], max_iter=10000, epoch_callback = ecb, perplexity=5,epoch=1000)
tsne_gg <- tsne(gg[,1:3], max_iter=10000, perplexity=5)
plot(tsne_gg,t='n',main="Larval M. peelii by Years"); text(tsne_gg,labels= row.names(gg),col=colors[gg$cluster])

# ####################################Incase we need 3d plots
# require(scatterplot3d)
# tsne_gg <- tsne(gg[,1:3], max_iter=1000,k=3, perplexity=5)
# scatterplot3d(tsne_gg[,1],tsne_gg[,2],tsne_gg[,3], pch=16, highlight.3d=TRUE,  main="3D Scatterplot")
# 
# require(rgl)
# #plot3d(tsne_gg[,1],tsne_gg[,2],tsne_gg[,3], col="red", size=3)
# plot3d(tsne_gg[,1],tsne_gg[,2],tsne_gg[,3], col=colors[gg$cluster], size=3)





