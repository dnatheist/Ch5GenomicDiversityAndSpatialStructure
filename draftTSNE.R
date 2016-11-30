#tSNE
require(tsne)

#munge
df<-ellarvalPeeliSnps; main<-"2011 Larval snps"
df<-df[-c(1:4),]

pic<-df[2:nrow(df),14] #check pic in order
plot(as.character(pic),main="Polymorphic Information Content (PIC)", ylab="PIC")# to show PIC is in order (ish) ideally cut off would be 1300 or so.
rm(pic)

#But choose only first 800 so PC does not crash
#df<-df[-c(1500:nrow(df)),] #800 rows (they should be in order of PIC so the highest first)

df<-data.frame(t(df))
df<-df[-(1:17),]

row.names(df)<-df[,1]

df[,1]<-NULL

indx <- sapply(df, is.factor)
df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))

# df<-dist(df)
# df<-data.frame(as.matrix(df))

i <- (colSums(df, na.rm=T) != 0) # T if colSum is not 0
df<-df[,i]
df<-df[, !apply(is.na(df) | df == 0, 2, all)]

#tsne go
colors = rainbow(length(row.names(df)))
names(colors) = unique( row.names(df))
ecb = function(x,y){ plot(x,t='n'); text(x,labels= row.names(df)) }
tsne_df = tsne(df[,1:13], max_iter=1000, epoch_callback = ecb, perplexity=5)

#tsne go
colors = rainbow(length(row.names(df)))
names(colors) = unique( row.names(df))
ecb = function(x,y){ plot(x,t='n'); text(x,labels= row.names(df), col=colors[ row.names(df)]) }
tsne_df = tsne(df[,1:13], max_iter=1000, epoch_callback = ecb, perplexity=5)


########################################

colors = rainbow(length(unique(iris$Species)))
names(colors) = unique(df$Species)
ecb = function(x,y){ plot(x,t='n'); text(x,labels=row.names(df), col=colors[df$Species]) }
tsne_iris = tsne(iris[,1:4], epoch_callback = ecb, perplexity=50)

# compare to PCA
dev.new()
pca_iris = princomp(iris[,1:4])$scores[,1:2]
plot(pca_iris, t='n')
text(pca_iris, labels=iris$Species,col=colors[iris$Species])
