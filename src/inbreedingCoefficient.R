#Using adegenet rather than related to calculate Inbreeding Coefficient
#Something goes wrong because of large numbers?

inbreedCoeff<-inbreeding(gi, N=100)
meanInbreedCoeff <- sapply(inbreedCoeff, mean)
hist(meanInbreedCoeff, col="firebrick", main="Mean inbreeding - Murray Cod in Upper Murrumbidgee")

par(mfrow=c(2,3))
#By Site - not very revealing
inbreedCoeffTharwaSandwash<-seppop(gi)$'Tharwa Sandwash'
inbreedCoeffTharwaSandwash<-inbreeding(inbreedCoeffTharwaSandwash, N=100)
meanInbreedCoeffTharwaSandwash <- sapply(inbreedCoeffTharwaSandwash, mean)
hist(meanInbreedCoeffTharwaSandwash, col="firebrick", main="Mean inbreeding - Murray Cod at Tharwa Sandwash, Upper Murrumbidgee")

inbreedCoeffLanyon<-seppop(gi)$Lanyon
inbreedCoeffLanyon<-inbreeding(inbreedCoeffLanyon, N=100)
meanInbreedCoeffLanyon <- sapply(inbreedCoeffLanyon, mean)
hist(meanInbreedCoeffLanyon, col="firebrick", main="Mean inbreeding - Murray Cod at Lanyon, Upper Murrumbidgee")

inbreedCoeffMurramore<-seppop(gi)$'Murramore'
inbreedCoeffMurramore<-inbreeding(inbreedCoeffMurramore, N=100)
meanInbreedCoeffMurramore <- sapply(inbreedCoeffMurramore, mean)
hist(meanInbreedCoeffMurramore, col="firebrick", main="Mean inbreeding - Murray Cod at Murramore, Upper Murrumbidgee")

inbreedCoeffKambahPool<-seppop(gi)$'Kambah Pool'
inbreedCoeffKambahPool<-inbreeding(inbreedCoeffKambahPool, N=100)
meanInbreedCoeffKambahPool <- sapply(inbreedCoeffKambahPool, mean)
hist(meanInbreedCoeffKambahPool, col="firebrick", main="Mean inbreeding - Murray Cod at Kambah Pool, Upper Murrumbidgee")


inbreedCoeffBullenRange<-seppop(gi)$'Bullen Range'
inbreedCoeffBullenRange<-inbreeding(inbreedCoeffBullenRange, N=100)
meanInbreedCoeffBullenRange <- sapply(inbreedCoeffBullenRange, mean)
hist(meanInbreedCoeffBullenRange, col="firebrick", main="Mean inbreeding - Murray Cod at Bullen Range, Upper Murrumbidgee")


inbreedCoeffNerreman<-seppop(gi)$'Nerreman'
inbreedCoeffNerreman<-inbreeding(inbreedCoeffNerreman, N=100)
meanInbreedCoeffNerreman <- sapply(inbreedCoeffNerreman, mean)
hist(meanInbreedCoeffNerreman, col="firebrick", main="Mean inbreeding - Murray Cod at Nerreman, Upper Murrumbidgee")

#But see http://lists.r-forge.r-project.org/pipermail/adegenet-forum/2014-May/000836.html

gi
res <- inbreeding(gi, res.type="function")
par(mfrow=c(8,10), mar=rep(.1,4))
for(i in 1:80) plot(res[[i]], yaxt="n", col="red")

par(mfrow=c(1,1)
plot(res[[1]])
#the y axis is sooooo small sampling fails in the first method. The reference link suggest res.type  = "function" is best. Interestingly the mean(meanInbreedCoeff) is 0.5021
