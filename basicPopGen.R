#gl to gind

x.mat <- as.matrix(gl) # x is a genlight object
x.mat[x.mat == 0] <- "1/1" # homozygote reference
x.mat[x.mat == 1] <- "1/2" # heterozygote
x.mat[x.mat == 2] <- "2/2" # homozygote alternate
x.gid <- df2genind(x.mat, sep = "/", ploidy = 2)
gi<-x.gid

gi
div <- summary(gi)
div
names(div)

plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")

gi2 <- genind2hierfstat(gi)

nAll(gi) # Number of alleles per locus

plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", main="Expected heterozygosity as a function of observed heterozygosity per locus")

bartlett.test(list(div$Hexp, div$Hobs)) # a test : H0: Hexp = Hobs

basicstat <- basic.stats(gi2, diploid = TRUE, digits = 2) 
names(basicstat) 

hw.test(gi, B = 1000)
