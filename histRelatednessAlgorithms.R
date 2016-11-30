load("output 3 Years Larval snps")

plot(output$relatedness$pair.no,output$relatedness[,5])

hist(output$relatedness$trioml, main="trioml")
hist(output$relatedness$wang, main="wang")
hist(output$relatedness$lynchli, main="lynchli")
hist(output$relatedness$lynchrd, main="lynchrd")
hist(output$relatedness$ritland, main="ritland")
hist(output$relatedness$quellergt, main="quellergt")
hist(output$relatedness$dyadml, main="dyadml")

summary(output$relatedness$dyadml)

