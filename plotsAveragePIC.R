## @knitr averagePIC

#A script to plot all M. peelii snps average PIC (PIC of reference and PIC of alternate state allele) against index number of the data frame. This shows that using the first X number of snps is using the alleles with the highest PIC.


df<-allPeeliSnps; main<-"All Peelii snps"
df<-df[-c(1:4),]

pic<-df[2:nrow(df),14] #check pic in order
plot(as.character(pic),main="Polymorphic Information Content (PIC)", ylab="PIC")# to show PIC is in order (ish) ideally cut off would be 1300 or so.
smoothScatter(pic,main="Polymorphic Information Content (PIC)", ylab="Average PIC")

# Place in strict order of PIC and have  alook
pico<-data.frame(pic)
rownames(pico) <- NULL
colnames(pico)[1]<-"PIC"

pico <- data.frame(pico[order(as.numeric(pico$PIC),decreasing = TRUE),])

pico[,2]<-pico[,1] #Order columns to facilitate plotting.
pico[,1]<-row.names(pico) #Replace name with index number
colnames(pico)[1]<-"Index"
colnames(pico)[2]<-"PIC"

plot(pico, main="Ordered Polymorphic Information Content (PIC) ")
smoothScatter(pico, main="Ordered Polymorphic Information Content (PIC)")

rm(pic)
rm(pico)