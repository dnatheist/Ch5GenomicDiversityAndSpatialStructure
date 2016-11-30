## @knitr simulation

# This script allows running of 'related' package from a two row DArT file. This is an R package that allows users to estimate pairwise relatedness for pairs of individuals based on codominant molecular markers (microsatellites, SNPs, etc.), and also has simulation capabilities for comparing the performance of different estimators and for testing the resolution of a data set. It can also test for relatedness patterns within groups. Relatedness can be estimated using any of seven different methods, and can incorporate inbreeding and genotyping errors. The underlying code for implementing these relatedness estimators is Jinliang Wang’s Fortran code for his COANCESTRY program (Wang 2011).
# Alan Couch, Institute for Applied Ecology, University of Canberra

####################
# Note that one needs to run 32 bit version of R
# (set this in Global options under Tools menu)
# because 'related' relies on some Fortran code that does not run under 64 bit R.
####################
 

#Simulation code direct from 'related' tutorial paper
fileName=paste0("./outData/",main," famSimInput",sep=" ")
load(file = fileName)
sim <- familysim(input$freqs,100)
outputSim <- coancestry ( sim , wang=1)
##output2 <- coancestry ( sim , trioml=1) # takes 22 hours so make an export for Bernd
simrel <- cleanuprvals (outputSim$relatedness , 100)
fileName=paste0("outData/",main," outputSim",sep=" ")
save(simrel,file = fileName)
relvalues <- simrel [, 6] #5 for trioml, 6 for wang, 7 for lynchli,8 lynchrd, 9 ritland,10quellergt, 11 dyadml
label1 <- rep ("PO", 100)
label2 <- rep (" Full ", 100)
label3 <- rep (" Half ", 100)
label4 <- rep (" Unrelated ", 100)
labels <- c( label1 , label2 , label3 , label4 )


#Plot as a density plot
Relationship <- labels
newdata <- as.data.frame ( cbind ( Relationship , relvalues ) )
newdata$relvalues <- as.numeric ( as.character ( newdata$relvalues ))
qplot ( relvalues , ..density.. , data = newdata , geom ="density", colour = as.factor ( Relationship ) , xlab =" Relatedness Value ", ylab ="Density",main=main)


