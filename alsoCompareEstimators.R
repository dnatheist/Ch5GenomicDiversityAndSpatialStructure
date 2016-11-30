#Compare estimators with each other 
plot(output$relatedness[,5:11]) # the wang estimator seems to give best correlation with the ML estimators.

plot(output$inbreeding[,2:3]) # corelation between inbreeding estimators. Should include the ML methods but does not output them for some reason despite the doco saying it should.

#Histogram of inbreeding coefficients
hist(output$inbreeding$LH)
hist(output$inbreeding$LR)
