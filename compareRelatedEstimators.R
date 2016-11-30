## @knitr compareRelEstimators

fileName=paste0("./outData/",main," famSimInput",sep=" ")
load(file = fileName) #this is the right file (familySim and CE can use same input file)
compareestimators (input , 100)

#Compare estimators with each other 
# plot(output$relatedness[,5:11]) # the wang estimator seems to give best correlation with the ML estimators.
# 
# plot(output$inbreeding[,2:3]) # corelation between inbreeding estimators. Should include the ML methods but does not output them for some reason despite the doco saying it should.
# 
# #Histogram of inbreeding coefficients
# hist(output$inbreeding$LH)
# hist(output$inbreeding$LR)

# # At some point this will be better to be change a la https://frasierlab.files.wordpress.com/2015/03/tutorial1.pdf custom comparisons (pg 16)
# 5.3.2 Custom Comparisons
# It is also relatively easy to compare a different set of relatedness estimators. Suppose you want to compare
# the two likelihood estimators with those of Queller & Goodnight and Wang. First, load the data into R and
# generate as many simulated individuals as you like (here we’ll use 100).
# input <- readgenotypedata (" GenotypeData . txt ")
# simdata <- familysim ( input$freqs , 100)
# Then, estimate relatedness using the desired estimators (note that this will take a long time if using the
#                                                          likelihood estimators).
# output <- coancestry ( simdata , dyadml =1 , trioml =1 , quellergt =1 , wang =1)
# This file will contain all pairwise estimates of relatedness, rather than just the ones that we are interested
# in. To reduce this to only the desired values, use the cleanuprvals function.
# simrel <- cleanuprvals ( output$relatedness , 100)
# Next, we need to parse out the data based on relatedness type and estimator used. In the file, the first
# 100 pairs are parent-offspring, the second set of 100 are full-sibs, the third set of 100 are half-sibs, and the
# fourth (last) set of 100 are unrelated. Obviously, these numbers would change if you simulated a different
# 16
# number of individuals. Moreover, the trioml estimates are in column 5, the wang estimates are in column
# 6, the quellergt estimates are in column 10, and the dyadml estimates are in column 11 (see table on p. 8).
# What the code below is doing is selecting the range of rows and columns that correspond to the appropriate
# relatedness value and estimator.
# triomlpo <- simrel [1:100 , 5]
# triomlfs <- simrel [(100 + 1) : (2 * 100) , 5]
# triomlhs <- simrel [((2 * 100) + 1) : (3 * 100) , 5]
# triomlur <- simrel [((3 * 100) + 1) : (4 * 100) , 5]
# wangpo <- simrel [1:100 , 6]
# wangfs <- simrel [(100 + 1) : (2 * 100) , 6]
# wanghs <- simrel [((2 * 100) + 1) : (3 * 100) , 6]
# wangur <- simrel [((3 * 100) + 1) : (4 * 100) , 6]
# quellergtpo <- simrel [1:100 , 10]
# quellergtfs <- simrel [(100 + 1) : (2 * 100) , 10]
# quellergths <- simrel [((2 * 100) + 1) : (3 * 100) , 10]
# quellergtur <- simrel [((3 * 100) + 1) : (4 * 100) , 10]
# dyadmlpo <- simrel [1:100 , 11]
# dyadmlfs <- simrel [(100 + 1) : (2 * 100) , 11]
# dyadmlhs <- simrel [((2 * 100) + 1) : (3 * 100) , 11]
# dyadmlur <- simrel [((3 * 100) + 1) : (4 * 100) , 11]
# Next, we need to create a list of labels for the different estimators, with each repeated the appropriate number
# of times (100 in this case).
# trioml <- rep ("tri", 100)
# wang <- rep ("W", 100)
# quellergt <- rep ("QG", 100)
# dyadml <- rep ("di", 100)
# estimator2 <- c( trioml , wang , quellergt , dyadml )
# Estimator <- rep ( estimator2 , 4)
# Create a list of labels for the different relatedness types
# po <- rep (" Parent - Offspring ", (4 * 100) )
# fs <- rep ("Full - Sibs ", (4 * 100) )
# hs <- rep ("Half - Sibs ", (4 * 100) )
# ur <- rep (" Unrelated ", (4 * 100) )
# relationship <- c(po , fs , hs , ur )
# Combine the different values for each estimator based on relatedness type, as lists.
# relatednesspo <- c( triomlpo , wangpo , quellergtpo , dyadmlpo )
# relatednessfs <- c( triomlfs , wangfs , quellergtfs , dyadmlfs )
# relatednesshs <- c( triomlhs , wanghs , quellergths , dyadmlhs )
# relatednessur <- c( triomlur , wangur , quellergtur , dyadmlur )
# Relatedness_Value <- c( relatednesspo , relatednessfs , relatednesshs , relatednessur )
# Combine the data.
# combineddata <- as . data . frame ( cbind ( Estimator , relationship , Relatedness_Value ))
# combineddata$Relatedness_Value <- as . numeric ( as . character (
#         combineddata$Relatedness_Value ))
# Plot the data. You may need to play around with the range of values of the y-axis (the ylim argument) so
# that it better captures your data.
# ggplot ( combineddata , aes ( x = Estimator , y = Relatedness_Value ) , ylim = c ( -0.5 ,
#                                                                                      1.0) ) +
#         geom_boxplot () +
#         facet_wrap (~ relationship )
# 17
# The resulting plot is shown below. Note that you are not limited to comparing only four estimators. You
# can compare as many as you wish, just make sure that the commands are changed accordingly.
