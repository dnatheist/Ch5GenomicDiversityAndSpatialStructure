# Clear the workspace, to be safe
rm(list=ls())
# load support functions
source("c:/R_generic_scripts/gl.read.dart.r")
source("c:/R_generic_scripts/is.fixed.r")
source("c:/R_generic_scripts/gl.fixed.diff.r")
source("c:/R_generic_scripts/gl.report.CallRate.r")
source("c:/R_generic_scripts/gl.report.RepAvg.r")
source("c:/R_generic_scripts/gl.filter.CallRate.r")
source("c:/R_generic_scripts/gl.filter.RepAvg.r")
source("c:/R_generic_scripts/gl2nhyb.r")

library(reshape2)
library(MASS)

# Set working directory
setwd("c://R_cod")

# READ DATA TO ADEGENET GENLIGHT FORMAT
gl<-gl.read.dart(datafile="allDArTsnps.csv", topskip=5, nmetavar=17, nas="-", ind.metafile="pop_assign_allDArTsnps.csv")
# Data identified for  285 individuals,  6364 loci

pop(gl)


  # Check the callrate
  gl.report.CallRate(gl)
  
  # FILTER STRINGENTLY ON CALLRATE
   gl <- gl.filter.CallRate(gl,0.983)

  # cHECK THE REPRODUCIBILITY
   gl.report.RepAvg(gl)
   
  # FILTER STRINGENTLY ON REPRODUCIBILITY
  gl <- gl.filter.RepAvg(gl,1)
  
  # DEFINE THE PARENTAL POPULATIONS
  p0 <- c("A21","A22","A23","A24","A37","A29","A30","A31","A38","A47","A39","A40","A41","A42","A43","A44","A45","A46","A26","A36","A33","A34","A35","A48","A28" )
  p1 <- c("Bend","Narra")

  # FILTER OUT FROM THE MASTER FILE, THOSE LOCI THAT ARE NOT FIXED BETWEEN PARENTALS, AND OUTPUT TO NEWHYBRIDS FORMAT
  gl.out <- gl2nhyb(gl, outfile="codtest.txt", p0, p1, threshold=0)
