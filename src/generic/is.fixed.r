
is.fixed<-function(s1, s2){
# A function to determine if 2 populations have a fixed difference at a given locus
  
# Takes as input two values allele frequencies to compare.
# Returns TRUE or FALSE, or NA if one or both of the values are missing or undefined.
  
# A fixed difference occurs when an allele is present in all individuals of one population and absent in the other.  
  
  if (is.na(s1) | is.na(s2)) {
    result<-NA}
  else{
    result<-0
    if ((s1==0) & s2==100) {result<-1}
    if (s1==100 & s2==0) {result<-1}
  }
  return(result)
}

# Test function
# is.fixed(0,100)
# is.fixed(100,0)
# is.fixed(80,0)
# is.fixed(100,NA)
# is.fixed(0,NA)
# is.fixed(NA,0)
# is.fixed(NaN,100)
# is.fixed(NaN,0)
# is.fixed(100,NaN)
# is.fixed(0,NaN)
# is.fixed(NaN,NaN)
# is.fixed(NA,NA)
