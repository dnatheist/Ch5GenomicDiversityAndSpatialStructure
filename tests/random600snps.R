## @knitr random600

#Randomly sample pairs of alleles from DArTsnps
#May wish to use 'df' rather than 'DArTsnp's and then insert in "runRelated.RMD" after loading data frame. NB: though these have an additional column 'allele' up front and so program below would need adjusting. 

set.seed(1) #in case we need to repeat
snp1<-sample(6:nrow(df),350) #get 50 extra to allow for removal of duplicate allelepairs later
#snp1<-as.integer(runif(300,18+1,nrow(df)-1)) #random sample from uniform distribution. Include -1 and +1 just to ensure it does not break if runif selects last or first row by (small) chance.

#grab the correct pair to the first random partner selected thus making 600 snps
snp2<-ifelse(snp1 %% 2 == 0, snp1+1, snp1-1)

allelePair<-c(snp1,snp2) #join the snp pairs up
allelePair<-sort(allelePair)# sort them
allelePair<-unique(allelePair)# remove any dups that may have arisen when concatenating snp1 and 2
hist(allelePair) #check them

DArTsnps600<-df[allelePair,] # grab the corresponding pairs of DArT snps.
df<-rbind(df[1,],DArTsnps600)
