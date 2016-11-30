## After munge genetics has run. This is to convert and output a file that can be used in the Pritchard labs. The original data is DArT output and should look like the following:
#################
alist<-c("loci",185,186,187,188,189,190,191,"A549",1,1,1,1,1,1,1,"A549",0,0,1,1,1,0,1,"A588",1,1,1,1,1,1,1,"A588",0,0,0,0,0,0,1,"A794",1,1,1,1,1,1,1,"A794",1,0,1,0,1,1,0,"A081",1,1,1,1,1,1,0,"A081",1,1,1,1,1,1,1)
df <- data.frame(matrix(unlist(alist), nrow=9, byrow=T),stringsAsFactors=FALSE)
colnames(df) = df[1, ]
df<-df[-1, ]  
df
##################
##There are two ways of doing this. The first starts at line 36, the second at line 63. Lines 11-30 are just for munging my untidy data.

##Make like the original input.

X<-as.vector(DArTsnps$X.)
require(stringr)
X1<-substr(X,1,8)
X2<-gsub(pattern = "\\|", replacement = "",  x = X1)
X3<-sub(pattern = "CloneID", replacement = "loci", x = X2)


df<-DArTsnps
df$X.<-X3


df<-df[,-c(2:17)] # 17 cols of DArT calculations and allele sequences
df<-df[-c(1:4),] # 4 rows of DArT descriptors

write.table(df,file = "c:/temp/temp.txt",col.names = FALSE)
df<-read.table(file = "c:/temp/temp.txt", header = TRUE)
df$X5<-NULL


## First Option
# TRANSPOSING DATA FRAME
tdf <- as.data.frame(t(df[,-1]))

# SETTING COLUMN NAMES
names(tdf) <- as.list(df$loci)    
# SETTING INDIVIDUAL COLUMN
tdf$individual <- rownames(tdf)

# STACK SAME COLUMNS (CHANGE 8 TO NUMBER OF COLS(307))
finaldf <- rbind(tdf[, c(ncol(tdf), seq(1, 12728, 2))],   # EVEN COLS
                 tdf[, c(ncol(tdf), seq(2, 12728, 2))])   # ODD COLS

# ORDER BY INDIVIDUAL COLUMN
finaldf <- finaldf[with(finaldf, order(individual)), ]
rownames(finaldf) <- 1:nrow(finaldf)

finaldf[is.na(finaldf)] <- "-9"

# CONVERT LOCI COLUMNS TO NUMERIC
finaldf[,-1] <- sapply(sapply(finaldf[,-1], as.character), as.numeric)


write.table(finaldf,file = "C:\\hybridsStructureProj\\structure.txt",sep = "\t", row.names = FALSE)

##Remember to delete "individual" manually in structure text file. Unless you can work out how to remove it in r.

### Second solution

library(dplyr)
library(tidyr)
df2<-df %>% gather(individual, val, -loci) %>%
        group_by(loci, individual) %>%
        mutate(row = row_number()) %>%
        spread(loci, val) %>%
        select(-row)

write.table(df2,file = "c:\\structure.txt",sep = "\t", row.names = FALSE)
