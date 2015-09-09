# BIOCONDUCTOR
CHEATSHEET FOR HARVARD PH525.4X
#WEEK 1

QUESTION 1.1.1
table(expr.meta$gender)
110
QUESTION 1.1.2
summary(expr.meta$pkyrs)
40
QUESTION 1.1.3
qqnorm(expr.meta$pkyrs)
FALSE
QUESTION 1.2.1
seqs = c("ATG", "TGA", "TAA","TAG")
n =sapply(seqs, function(x) countPattern(x,chr11seq))
which.max(n)
TAA & 2624324
QUESTION 1.2.2
chr7seq = BSgenome.Hsapiens.UCSC.hg19$chr7
alphabetFrequency(chr7seq, as.prob=TRUE)
0.19901933
QUESTION 1.2.3  
s17$loc[which(s17$RefSNP_id=="73971683")]
##or with dplyr
library(dplyr)
s17 %>% filter(RefSNP_id=="73971683") %>% select(loc)
135246
QUESTION 1.3.1
boxplot(e["209169_at",]~tissue,las=2)
This gene is expressed in the brain but not the other tissues 
QUESTION 1.3.2  
IDs = c("201884_at","209169_at", "206269_at","207437_at","219832_s_at","212827_at")
library(rafalib)
mypar2(3,2)
for(i in IDs){
 boxplot(e[i,]~tissue,las=2)
"206269_at" 
QUESTION 1.3.3
hgu133a.db 
QUESTION 1.3.4
 In an eSet object with the e in the assayData accessible with exprs() and tissue as one of the columns in the phenoData
 

