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
 
Q 1.4.1 
The number of copies of RNA depends on how much the gene is being transcribed. “Housekeeping genes” such as those used to make ribosomal or transfer RNA are transcribed at a high rate, while others, such as mRNA for some transcription factors, are transcribed less frequently. Some genes are not transcribed at all in certain cell types.
"VARIES"

Q 1.4.2
Humans are a diploid species, meaning the somatic cells typically contain two copies of the autosomal (not X or Y) chromosomes, before S phase in which the chromosomes are duplicated.
"2 copies"

Q 1.4.3
Genetic information, such as single nucleotide polymorphisms, have a chance of being transmitted across generations. mRNA transcripts and proteins such as transcription factors degrade over time and, most importantly, do not replicate themselves. In other words, DNA (not proteins or RNA) is known as the main molecule of genetic inheritance. (Side note: There are cases of mRNA and proteins being temporarily inherited, for example the mRNA which are in the egg cell at the moment it is fertilized by the sperm cell.)
"a SNP (single nucleotide polymorphism) a SNP (single nucleotide polymorphism)"

Q 1.4.4
The existence or nonexistence of the bud is a low-dimensional characteristic of the organism.
"phenotypic states phenotypic states"



