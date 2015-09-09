# BIOCONDUCTOR

CHEATSHEET FOR HARVARD PH525.4X

#WEEK 1

QUESTION 1.1.1
table(expr.meta$gender)
"110"

QUESTION 1.1.2
summary(expr.meta$pkyrs)
"40"

QUESTION 1.1.3
qqnorm(expr.meta$pkyrs)
"FALSE"

QUESTION 1.2.1
seqs = c("ATG", "TGA", "TAA","TAG")
n =sapply(seqs, function(x) countPattern(x,chr11seq))
which.max(n)
"TAA & 2624324"

QUESTION 1.2.2
chr7seq = BSgenome.Hsapiens.UCSC.hg19$chr7
alphabetFrequency(chr7seq, as.prob=TRUE)
"0.19901933"

QUESTION 1.2.3  
s17$loc[which(s17$RefSNP_id=="73971683")]
##or with dplyr
library(dplyr)
s17 %>% filter(RefSNP_id=="73971683") %>% select(loc)
"135246"

QUESTION 1.3.1
boxplot(e["209169_at",]~tissue,las=2)
"This gene is expressed in the brain but not the other tissues"

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
"a SNP (single nucleotide polymorphism)"

Q 1.4.4
The existence or nonexistence of the bud is a low-dimensional characteristic of the organism.
"phenotypic states phenotypic states"

QUESTION 1.5.1
length(unique(keys(Homo.sapiens, keytype="ENTREZID")))
"47721"

QUESTION 1.5.2
length(unique(keys(Homo.sapiens, keytype="ENSEMBL")))
"28553"

QUESTION 1.5.3
select(Homo.sapiens, key="9575", keytype="ENTREZID", columns=c("SYMBOL", "ENSEMBL", "ENTREZID", "CHR"))
"ENSG00000134852"

QUESTION 1.5.4
length(unique(tab$ENTREZID))
"56"

QUESTION 1.6.1
female = samp[,samp$sex == "Female"]
sum(exprs(female[1,]))
"1231.46"

QUESTION 1.6.2
experimentData(samp)
"pfermat@lab.not.exist"

QUESTION 1.6.3
annotation(samp)

library(BiocInstaller)
biocLite("hgu95av2.db")

library(hgu95av2.db)
fns = featureNames(samp)
annot = select(hgu95av2.db, keys=fns, keytype="PROBEID", columns="SYMBOL")
## this map is not one to one, so pick one:
geneSymbols = annot [ match(fns, annot$PROBEID), "SYMBOL"]
"hgu95av2"

QUESTION 1.6.4
cor(samp$score, exprs(samp)["31489_at",])
"0.1376892"


##WEEK 2


QUESTION 2.1.1
biocVersion()
"3.0"

QUESTION 2.1.2
library(ERBS)
data(HepG2)
class(HepG2)
"GRanges"

class(HepG2)
##More specifically you can see it here:
attributes(class(HepG2))
help("GRanges-class")
"GenomicRanges"

QUESTION 2.1.3
##You can see it by just typing and reading the provided information
HepG2
## or you can type
length(HepG2)
"303"

QUESTION 2.2.1
library(ERBS)
data(HepG2)
median( mcols(HepG2)$signalValue )
"7.024"

QUESTION 2.2.2
i=which.max(  mcols(HepG2)$signalValue  )
seqnames(HepG2[i])
"chrX"

QUESTION 2.2.3
chr = seqnames(HepG2)
table(chr)[16]
"31"

QUESTION 2.2.4
median( width(HepG2) )
##You can see the histogram
hist(width(HepG2),nclass=25)
"560"

QUESTION 2.3.1
IRanges(101,200) * 2
"126"

QUESTION 2.3.2
narrow(IRanges(101,200), start=20)
"120"

QUESTION 2.3.3
IRanges(101,200)+25
"150"

QUESTION 2.3.4
sum(width(IRanges(start=c(1,11,21),end=c(3,15,27))))
"15"

QUESTION 2.3.5
range(x)
sum(width(gaps(x)))
"130"

QUESTION 2.3.6
length(disjoin(x))
"17"

QUESTION 2.3.7
From the man page: resize resizes the ranges to the specified width where either the start, end, or center is used as an anchor. The default is fix="start", so resize(x,1) gives you the starting integer of each range in x.

"It gives you just the starting point of each range. It gives you just the starting point of each range."




