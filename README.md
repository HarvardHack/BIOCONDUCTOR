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

QUESTION 2.4.1  
par(mfrow=c(2,1))
plotGRanges(x)
plotGRanges(resize(x,1))
"at the left-most point of the "+" strand ranges in x, and the right-most point of the "-" strand ranges in x"

QUESTION 2.4.2
disjoined = disjoin(c(x,y))
in.both = disjoined %over% x & disjoined %over% y
sum(width(disjoined[ in.both ]))
"140"

QUESTION 2.4.3
disjoined = disjoin(c(x,y))
not.in.both = !(disjoined %over% x & disjoined %over% y)
sum(width(disjoined[ not.in.both ]))
"130"

QUESTION 2.4.4
z = GRanges("chr1", range(ranges(x)), strand="-")
sum(x %over% z)
"0"

QUESTION 2.5.1  
start(HepG2)[17]
"46528596"

QUESTION 2.5.2
start(HepG2[17])
"46528596"

QUESTION 2.5.3
d = distanceToNearest(HepG2[17],GM12878)
i = subjectHits(d)
start(GM12878[i])
"46524762"

QUESTION 2.5.4
d = distanceToNearest(HepG2[17],GM12878)
mcols(d)$distance
"2284"

QUESTION 2.5.5
d = distanceToNearest(HepG2,GM12878)
mean( mcols(d)$distance < 2000)
"0.2673267"

QUESTION 2.6.1
genome(ghs)
"hg19"
length(ghs)
"23056"

QUESTION 2.6.2
which.max( table( seqnames( ghs ) ))
"chr1"

QUESTION 2.6.3
w = width( ghs )
hist( w )
## the larger values are so much larger than the bulk of the data
## that we can barely see the frequencies of large values in the histogram 
## the log transformation better shows how fat the right tail is:
hist( log10(w))
"A distribution with a very fat right tail A distribution with a very fat right tail"

QUESTION 2.6.4
w = width( ghs )
median(w)
"20115.5"

QUESTION 2.7.1  
## first order them
erbs3 = erbs[order(erbs),]
##confirm same chr
all( seqnames(erbs2)==seqnames(erbs3) )
mean( start(erbs2)==start(erbs3) & end(erbs2)==end(erbs3) )
##the intersection should be smaller
all( width(erbs2) <= width(erbs3) )
"Over 90% of these regions in these two objects are the same with the different regions being smaller in erbs2."

QUESTION 2.7.2
tssgr= resize(ghs,1)
start(tssgr["100113402"])
"70563402"

QUESTION 2.7.3
library(Homo.sapiens)
ghs = genes(Homo.sapiens)
tssgr= resize(ghs,1)
i = nearest(erbs[4],tssgr)
mcols(tssgr)$GENEID[i]
"2101"

QUESTION 2.7.4
library(Homo.sapiens)
ghs = genes(Homo.sapiens)
tssgr= resize(ghs,1)
i = nearest(erbs[4],tssgr)
gene = as.character(mcols(tssgr)$GENEID[i])

select(Homo.sapiens,key=gene,column="SYMBOL",keytype="GENEID")
"ESRRA"

QUESTION 2.8.1
genome(erbs)
"hg19"

QUESTION 2.8.2
seqs = getSeq(Hsapiens,erbs)
gc = alphabetFrequency(seqs)[,2:3]
n = width(erbs)
gccontent = rowSums(gc)/n
"0.652568"

QUESTION 2.8.3
control = shift(erbs,10000)
controlseqs = getSeq(Hsapiens,control)
gc = alphabetFrequency(controlseqs)[,2:3]
n = width(control)
controlgccontent = rowSums(gc)/n
median(controlgccontent)
"0.4860174"

QUESTION 2.9.1
library(BSgenome)
grep("Drerio", available.genomes(), value=TRUE) # exclude masked
"3"

QUESTION 2.9.2
class(c17m)
"MaskedDNAString"

QUESTION 2.9.3
c22m = BSgenome.Hsapiens.UCSC.hg19.masked$chr22
  round(100*sum(width(masks(c22m)$AGAPS))/length(c22m),0)
  "31"
  
QUESTION 2.10.1
start(HepG2[1])-start(nHepG2[1])
"199761"

QUESTION 2.11.1
The linear structures are transcripts, composed of untranslated regions at the ends, exons (small yellow rectangles), and introns (arrows between exons). The collection of transcripts makes up a "gene model".
"TRANSCRIPT"

QUESTION 2.11.2
"count them in the display" "27"

QUESTION 2.11.3
"27" 
 
Q 2.12.1
class(targets)
"data.frame"

Q 2.12.2
 "it cannot be determined without forensic work it cannot be determined without forensic work"

Q 2.12.3
 "use a descriptive filename use a descriptive filename"

Q 2.13.1
"8"

Q 2.14.1
resize(g,1)
"kp is the set of genes whose transcription start site lies in a HepG2 binding site"

Q 2.14.2
You enter the string in the search box and look at the bottom of the page which tells how many entries are found
"30"


##WEEK 3


Q 3.1.1
The output of a microarray experiment is light intensity, which is recorded as rational valued numbers.
"More copies of a gene in a sample results in higher light intensity at the corresponding probe"

Q 3.1.2
Microarrays do not directly measure gene expression or number of RNA molecules, but the light intensity of labeled cDNA hybridizing to probes printed on the array. For this reason, probe-specific effects make it difficult to compare gene expression/RNA abundance across probes for a given sample.
"the hybridization of RNA or cDNA to probes the hybridization of RNA or cDNA to probes"

Q 3.1.3
Golub, et al. 1999 (PMID: 10521349) in particular used gene expression microarrays to cluster tumors into classes for acute myeloid leukemia (AML) and acute lymphoblastic leukemia (ALL). Rare variants underlying monogenic disorder and transcription factor binding sites cannot be directly inferred by the output of gene expression microarrays.
 "classifying tumor samples from various patients into distinct phenotypic classes"
 
 Q 3.2.1
 The main advance of NGS technology is that we can process many more sequences in parallel. However, to achieve this we are limited (as of 2015) in the size of these sequences.
 " NGS produces larger volumes of shorter sequenceS"
 
 Q 3.2.2
 When the signal intensity for a given cluster is about equal for each nucleotide, the base calling algorithms cannot say for certain which nucleotide was incorporated.
 "the signal intensity for each nucleotide (A,C,G or T) was nearly equal "
 
 Q 3.2.3
 After aligning the sequencing reads to the genome, we can look for columns of reads where there are nucleotides different than the reference genome. However, for microarray genotyping, the locations of SNPs must already be known, in order to be printed on the array.
 "we can find new SNPs"
 
 Q 3.2.4
 All of the first three choices can and do occur in genomics. For examples of each: genomes of cancerous cells involve translocations, which can create regions of the genome which do not occur in the reference genome. There are regions which occur 1000s of times in a genome, called transposons. Many software pipelines discard reads which map to many locations, and then true biology (e.g. a bound protein) which can occur at those regions is ignored. Also, in sequencing RNA or DNA from cancer cells, there can be so many mutations that the aligning programs have difficulty finding the match in the reference genome. Or, this can happen if we align reads from one species to the genome of another species (which does not yet have a reference genome constructed).
"The region we care about is not in the reference genome"
"The region we care about occurs 1000s of times in the genome, and we are ignoring reads which align to so many places"
"For the region we care about, the organism's genome is so different to the reference genome that the aligning program can't find a match"
 
 Q 3.2.5
 We can generally expect linear scaling with the number of reads which align to a genomic location for sequencing experiments. However, note that the question insisted we sequencing more DNA from the same DNA "library", which is called a "technical replicate". If we perform a new experiment with a new organism/tissue/population of cells, which is called a "biological replicate", we know there is biological variation in the underlying quantity: all organisms/tissues/cells do not have equal level of mRNA.
"around 3000" 
 
 
