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
"In an eSet object with the e in the assayData accessible with exprs() and tissue as one of the columns in the phenoData"
 
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
 
 Q 3.3.1
 "ExpressionSet"
 
 Q 3.3.2
"geneExpression will be the assayData, sampleInfo will be the phenoData, and expressionAnnotation will be the featureData" 
 
 Q 3.3.3
 pd = AnnotatedDataFrame(sampleInfo)
rownames(pd) = colnames(geneExpression)
pd = pd[, colnames(pd)!="filename"] ##redundant
The reason we use the class AnnotatedDataFrame, as opposed to just using data frames, is to encourage users to describe the variables represented in the table. Learn more here:

?AnnotatedDataFrame
"2005-06-27"

 Q 3.3.4
varLabels(pd)[1]
"ethnicity"

 Q 3.3.5
fd = AnnotatedDataFrame(geneAnnotation)
rownames(fd) = geneAnnotation$PROBEID
fd = fd[,colnames(fd)!="PROBEID"] ##redundant
pData(fd)["204810_s_at","CHR"]
[explanation]
"chr19"
 
 Q 3.3.6
 pd = AnnotatedDataFrame(sampleInfo)
rownames(pd) = pd$filename

fd = AnnotatedDataFrame(geneAnnotation)
rownames(fd) = geneAnnotation$PROBEID

eset = ExpressionSet(assayData=geneExpression,
              phenoData=pd,
              featureData=fd)
Now that you have the data organized in such a way that you can access the different components from the same object:

dim (pData(eset))
dim( exprs(eset) )
dim( featureData( eset ))
We also know (because we created the original datasets) that these measurements were made with the hgfocus array, thus we can add this as well:

annotation(eset) = "hgfocus"
"0.8049265"
 
 Q 3.3.7
  eset[1:10,1:5] 
 
 
 Q 3.3.8
tss = start(resize( granges(se),1))
sum(  tss < 50*10^6 & seqnames( se)=="chr1" )

### we will re-order se
se = se[order(granges(se)),]
ind = se$group==1
de = rowMeans( assay(se)[,ind])-rowMeans( assay(se)[,!ind])
chrs = unique( seqnames(se))
library(rafalib)
mypar2(3,2)
for(i in c(1:4)){
  ind = which(seqnames( se) == chrs[i])
  plot(start(se)[ind], de[ind], ylim=c(-1,1),main=as.character(chrs[i]))
  abline(h=0)
  }
##now X and Y
for(i in 23:24){
  ind = which(seqnames( se) == chrs[i])
  ##note we use different ylims
  plot(start(se)[ind], de[ind], ylim=c(-5,5),main=as.character(chrs[i]))
  abline(h=0)
  }
  "265"

 Q 3.4.1
 datadir = "/your/path/here"
basedir = paste0(datadir, "/celfiles")
setwd(basedir)
tab = read.delim("sampleinfo.txt",check.names=FALSE,as.is=TRUE)
rownames(tab) = tab$filenames
tab["1521a99hpp_av06.CEL.gz","36311_at"]
"4"
 
Q 3.4.2
fns = list.celfiles()
fns
all(fns %in% tab[,1])##check
ab = ReadAffy(phenoData=tab)
sum( probeNames(ab)=="36311_at")
length( featureNames(ab) )
length( probeNames(ab))
"16"

Q 3.4.3
pid = "36085_at"
##which columns should we use?
ind = which(pData(ab)[,1]%in%c("1532a99hpp_av04.CEL.gz","1532b99hpp_av04.CEL.gz"))

##extract the correct rows
mat = pm(ab) [ probeNames(ab)==pid, ind] 

##what are the intended conc
conc = pData(ab)[ind, pid]

##make the plots
mypar2(1,1)
matplot(conc, t(mat), log="y", type="l")

##now comput log fold changesa
lfc = log2(mat[,2]/mat[,1])
stripchart(lfc,vertical=TRUE,ylim=c(-0.5,1.5))
abline(h=log2(conc[2]/conc[1])) #intended log fold
abline(h=0)
"There is wide variation across probes representing the same gene but these mostly cancel out and the log-ratios are close to the intended log fold change of 1 (between 0.3 and 1.2)"

Q 3.4.4
library(genefilter)
g = factor(pData(ab)[,2])
e = rma(ab)
tt = rowttests(exprs(e),g)
tt["36085_at","p.value"]
setwd(basedir)
ejust = justRMA(filenames=tab[,1],phenoData=tab)
"0.002611949"

Q 3.4.5
g = factor(pData(e)[,2])
tt = rowttests(exprs(e),g)
lfc = -tt$dm

sig = colnames(pData(ab))[-1]
boxplot( split(lfc, rownames(tt)%in%sig))
##close up 
boxplot( split(lfc, rownames(tt)%in%sig),ylim=c(-1,1))
" With one exception, the spiked-in genes show absolute value for log fold changes of about 0.6, while the rest of the genes are mostly between -0.25 and 0.25 "

Q 3.4.6
i = which(RG$genes$ID=="H200015482")
j = which(rownames(RG$targets)=="6Hs.168")
log2(RG$R[i,j]/RG$G[i,j])
MA = MA.RG(RG,bc.method="none")
MA$M[i,j]
setwd(wd)
"0.6456647"

QUESTION 3.5A.1  
bf = BamFile(filename)

gr = GRanges("chr4",IRanges(440000, 470000))

countBam(bf, param=ScanBamParam(which=gr))

# or, alternatively:

reads = scanBam(bf, param=ScanBamParam(what="pos", which=gr))

length(reads[[1]]$pos)
"1665"

QUESTION 3.5A.2
bf = BamFile(filename)

gr = GRanges("chr4",IRanges(440000, 470000))

reads = scanBam(bf, param=ScanBamParam(what="seq", which=gr))

mean(letterFrequency(reads[[1]]$seq, "GC", as.prob=TRUE))
"0.4452773"

QUESTION 3.5A.3
countOverlaps(g2[ "FBgn0039890" ], ga)
"1085"

QUESTION 3.5B.1
g.so = summarizeOverlaps(g, bf, ignore.strand=TRUE)

grl.so = summarizeOverlaps(grl, bf, ignore.strand=TRUE)

# this generates a warning about 0 counts

plot(assay(g.so),assay(grl.so),log="xy");abline(0,1)

# adding a pseudocount fixes the log(0) issue

plot(assay(g.so)+1,assay(grl.so)+1,log="xy");abline(0,1)

ratio = assay(grl.so) / assay(g.so)

mean(ratio[assay(g.so) > 0])
"0.8929273"

QUESTION 3.5B.2
count = assay(grl.so)

# here, we can just use sum() because there is only one column

# otherwise we would use: sweep(count, 2, colSums(count), "/")

fpm = (count/sum(count)) * 1e6

head(fpm,1)
"4275.60731"

QUESTION 3.5B.3
ebp = sum(width(reduce(grl)))

count = assay(grl.so)

fpm = (count/sum(count)) * 1e6

fpkm = (fpm/ebp) * 1e3

head(fpkm)
"1782.245648"

QUESTION 3.6.1  
p0s = colMeans( exprs(bottomly.eset) == 0)
median( p0s )
"0.687021"

QUESTION 3.6.2
p0s = colMeans( exprs(bottomly.eset) == 0)
boxplot(split( p0s , pData(bottomly.eset)$experiment.number))
"The proportion of 0s varies by about 2% across the batches and there appears to be statistically significant differences, with experiment 7 having the lowest values." 

QUESTION 3.6.3
d = dist( t(y) )
mds = cmdscale(d)
batch = pData( bottomly.eset)$experiment.number - 3
strain = as.numeric(pData (bottomly.eset)$strain)
library(rafalib)
mypar2(1,1)
plot(mds,col=batch,pch=strain)
legend("topleft",col=unique(batch),legend=unique(batch)+3,pch=1)
legend("bottomleft",pch=unique(strain),legend=unique(strain))
"Batch seems to explain more variability that strain. The first PC splits batch 7 from 3 and 4 into two groups. The second PC splits batches 3 and 4 into two groups (with one sample as an exception)" 

QUESTIONS 3.7
EXPLANATION
-The pair of u and v is more likely to have a high correlation than differences or log fold changes across pairings. This is because the differences and log fold changes involve moving the values of the pair closer to zero in which the individual samples had large values. High correlations in microarray data are driven by rows where both of the pair have large values.
-The raw values reported from a microarray are not directly translatable to due to background noise, sample- and probe-specific effects. In particular, a value of 0 does not correspond to “no expression”, and a value above 0 does not correspond to "expressed".
-The raw values reported from a microarray are not directly translatable to inference about expression, therefore direct relative comparisons of these raw values across different sets of proves are not meaningful.
"Poisson with lambda = 1/2"   "u and v"  "Because different probes have different background levels, we cannot say for sure if gene A is expressed or not"  "Because not all probes are as good at detecting real gene expression, we cannot say for sure if gene B has higher expression than gene A" 

QUESTION 3.7.1
ind=which(pns==colnames(pd)[1]) ##probes in gene 1
concentration=pd[,1]
i = which(concentration==0)
max( pms[ind,i] )
"26"   "468.3"

QUESTION 3.7.2
ind = which(pns==colnames(pd)[j]) ##probes in gene 1
y = as.vector(pms[ind,i])
x = as.vector(mms[ind,i])
plot(x,y)
plot(log2(x),log2(y))
cor(log2(x),log2(y))
hist(log2(x)-log2(y))
"0.9305712"

QUESTION 3.7.3
ind = c(1,15,29)
pm1 = log2( pm(bg1)[,ind])
pm2 = log2( pm(bg2)[,ind])

SD1 = rowSds(pm1)
A1 = rowMeans(pm1)
SD2 = rowSds(pm2)
A2 = rowMeans(pm2)
mypar2(2,1)
splot(A1,SD1,ylim=c(0,3),cex=.25)
splot(A2,SD2,ylim=c(0,3),cex=.25)
 "Method 1 shows larger variability especially for the lower values of A."
 
 QUESTION 3.8.1
 EXPLANATION
The plot depicts looking up the quantile, or propotion of genes with a smaller value, of a gene on a single array (top left to bottom left panel), and matching this quantile with the quantile of a multi-array average of sorted gene expression vectors (bottom left to bottom right panel), then outputting the corresponding value of the multi-array average (bottom right to top right panel). This is the quantile normalization method.
Quantile normalization produces a set of vectors which have the exact same empirical distribution function.
If distributions of the gene expression across samples are known to be different, e.g. lots of highly expressed genes in one group and only few in another group, then the assumptions of quantile normalization – that the differences in values of the i-th highest gene of each array are due to technical artifacts – are not met.
Bigger window sizes corresponds to more smooth curves from averaging the lines from many observations. Very small window sizes leads to the curve fitting around each point, leading to a spikier, jagged curve.
"Quantile normalization"  "The distributions"  "You expect the distributions to be different"  "Smoother curve"

QUESTION 3.8.2  
median(sds[spikeinIndex])
boxplot(sds[-spikeinIndex],sds[spikeinIndex],range=0)
"0.1990065"  "1.820172"

QUESTION 3.8.3  
library(preprocessCore)
npm = normalize.quantiles(pm)
nsds=rowSds(log2(npm))
median(nsds[-spikeinIndex])
boxplot(sds[-spikeinIndex],nsds[-spikeinIndex],range=0,ylim=c(0,0.5))
"0.1274938"

QUESTION 3.8.4  
library(preprocessCore)
npm = normalize.quantiles(pm)
nsds=rowSds(log2(npm))
median(nsds[spikeinIndex])
boxplot(sds[-spikeinIndex],nsds[-spikeinIndex],sds[spikeinIndex],nsds[spikeinIndex],range=0)
"1.813862"

QUESTION 3.8.5
mypar2(2,1)
boxplot(M[erccIndex,],range=0,ylim=c(-2,2))
abline(h=0,lty=2)
boxplot(M[-erccIndex,],range=0,ylim=c(-2,2))
abline(h=0,lty=2)
 "Quantile normalization has centered the majority of genes. The controls are only a minority so this is the preferred result"
 
QUESTION 3.9.1
x = getSeq(Mmusculus, reduced.exons)

dss = DNAStringSet(lapply(x, unlist))

gc = letterFrequency(dss, "GC", as.prob=TRUE)
"0.4151355"

QUESTION 3.9.2
Examine the plot.
"7"  "4"

QUESTION 3.10.1
summary(width(h2bw))
4

QUESTION 3.10.2
sum(duplicated(inpeak))
## no duplicated indices here, so:
median(h2bw[ inpeak ]$score)
"45.7"

QUESTION 3.10.3
median(h2bw[ -inpeak ]$score)
"4.2"

QUESTION 3.10.4
library(Homo.sapiens)

select(Homo.sapiens, keys="ESRRA", keytype="SYMBOL", columns="CHRLOC")

narrind = queryHits(findOverlaps(HepG2, GRanges("chr11", IRanges(64073044, width=1))))

bwind = queryHits(fo)[ subjectHits(fo)==narrind]

max( h2bw$score[ bwind ] )
"202.1"

QUESTION 3.10.5
peakcov[ which.max(peakcov$score) ]
start(HepG2[5]) + HepG2[5]$peak
"64072320"
 



