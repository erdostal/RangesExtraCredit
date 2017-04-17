# Genomic Ranges Extra Credit
## Emma Miller

### Set up R environment in order to complete assignment
```{r}
#Load required packages for this pipeline
library(IRanges)
library(GenomicRanges)
library("GenomicFeatures")
library(rtracklayer)
```

### Find and import proper files
```{r}
#Import Mouse annotation library and save it as "txdb"
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene

#Import mouse data as a GRanges object
mm_gtf <- import('Mus_musculus.GRCm38.75_chr1.gtf.gz')
```
### Draw out only exons involved in coded proteins and double check that dataframe
```{r}
#Show number of genes per category of exon
table(mm_gtf$gene_biotype)
#create new object with only protein coding genes
chr1_pcg <- mm_gtf[mm_gtf$type == "gene" &
                     mm_gtf$gene_biotype == "protein_coding"]
#Making sure that this data makes since given that this is only from chromosome1
summary(width(chr1_pcg))
length(chr1_pcg)
```
### Take exons involved in protein coding regions, and generate promoter regions with two strategies
```{r}
#Creates a range width of 3000 upstream of the gene range for "promoter" region
chr1_pcg_3kb_up <- flank(chr1_pcg, width = 3000)
#Same idea using the promoters command instead of flank
chr1_pcg_3kb_up2 <- promoters(chr1_pcg, upstream = 3000, downstream=0)
#Both functions give you the same flanking region 
identical(chr1_pcg_3kb_up, chr1_pcg_3kb_up2)
```

### Install and load packages and sequence data, and double check some features of that data
```{r}
#Install Biostrings packages 
library(BiocInstaller)
biocLite("BSgenome")
bioClite("BSgenome.Mmusculus.UCSC.mm10")
#load genome and give new shortcut name
library(BSgenome.Mmusculus.UCSC.mm10)
mm_gm <- BSgenome.Mmusculus.UCSC.mm10
#Check file
organism(mm_gm)
providerVersion(mm_gm)
provider(mm_gm)
#Look at sequence information, stored like a list so can be accessed with indeces
seqinfo(mm_gm)
mm_gm$chrM
mm_gm[[22]]
```

### Not needed for this particular assignment, but Biostrings can be used to search for a short particular string
```{r}
#Search sequence for a particular string
library(Biostrings)
matchPattern("GGcGCGCC", mm_gm$chr1)
```

### Or it can be used to adjust  labels within objects in order to match certain databases
### Here we are making sure chromosome labels match between our range data and NCBI sequence data
```{r}
#Adjust sequence level names from Chr1 and Chr2 to "1" and "2" to match NCBI
#Checking that sequences are BSgenome object
all(seqlevels(chr1_pcg_3kb_up) %in% seqlevels(mm_gm))
#Creating GRanges object to change chromosome names
gr <- GRanges(c("chr1", "chr2"), IRanges(start=c(3,4),width= 10))
seqlevels(gr)
#Change sequence levels
seqlevels(gr) <- c("1", "2")
seqlevels(gr)
seqlevelsStyle(chr1_pcg_3kb_up)
seqlevelsStyle(mm_gm)
#function provied to change style
seqlevelsStyle(chr1_pcg_3kb_up) <- "UCSC"
#Demonstrate that items are consistent
all(seqlevels(chr1_pcg_3kb_up) %in% seqlevels(mm_gm))
```

### Get the sequences for the promoter regions we are intereseted in
```{r}
#Takes BSgenome and the GRanges and give back the sequences we are looking for
chr1_pcg_3kb_seqs <- getSeq(mm_gm, chr1_pcg_3kb_up)
chr1_pcg_3kb_seqs
```

### Write these sequences to a new fasta file for downstream application
```{r}
writeXStringSet(chr1_pcg_3kb_seqs,file='mm10_chr1_3kb_promoters.fasta')
```
