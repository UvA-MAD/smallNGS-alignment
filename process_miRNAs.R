#!/usr/bin/env Rscript
#### Processing miRNA aligmnents
#### 2014-05-12 Inez Terpstra
#### Alignments are filtered in R
#### - 10% mismatches are allowed
#### - each InDel counts as 1 mismatch
#### - a 3p overhang of 3 nucleotides is allowed
#### - aligned reads are assigned to mature miRNA if the read extends the mature sequence with 3 or less nucleotides
#### unmapped reads are written to new fastq file
#### total reads, reads lengths and reads counts (in getMappingData) are saved
#----------------------------------------------------------------------------------------------------
print("Processing miRNA alignments")
#----------------------------------------------------------
# Analyze NGS miRNA data per barcode
#----------------------------------------------------------
args <- commandArgs(TRUE)

base.dir <- args[1]
experiment <- args[2]
filename <- args[3]
chip <- args[4]
barcode <- args[5]
organism <- args[6]
mapped.file <- args[7]
dir.db <- args[8]
db <- args[9]
fgp <- args[10]

#----------------------------------------------------------
suppressPackageStartupMessages(library(ShortRead))
#----------------------------------------------------------------------------------------------------------

dir.fq <- paste(base.dir,"/Scratch/Fastq",sep="")
dir.map <- paste(base.dir,"/Scratch/",experiment,"/Data/miRNA_aligned_data",sep="")
dir.umap <- paste(base.dir,"/Scratch/",experiment,"/Fastq/unmapped_miRNA",sep="")
dir.res <- paste(base.dir,"/Results/",experiment,"/Output",sep="") 
dir.plot <- paste(base.dir,"/Images",experiment,sep="")
dir.deseq <- paste(base.dir,"/Results/",experiment,"/DESeq",sep="")

dir.src <- "/mad/MAD-RBAB/10_Protocols/DryLab/smallNGS_alignment" 

#==========================================================================================
# functions
#---------------------------------------------------------------------------------------------------------
source(paste(dir.src,"spike_and_miRNA_functions.R",sep="/"))

#-------------------------------------------------------------------------------------------------------------------------
#retrieve miRNA data
#-------------------------------------------------------------------------------------------------------------------------
setwd(dir.db)
load(file="miRNAdata.RData")

#----------------------------------------------------------------------------------------------------------
# read in alignments of truncated reads to miRNAs mapped with SHRiMP2 and 
# filter alignments 
# read in mapping results and raw unmapped data (to remove mapped reads for further alignments)
#----------------------------------------------------------------------------------------------------------
run.name <- paste(chip,barcode,sep="_")

total.reads <- matrix(0,1,19)
colnames(total.reads) <- c("unmapped", "nt20", "nt21", "nt22", "nt23", "nt24", "nt25", "nt26", "nt27", "nt28", "nt29", "nt30", "nt31", "nt32", "nt33", "nt34", "nt35", "hairpin", "mature")
rownames(total.reads) <- run.name

bam.file.global <- paste("sorted",mapped.file,"global.bam",sep="_")
bam.file.local <- paste("sorted",mapped.file,"local.bam",sep="_")

mapped.mature <- matrix(0,length(unique(miRNA$mature)),1)
rownames(mapped.mature) <- unique(miRNA$mature)
colnames(mapped.mature) <- run.name

mapped.hairpin <- matrix(0,length(hairpin),1)
rownames(mapped.hairpin) <- pre.miRNA$ID
colnames(mapped.hairpin) <- run.name

read.lengths <- NULL
#list read lengths of mapped reads before trimming
aln.miRNA.read.lengths <- NULL

#----------------------------------------------------------------------------------------------------------
## parameters
#filter on valid alignments
#percent mismatches allowed
percent.mismatches <- 10
#number of unaligned nucleotides allowed at the ends
unmatched5p <- 0  
unmatched3p <- 3
# read is allowed to extend the mature sequence with x nucleotides
allowMatureExtension <- 3

#----------------------------------------------------------------------------------------------------------
## read in alignments
setwd(dir.map)
bam.mapped <- getMappingData(run.name, bam.file.local, bam.file.global, percent.mismatches, unmatched5p, unmatched3p, allowMatureExtension, pre.miRNA, miRNA, fgp)

print("bam read")

load(file=paste("obj",run.name,"mapped.mature.out",sep="_"))
load(file=paste("obj",run.name,"mapped.hairpin.out",sep="_"))

total.reads[18] <- colSums(mapped.hairpin)
total.reads[19] <- colSums(mapped.mature)

fq.new <- paste(dir.umap,"/",run.name,"_unm.fastq",sep="")
if (file.exists(fq.new)) file.remove(fq.new)

#read in fastq file and write unmapped reads to new fastq file
strm <- FastqStreamer(filename)

chunk <- 0
repeat {
  fq <- yield(strm)
  chunk <- chunk+1
  if (length(fq) == 0)
    break
  ## process chunk
  read.lengths <- c(read.lengths,width(sread(fq)))
  total.reads[1] <- total.reads[1]+length(fq)
  total.reads[2] <- total.reads[2]+sum(read.lengths==20)
  total.reads[3] <- total.reads[3]+sum(read.lengths==21)
  total.reads[4] <- total.reads[4]+sum(read.lengths==22)
  total.reads[5] <- total.reads[5]+sum(read.lengths==23)
  total.reads[6] <- total.reads[6]+sum(read.lengths==24)
  total.reads[7] <- total.reads[7]+sum(read.lengths==25)
  total.reads[8] <- total.reads[8]+sum(read.lengths==26)
  total.reads[9] <- total.reads[9]+sum(read.lengths==27)
  total.reads[10] <- total.reads[10]+sum(read.lengths==28)
  total.reads[11] <- total.reads[11]+sum(read.lengths==29)
  total.reads[12] <- total.reads[12]+sum(read.lengths==30)
  total.reads[13] <- total.reads[13]+sum(read.lengths==31)
  total.reads[14] <- total.reads[14]+sum(read.lengths==32)
  total.reads[15] <- total.reads[15]+sum(read.lengths==33)
  total.reads[16] <- total.reads[16]+sum(read.lengths==34)
  total.reads[17] <- total.reads[17]+sum(read.lengths==35)
  rm.idx <- id(fq) %in% bam.mapped$qname
  aln.miRNA.read.lengths <- c(aln.miRNA.read.lengths,width(sread(fq[rm.idx])))
  fq.unmapped <- fq[rm.idx==FALSE]
  writeFastq(fq.unmapped,fq.new,full=TRUE,"a")
  print(paste("unmapped fastq chunk",chunk,"written",sep=" "))
}
close(strm)

save(read.lengths,file=paste("obj",run.name,"read.lengths.out",sep="_"))
save(aln.miRNA.read.lengths,file=paste("obj",run.name,"aln.miRNA.read.lengths.out",sep="_"))
save(total.reads,file=paste("obj",run.name,"total.reads.out",sep="_"))


