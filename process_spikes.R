#!/usr/bin/env Rscript
#### Processing spike aligmnents
#### 2014-05-12 Inez Terpstra
#### Alignments are filtered in R
#### - aligned reads should be 80% of spike length
#### unmapped reads are written to new fastq file
#### total reads, reads lengths and reads counts are saved
#----------------------------------------------------------------------------------------------------
print("Processing spike alignments")
#----------------------------------------------------------
# Analyze NGS spike data per barcode
#----------------------------------------------------------
args <- commandArgs(TRUE)

base.dir <- args[1]
experiment <- args[2]
filename <- args[3]
chip <- args[4]
barcode <- args[5]
mapped.file <- args[6]
dir.db <- args[7]
db <- args[8]

#----------------------------------------------------------
suppressPackageStartupMessages(library(ShortRead))
#----------------------------------------------------------------------------------------------------------

dir.fq <- paste(base.dir,"/Scratch/Fastq",sep="")
dir.map <- paste(base.dir,"/Scratch/",experiment,"/Data/spike_aligned_data",sep="")
dir.umap <- paste(base.dir,"/Scratch/",experiment,"/Fastq/unmapped_spikes",sep="")
dir.res <- paste(base.dir,"/Results/",experiment,"/Output",sep="") 
dir.plot <- paste(base.dir,"/Images",experiment,sep="")
dir.deseq <- paste(base.dir,"/Results/",experiment,"/DESeq",sep="")

dir.src <- "/mad/MAD-RBAB/10_Protocols/DryLab/smallNGS_alignment" 

#==========================================================================================
# functions
#---------------------------------------------------------------------------------------------------------
source(paste(dir.src,"spike_and_miRNA_functions.R",sep="/"))

#----------------------------------------------------------------------------------------------------------------------------------------
# retrieve spike data
#-------------------------------------------------------------------------------------------------------------------------
spike.data <- readDNAStringSet(paste(dir.db,"spikes.fasta",sep="/"))
spikes <- names(spike.data)
NS.spikes <- spikes[grep("NSPK",spikes)]
SS.spikes <- spikes[grep("SSPK",spikes)]

#----------------------------------------------------------------------------------------------------------
# read in alignments of truncated reads to spikes mapped with Bowtie2 and 
# filter alignments on read length compared to spike length
# read in mapping results and raw data (to remove mapped reads for further alignments)
#----------------------------------------------------------------------------------------------------------
run.name <- paste(chip,barcode,sep="_")

total.reads <- matrix(0,1,20)
colnames(total.reads) <- c("total", "nt20", "nt21", "nt22", "nt23", "nt24", "nt25", "nt26", "nt27", "nt28", "nt29", "nt30", "nt31", "nt32", "nt33", "nt34", "nt35", "size", "norm", "spike")
rownames(total.reads) <- run.name

bam.file <- paste("sorted_",mapped.file,".bam",sep="")

mapped.spike <- matrix(0,length(spikes),1)
rownames(mapped.spike) <- spikes
colnames(mapped.spike) <- run.name

read.lengths <- NULL
#list read lengths of mapped reads before trimming
aln.spike.read.lengths <- NULL

## read in alignments
setwd(dir.map)
bam <- getBam(bam.file)
  
#select valid alignments only:
#read length is 80% of spike length for small spikes and 32 for long spikes (>= 40)
bam$filter <- NULL
long.spikes <- spikes[width(spike.data)>=40]
long.idx <- bam$rname %in% long.spikes
bam$filter[long.idx] <- bam$qwidth[long.idx] >= 32
small.spikes <- spikes[width(spike.data)<40]
small.idx <- bam$rname %in% setdiff(small.spikes,NS.spikes)
bam$filter[small.idx] <- bam$qwidth[small.idx] >= 0.8*as.double(substring(bam$rname[small.idx],10))
#normalisation spikes are 25 nucleotides long, but length cannto be retrieved from spike name
norm.idx <- bam$rname %in% NS.spikes
bam$filter[norm.idx] <- bam$qwidth[norm.idx] >= 0.8*25
  
#mapped.reads when bam$filter==TRUE
bam.nm <- bam[bam$filter,]

total.reads[18] <- sum(as.vector(bam.nm$rname) %in% SS.spikes)
total.reads[19] <- sum(as.vector(bam.nm$rname) %in% NS.spikes)
total.reads[20] <- dim(bam.nm)[1]

fq.new <- paste(dir.umap,"/",run.name,"_unm.fastq",sep="")
if (file.exists(fq.new)) file.remove(fq.new)

#read in fastq file and write unmapped reads to new fastq file
strm <- FastqStreamer(paste(dir.fq,"/",filename,sep=""))

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
  rm.idx <- id(fq) %in% bam.nm$qname
  aln.spike.read.lengths <- c(aln.spike.read.lengths,width(sread(fq[rm.idx])))
  fq.unmapped <- fq[rm.idx==FALSE]
  writeFastq(fq.unmapped,fq.new,full=TRUE,"a")
  print(paste("unmapped fastq chunk",chunk,"written",sep=" "))
}
close(strm)

mapped.spike[names(table(bam.nm$rname)),1] <- table(bam.nm$rname)

NM <- table(bam.nm$NM)
mapped.spike.NM <- matrix(0,length(spikes),length(NM))
colnames(mapped.spike.NM) <- paste("NM",names(NM),sep="")
rownames(mapped.spike.NM) <- spikes
for (j in names(NM)){
  mapped.spike.NM[names(table(bam.nm$rname[bam.nm$NM==j])),paste("NM",j,sep="")] <- table(bam.nm$rname[bam.nm$NM==j])
}

save(mapped.spike,file=paste("obj",run.name,"mapped.spike.out",sep="_"))
save(mapped.spike.NM,file=paste("obj",run.name,"mapped.spike.NM.out",sep="_"))
save(read.lengths,file=paste("obj",run.name,"read.lengths.out",sep="_"))
save(aln.spike.read.lengths,file=paste("obj",run.name,"aln.spike.read.lengths.out",sep="_"))
save(total.reads,file=paste("obj",run.name,"total.reads.out",sep="_"))


