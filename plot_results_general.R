#!/usr/bin/env Rscript
#### Plot spike alignment results
#### 2014-05-15 Inez Terpstra
#----------------------------------------------------------------------------------------------------
print("Start visualization")
#----------------------------------------------------------------------------------------------------
args <- commandArgs(TRUE)

base.dir <- args[1]
experiment <- args[2]

#----------------------------------------------------------
suppressPackageStartupMessages(library(ShortRead))
#----------------------------------------------------------------------------------------------------------
dir.fq <- paste(base.dir,"/Scratch/Fastq",sep="")
dir.res <- paste(base.dir,"/Results/",experiment,"/Output",sep="") 
dir.plot <- paste(base.dir,"/Images/",experiment,sep="") 

#-------------------------------------------------------------------------------------
#base quality boxplot
#-----------------------------------------------------------------------------------------
setwd(dir.fq)
fq.files <- dir(pattern="fastq")
fq.file <- sample(fq.files,1)
chip <- sub(".+-","",fq.file)
chip <- sub("\\.fastq","",chip)
barcode <- strsplit(fq.file,"_")[[1]][2]

fq <- FastqSampler(fq.file,5e5) 
fq.sample <- yield (fq)

setwd(dir.plot)

titel <- paste("quality scores sample",chip,barcode,sep=" ")
fq.qa <- as(quality(fq.sample), "matrix")
png(filename = paste(titel,".png",sep=""), width=635, height=460, units="px", bg="white")
  boxplot(as.data.frame((fq.qa)), outline=FALSE, names=c(1:dim(fq.qa)[2]), xlab="position in read bp", main=titel)
dev.off()


