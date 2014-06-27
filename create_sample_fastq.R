base.dir <- "/mad/MAD-RBAB/02_Collaborators/MAD1022-AMC_Peter/MAD1022-P007-Sequencing/MAD1022-P007-E001_2014_RNA_Seq_svleeuw1"
experiment <- "rRNA_test"

#----------------------------------------------------------
suppressPackageStartupMessages(library(ShortRead))
#----------------------------------------------------------------------------------------------------------
dir.fq <- paste(base.dir,"/Scratch/Fastq",sep="")
dir.res <- paste(base.dir,"/Results/",experiment,"/Output",sep="") 
dir.plot <- paste(base.dir,"/Images/",experiment,sep="") 

#-------------------------------------------------------------------------------------
#take sample of fastq file and write to file
#-----------------------------------------------------------------------------------------
setwd(dir.fq)
fq.files <- dir(pattern="fastq")
fq.file <- sample(fq.files,1)
chip <- sub(".+-","",fq.file)
chip <- sub("\\.fastq","",chip)
barcode <- strsplit(fq.file,"_")[[1]][2]
fq <- FastqSampler(fq.file,1000000) 
fq.sample <- yield(fq)
new_fq <- sub(chip,paste("sub",chip,sep=""),fq.files)
writeFastq(fq.sample,new_fq)
