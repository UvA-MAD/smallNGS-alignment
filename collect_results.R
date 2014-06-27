#!/usr/bin/env Rscript
#### Collect results from aligmnents of different barcodes
#### 2014-05-15 Inez Terpstra
#### total reads and reads counts are combined per experiment
#----------------------------------------------------------
args <- commandArgs(TRUE)

base.dir <- args[1]
experiment <- args[2]
small <- args[3]

#----------------------------------------------------------------------------------------------------
print(paste("Collecting",small,"count data",sep=" "))

#==========================================================================================
# environment
#---------------------------------------------------------------------------------------------------------
dir.map <- paste(base.dir,"/Scratch/",experiment,"/Data/",small,"_aligned_data",sep="") 
dir.res <- paste(base.dir,"/Results/",experiment,"/Output",sep="") 
dir.deseq <- paste(base.dir,"/Results/",experiment,"/DESeq",sep="")

#---------------------------------------------------------------------------------------------------------
### read in mapped data per barcode
setwd(dir.map)
total.reads.files <- dir(pattern="total.reads.out")
mapped.NM.files <- dir(pattern="NM.out")
mapped.files <- sub("NM.","",mapped.NM.files)
read.lengths.files <- dir(pattern="_read.lengths.out")
aln.read.lengths.files <- dir(pattern=paste("aln.",small,".read.lengths.out",sep=""))

## total reads
total.reads.all <- NULL
for (i in 1:length(total.reads.files)){
  total.reads.all <- rbind(total.reads.all,get(load(total.reads.files[i])))
}
save(total.reads.all,file=paste(dir.res,"/obj_",experiment,"_total.reads_",small,".out",sep=""))
write.table(data.frame(Chip_barcode=rownames(total.reads.all),total.reads.all),file=paste(dir.res,"/",experiment,"_total.reads_",small,".txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE)

## mapped
mapped.all <- NULL
for (i in 1:length(mapped.files)){
  mapped.all <- cbind(mapped.all,get(load(mapped.files[i])))
}
save(mapped.all,file=paste(dir.res,"/obj_",experiment,"_mapped_",small,".out",sep=""))

## mapped NM
for (i in 1:length(mapped.NM.files)){
  write.csv(get(load(mapped.NM.files[i])),file=paste(dir.res, "/",substring(mapped.NM.files[i],5,nchar(mapped.NM.files[i])-4),".csv", sep=""))
}

## read lengths
read.lengths.all <- NULL
for (i in 1:length(read.lengths.files)){
  tmp <- strsplit(read.lengths.files[i],"_")[[1]]
  listname <- paste(tmp[2],tmp[3],sep="_")
  read.lengths.all[[listname]] <- get(load(read.lengths.files[i]))
}
save(read.lengths.all,file=paste(dir.res,"/obj_",experiment,"_read.lengths_",small,".out",sep=""))

## aln read lengths
aln.read.lengths.all <- NULL
for (i in 1:length(aln.read.lengths.files)){
  tmp <- strsplit(aln.read.lengths.files[i],"_")[[1]]
  listname <- paste(tmp[2],tmp[3],sep="_")
  aln.read.lengths.all[[listname]] <- get(load(aln.read.lengths.files[i]))
}
save(aln.read.lengths.all,file=paste(dir.res, "/obj_", experiment, "_aln.read.lengths_", small, ".out",sep=""))

#---------------------------------------------------------------------------------------------------------
# write collected data to CountTable
#---------------------------------------------------------------------------------------------------------
setwd(dir.deseq)
write.table(cbind(id=rownames(mapped.all),mapped.all), file=paste("CountTable_",small,".txt",sep=""), row.names=F, quote=F, sep="\t")


