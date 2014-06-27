#!/usr/bin/env Rscript
#### Normalize results 
#### 2014-05-15 Inez Terpstra
#### all alignments are normalized on all counts and on spikes only, if spikes are mapped
#----------------------------------------------------------
args <- commandArgs(TRUE)

base.dir <- args[1]
experiment <- args[2]
organism <- args[3]

#----------------------------------------------------------------------------------------------------
print("Normalizing count data")

#==========================================================================================
# environment
#---------------------------------------------------------------------------------------------------------
dir.res <- paste(base.dir,"/Results/",experiment,"/Output",sep="") 
dir.deseq <- paste(base.dir,"/Results/",experiment,"/DESeq",sep="")
dir.seq <- paste(base.dir,"/DataSeq",sep="")
dir.plot <- paste(base.dir,"/Images/",experiment,sep="") 
dir.db.base <- "/mad/MAD-RBAB/05_Reference-db"
dir.src <- "/mad/MAD-RBAB/10_Protocols/DryLab/smallNGS_alignment"

#----------------------------------------------------------
#suppressPackageStartupMessages(library(ShortRead))
library(calibrate)
#----------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------
### read in Design
#---------------------------------------------------------------------------------------------------------
setwd(dir.seq)
run.design <- read.table(file="design.txt",sep="\t",header=TRUE)
run.design <- run.design[with(run.design,order(Chip,Barcode)),]

#---------------------------------------------------------------------------------------------------------
### read in CountTables
#---------------------------------------------------------------------------------------------------------
setwd(dir.deseq)
count.files <- dir(pattern="CountTable")

#==========================================================================================
# functions
#---------------------------------------------------------------------------------------------------------
source(paste(dir.src,"normalization_functions.R",sep="/"))
source(paste(dir.src,"plot_functions.R",sep="/"))

#-------------------------------------------------------------------------------------
#### I - spikes
## check if spikes are present:
idx <- grep("spike",count.files)
if (length(idx) > 0){
  spikes <- TRUE
  spike.data <- read.delim(paste(dir.deseq,count.files[idx],sep="/"), row.names=1)
  spike.idx <- grep("NSPK",rownames(spike.data))                  
  sf.spike <- sizeFactors.mad(spike.data[spike.idx,], locfunc = median)

  #create annotation file
  x.id <- rownames(spike.data)
  anno.spike <- data.frame(id=x.id)
  write.table(anno.spike, paste(dir.deseq,"/Info/Annotation_spike.txt",sep=""), row.names=F, quote=F, sep="\t")
} else spikes <- FALSE

#-------------------------------------------------------------------------------------
#### small RNAs

if (spikes) count.files <- count.files[-idx]
for (filenr in 1:length(count.files)){
  data <- read.delim(paste(dir.deseq,count.files[filenr],sep="/"), row.names=1)
  small <- strsplit(count.files[filenr],"_")[[1]][2]
  small <- substring(small,1,nchar(small)-4)

  #create annotation file
  x.id <- rownames(data)
  if (small == "miRNA" | small == "mature"){
    load(paste(dir.db.base,"external",organism,"/RNA/miRNA/miRNAdata.RData",sep="/"))
    if (small == "mature"){
      x.hairpin <- NULL
      for (i in 1:length(x.id)){
        x.hairpin <- c(x.hairpin,miRNA$ID[miRNA$mature==x.id[i]][1])  
      }
      x.acc <- NULL
      for (i in 1:length(x.id)){
        x.acc <- c(x.acc,miRNA$acc[miRNA$mature==x.id[i]][1])  
      }
      anno <- data.frame(id=x.id,hairpin=x.hairpin,acc=x.acc)
    } else {
      x.acc <- NULL
      for (i in 1:length(x.id)){
        x.acc <- c(x.acc,pre.miRNA$acc[pre.miRNA$ID==x.id[i]][1])  
      }
      anno <- data.frame(id=x.id,acc=x.acc)
    }
  } else 
  if (small == "isomirs"){
    anno <- as.data.frame(data[,1:5])
    data <- data[,-c(1:5)]
  } else anno <- data.frame(id=x.id)
  write.table(anno, paste(dir.deseq,"/Info/Annotation_",small,".txt",sep=""), row.names=F, quote=F, sep="\t")

  if (spikes){
    #include spikes in dataset 
    all.data <- rbind(spike.data,data)
    nr.spike <- dim(anno.spike)[1]
    nr.small <- dim(anno)[1]
    nr.anno <- dim(anno)[2]
    all.anno <- data.frame(id=c(as.vector(anno.spike[,1]),as.vector(anno[,1])))
    if (nr.anno > 1) for (i in 2:nr.anno) all.anno <- cbind(all.anno,c(rep("",nr.spike),as.vector(anno[,i])))
    colnames(all.anno) <- names(anno)
    #normalized on spikes only      
    data.spike.norm <- t(t(all.data)/sf.spike)
    #export spike normalized data
    write.table(cbind(all.anno,data.spike.norm), file=paste(dir.deseq,"/",small,"_spike_norm.txt",sep=""),  row.names=FALSE, sep="\t")
    #plot boxplots 
    setwd(paste(dir.plot,"Box",sep="/"))
    Box.mad(small,data.spike.norm,"spike_norm",".")
    #plot MA plots
    setwd(paste(dir.plot,"MA",sep="/"))
    MA.mad(small,data.spike.norm,"spike_norm")
    #plot PCA plots
    setwd(paste(dir.plot,"PCA",sep="/"))
    PCA.mad(small,data.spike.norm,"spike.norm",run.design)
  } else {     
    all.data <- data
    all.anno <- anno
  }
   
  #normalize on all counts
  sf <- sizeFactors.mad(all.data, locfunc = median)
  data.norm <- t(t(all.data)/sf)
  #export raw and normalized data
  write.table(cbind(all.anno,all.data), file=paste(dir.deseq,"/",small,"_raw.txt",sep=""),  row.names=FALSE, sep="\t")
  write.table(cbind(all.anno,data.norm), file=paste(dir.deseq,"/",small,"_norm.txt",sep=""),  row.names=FALSE, sep="\t")
  #plot boxplots 
  setwd(paste(dir.plot,"Box",sep="/"))
  Box.mad(small,all.data,"raw",".")
  Box.mad(small,data.norm,"norm",".")
  #plot MA plots
  setwd(paste(dir.plot,"MA",sep="/"))
  MA.mad(small,all.data,"raw")
  MA.mad(small,data.norm,"norm")
  #plot PCA plots
  setwd(paste(dir.plot,"PCA",sep="/"))
  PCA.mad(small,all.data,"raw",run.design)
  PCA.mad(small,data.norm,"norm",run.design)
}


