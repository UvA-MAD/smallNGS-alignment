#!/usr/bin/env Rscript
#### Collect results from miRNA aligmnents of different barcodes
#### 2014-05-15 Inez Terpstra
#### mature miRNA reads counts are combined per experiment
#----------------------------------------------------------------------------------------------------
print("Collecting mature miRNA count data")
#----------------------------------------------------------
args <- commandArgs(TRUE)

base.dir <- args[1]
experiment <- args[2]

#==========================================================================================
# environment
#---------------------------------------------------------------------------------------------------------
dir.map <- paste(base.dir,"/Scratch/",experiment,"/Data/miRNA_aligned_data",sep="") 
dir.res <- paste(base.dir,"/Results/",experiment,"/Output",sep="") 
dir.deseq <- paste(base.dir,"/Results/",experiment,"/DESeq",sep="")

#---------------------------------------------------------------------------------------------------------
### read in mapped data per barcode
setwd(dir.map)
mapped.files <- dir(pattern="mature.out")
isomir.files <- dir(pattern=".part")
nr.runs <- length(isomir.files)

## mapped mature
mapped.all <- NULL
for (i in 1:nr.runs){
  mapped.all <- cbind(mapped.all,get(load(mapped.files[i])))
}
save(mapped.all,file=paste(dir.res,"/obj_",experiment,"_mapped_mature.out",sep=""))

## mapped isomirs - collect all in 1 list
hairpin.parts <- NULL
hairpin.part.names <- NULL
for (i in 1:nr.runs){
  load(file=isomir.files[i])
  barcode <- strsplit(isomir.files[i],"_")[[1]][3]
  chip <- strsplit(isomir.files[i],"_")[[1]][2]
  run.name <- paste(chip,barcode,sep="_") 
  hairpin.parts[[run.name]] <- mapped.hairpin.parts[,c(1:2,4:6,3)]
  hairpin.part.names <- c(hairpin.part.names,with(hairpin.parts[[run.name]],paste(ID,part,start,stop,len)))
}

#length(hairpin.part.names)
#length(unique(hairpin.part.names))

hairpin.merge <- hairpin.parts[[1]]
run.names <- names(hairpin.parts)
colnames(hairpin.merge) <- c(head(colnames(hairpin.parts[[1]]),-1),run.names[1])
if (nr.runs>1){
  for (run.nr in 2:nr.runs){
    merge.names <- colnames(hairpin.merge)
    hairpin.merge <- merge(hairpin.merge,hairpin.parts[[run.nr]],by=c("ID", "part", "start", "stop", "len"), all.x=TRUE, all.y=TRUE)
    colnames(hairpin.merge) <- c(merge.names,run.names[run.nr])
    print(dim(hairpin.merge))
  }
}

combID <- function(x) paste(x[1],x[2],as.integer(x[3]),as.integer(x[4]),as.integer(x[5]),sep="_")
isomir.id <- apply(hairpin.merge[,1:5],1,combID)
rownames(hairpin.merge) <- isomir.id

#replace NAs with 0
hairpin.merge[is.na(hairpin.merge)] <- 0

save(hairpin.merge,file=paste(dir.res,"/obj_",experiment,"_mapped.part_hairpin.out",sep=""))
#load(file=paste(dir.res,"/obj_",experiment,"_mapped.part_hairpin.out",sep=""))

write.csv(hairpin.merge,file=paste(dir.res,"/",experiment,"_mapped_isomirs.csv",sep=""))

#---------------------------------------------------------------------------------------------------------
# write collected data to CountTable
#---------------------------------------------------------------------------------------------------------
setwd(dir.deseq)
write.table(cbind(id=rownames(mapped.all),mapped.all), file="CountTable_mature.txt", row.names=F, quote=F, sep="\t")
write.table(cbind(id=rownames(hairpin.merge),hairpin.merge), file="CountTable_isomirs.txt", row.names=F, quote=F, sep="\t")


