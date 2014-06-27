#### Prepare miRNA database R objects per organism from miRBase 20
#### 2014-05-21 Inez Terpstra
#----------------------------------------------------------------------------------------------------

library(Biostrings)
library(ShortRead)

#----------------------------------------------------------------------------------------------------------
#location of mirBase 20 miRNA data
dir.db.base="/mad/MAD-RBAB/05_Reference-db"
dir.miRNAdb <- paste(dir.db.base,"external/RNA/mirBase20",sep="/")

#----------------------------------------------------------------------------------------------------------
setwd(dir.miRNAdb)
miRNA.dat <- readLines("miRNA.dat")

#for each pre-miRNA list mature forms
pre.miRNA <- data.frame(ID=NA, acc=NA,Desc=NA, SQ=NA,LEN=NA)
miRNA <- NULL

#end of each miRNA marked with //
idx <- grep("^//",miRNA.dat)
idx <- c(1,idx)

#for each pre-miRNA extract mature miRNA info
for (j in 1:(length(idx)-1)){
  pre.miRNA.tmp <- data.frame(ID=NA, acc=NA,Desc=NA, SQ=NA,LEN=NA,start=NA,stop=NA)
  ss.miRNA <- miRNA.dat[idx[j]:idx[j+1]]
  #not all miRNA entries have ID or ACC
  if (length(grep("^ID",ss.miRNA))>0){
    pre.miRNA.tmp$ID <- strsplit(ss.miRNA[grep("^ID",ss.miRNA)]," +")[[1]][2]
    pre.miRNA.tmp$LEN <- as.double(tail(strsplit(sub(" BP.+","",ss.miRNA[grep("^ID",ss.miRNA)]),"; ")[[1]],1))
  }
  if(nchar(pre.miRNA.tmp$ID)>40) print(j)
  if (length(grep("^AC",ss.miRNA))>0){ 
    pre.miRNA.tmp$acc <- strsplit(ss.miRNA[grep("^AC",ss.miRNA)],"   ")[[1]][2]
    pre.miRNA.tmp$acc <- sub(";","",pre.miRNA.tmp$acc)
  }
  #not all miRNA entries have SQ
  i <- grep("^SQ +Sequence",ss.miRNA,value=FALSE)
  if (length(i) > 0){
    sq <- NULL
    for (s in (i+1):(length(ss.miRNA)-1)){
      sq <- paste(sq,ss.miRNA[s])
    }
    #clean up sequence from spaces and numbers
    a <- unlist(strsplit(sq," "))
    pre.miRNA.tmp$SQ <- paste(a[grep("[a-z]",a)],collapse="")
    #length between Sequence and BP
    pre.miRNA.tmp$LEN <- as.double(sub("^SQ +Sequence ","",sub(" BP.+","",ss.miRNA[i])))
  }
  if (length(grep("^DE",ss.miRNA))>0) pre.miRNA.tmp$Desc <-  strsplit(ss.miRNA[grep("^DE",ss.miRNA)],"   ")[[1]][2]
  
  #get mature info
  FT.miRNA <- ss.miRNA[grep("^FT   miRNA",ss.miRNA)]
  #not all miRNA have mature info
  m <- length(FT.miRNA)
  if (m > 0) for (i in 1:m){
    miRNA.tmp <- data.frame(ID=NA, mature=NA,acc=NA,start=NA,stop=NA)
    miRNA.tmp$ID <- pre.miRNA.tmp$ID
    miRNA.tmp$mature <-  strsplit(tail(strsplit(ss.miRNA[grep("/product",ss.miRNA)]," ")[[i]],1),"\\\"")[[1]][2]
    miRNA.tmp$acc <-  strsplit(tail(strsplit(ss.miRNA[grep("/accession",ss.miRNA)]," ")[[i]],1),"\\\"")[[1]][2]  
    miRNA.tmp$start <- as.double(strsplit(tail(strsplit(FT.miRNA," ")[[i]],1),"\\..")[[1]][1])
    miRNA.tmp$stop <- as.double(strsplit(tail(strsplit(FT.miRNA," ")[[i]],1),"\\..")[[1]][2])
    miRNA <- rbind(miRNA,miRNA.tmp)
  }
  pre.miRNA[j,] <- pre.miRNA.tmp
}

		#check entries:

		sum(is.na(pre.miRNA$ID))
		sum(is.na(pre.miRNA$acc))
		sum(is.na(miRNA$ID))
		sum(is.na(miRNA$mature))
		sum(is.na(miRNA$acc))
		miRNA[grep("BP.",miRNA$ID),]

		table(pre.miRNA$LEN)
		table(miRNA$start)
		table(miRNA$stop)

		#lookup a specific miRNA:
		j <- grep("cpa-MIR162a",miRNA.dat[idx+1])
		k <- grep("cpa-miR162a",miRNA.dat)
		miRNA.dat[(k[1]-50):(k[1]+50)]

# save objects:
save(pre.miRNA,file="obj_premiRNA.out")
load("obj_premiRNA.out")

save(miRNA,file="obj_miRNA.out")
load("obj_miRNA.out")



