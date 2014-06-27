#===============================================================
# Create piRNA database
#===============================================================
#
#  from all text files downloaded from piRNABank create single FASTA file
#
#------------------------------------------------------------------

library(Biostrings)
library(ShortRead)

organisms <- list(Zebrafish="dre", Mouse="mmu", Human="hsa", Arabidopsis="ath", Worm="cel", Drosophila="dme", Rat="rno")

org <- "Drosophila"

  dir.db.base <- "/mad/MAD-RBAB/05_Reference-db"
  dir.db <- paste(dir.db.base,"external",organisms[[org]],"RNA/piRNA",sep="/")
  dir.piRNAdb <- paste(dir.db,org,sep="/")

  setwd(dir.piRNAdb)

  pi.files <- dir()
  j <- 0
  file.exists <- FALSE
  DNA <- TRUE
  #for each piRNA file extract unique sequences and add to new piRNA file
  for (i in 1:length(pi.files)){
    j <- j+1
    if (i==1){ 
      DNA <- tryCatch({
        readDNAStringSet(pi.files[i])
        DNA <- TRUE
      }, warning=function(w){FALSE})  
    }
    if (DNA) piseq <- readDNAStringSet(pi.files[i]) else piseq <- DNAStringSet(readRNAStringSet(pi.files[i]))
    pi.new <- piseq[!srduplicated(piseq)]
    pi.nam <- paste(substring(names(pi.new),1,14),1:length(pi.new),sep=":")
    names(pi.new) <- pi.nam
    if (j==1) {
      pi.unique <- pi.new
      print("reset pi.unique")
    } else pi.unique <- c(pi.unique,pi.new)
    if (j==1000){
      if (file.exists) {
        writeXStringSet(pi.unique,paste(dir.db,"piRNA.fa",sep="/"),format="fasta",append=TRUE)
        print("fasta appended")
      } else {
        writeXStringSet(pi.unique,paste(dir.db,"piRNA.fa",sep="/"),format="fasta")
        file.exists <- TRUE
        print("fasta written")
      }
      j <- 0
    }
  }
  writeXStringSet(pi.unique,paste(dir.db,"piRNA.fa",sep="/"),format="fasta",append=TRUE)
  print("last fasta appended")

  setwd(dir.db)
  piseq <- readDNAStringSet("piRNA.fa")  
  print(piseq)
  print(summary(width(piseq)))

#---------------------------------------------------------
#Zebrafish
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  16.00   26.00   27.00   27.01   28.00   35.00 

#Mouse
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  16.00   28.00   29.00   28.89   30.00   38.00 

#Human
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  26.00   27.00   29.00   28.69   30.00   32.00 

#Rat
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    18.00   27.00   29.00   28.41   30.00   32.00   

#Fly
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  15.00   24.00   25.00   25.25   26.00   31.00 

