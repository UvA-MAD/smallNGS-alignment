#===============================================================
# Create database
#===============================================================
#
#  from fasta files create single new file formatted properly (CAPS and short lines)
#
#------------------------------------------------------------------

library(Biostrings)
library(ShortRead)

dir.db.base <- "/mad/MAD-RBAB/05_Reference-db"

## convert human rRNA database file from Rob (20140614) hs_45S_rDNA.fasta
org <- "hsa"
small <- "rRNA"
dir.db <- paste(dir.db.base,"external",org,"RNA",small,sep="/")
setwd(dir.db)
fasta.files <- dir(pattern="fasta")

v <- "1"
rRNA <- readDNAStringSet(fasta.files[2])
writeXStringSet(rRNA,paste(dir.db,"/",small,"_",v,".fa",sep=""),format="fasta")

## convert human rRNA database files from BioMart (20140623) BioMart_rRNA.fasta and BioMart_Mt_rRNA.fasta
org <- "hsa"
small <- "rRNA"
dir.db <- paste(dir.db.base,"external",org,"RNA",small,sep="/")
setwd(dir.db)
fasta.files <- dir(pattern="BioMart")

v <- "bm20140623"
mt_rRNA <- readDNAStringSet(fasta.files[1])
rRNA <- readDNAStringSet(fasta.files[2])
rRNA <- c(rRNA,mt_rRNA)
writeXStringSet(rRNA,paste(dir.db,"/",small,"_",v,".fa",sep=""),format="fasta")


