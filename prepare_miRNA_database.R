#### Extract organism specific sequences from miRNA hairpin database (miRBase20).

library(Biostrings)
library(ShortRead)

#----------------------------------------------------------------------------------------------------------
#location of mirBase 20 miRNA data
dir.db.base="/mad/MAD-RBAB/05_Reference-db"
dir.db <- paste(dir.db.base,"external",org,"RNA/miRNA",sep="/")
dir.miRNAdb <- paste(dir.db.base,"external/RNA/mirBase20",sep="/")

#-----------------------------------------------------------------------------------
##### load miRNA data from miRBase20
setwd(dir.miRNAdb)

#-----------------------------------------------------------------------------------
# load objects:
load("obj_premiRNA.out")
load("obj_miRNA.out")

pre.miRNA.tmp <- pre.miRNA
miRNA.tmp <- miRNA

#----------------------------------------------------------------------------------------------------------
# retrieve reference - mature and hairpin miRNA sequences
# fasta file should be in DNA format for mapping with SHRiMP2
mature.RNA <- readRNAStringSet("mature.fa",format="fasta")
hairpin.RNA <- readRNAStringSet("hairpin.fa",format="fasta")
# done --
writeXStringSet(DNAStringSet(hairpin.RNA),file="hairpin20.fa",format="fasta")
writeXStringSet(DNAStringSet(mature.RNA),file="mature20.fa",format="fasta")

mature <- readDNAStringSet("mature20.fa",format="fasta")
hairpin <- readDNAStringSet("hairpin20.fa",format="fasta")

mature.tmp <- mature
hairpin.tmp <- hairpin

mature.ID <- sapply(strsplit(names(mature.tmp)," "),"[[",1)
hairpin.ID <- sapply(strsplit(names(hairpin.tmp)," "),"[[",1)

#----------------------------------------------------------------------------------------------------------
### select miRNA data of organism of interest 

#zebrafish
org <- "dre"
pre.miRNA <- pre.miRNA.tmp[grep(org,pre.miRNA.tmp$ID),]
miRNA <- miRNA.tmp[grep(org,miRNA.tmp$ID),]
save(pre.miRNA,file="obj_premiRNAdre.out")
save(miRNA,file="obj_miRNAdre.out")
mature <- mature.tmp[grep(org,mature.ID)]
hairpin <- hairpin.tmp[ grep(org,hairpin.ID)]
writeXStringSet(DNAStringSet(mature), file=paste(dir.db.base, "external", org, "RNA/miRNA/mature20.fa",sep="/"), format="fasta")
writeXStringSet(DNAStringSet(hairpin), file=paste(dir.db.base, "external", org, "RNA/miRNA/hairpin20.fa",sep="/"), format="fasta")
save(hairpin,mature,pre.miRNA,miRNA, file=paste(dir.db.base,"external",org,"RNA/miRNA/miRNAdata.RData",sep="/"))
hairpin.id <- sapply(strsplit(names(hairpin)," "),"[[",1)
if (length(which(pre.miRNA$ID != hairpin.id)) > 0) print(paste(org,"check hairpin an pre.miRNA, id discrepancy!!",sep=" "))

#mouse
org <- "mmu"
pre.miRNA <- pre.miRNA.tmp[grep(org,pre.miRNA.tmp$ID),]
miRNA <- miRNA.tmp[grep(org,miRNA.tmp$ID),]
save(pre.miRNA,file="obj_premiRNAmmu.out")
save(miRNA,file="obj_miRNAmmu.out")
mature <- mature.tmp[grep(org,mature.ID)]
hairpin <- hairpin.tmp[ grep(org,hairpin.ID)]
writeXStringSet(DNAStringSet(mature), file=paste(dir.db.base, "external", org, "RNA/miRNA/mature20.fa",sep="/"), format="fasta")
writeXStringSet(DNAStringSet(hairpin), file=paste(dir.db.base, "external", org, "RNA/miRNA/hairpin20.fa",sep="/"), format="fasta")
save(hairpin,mature,pre.miRNA,miRNA, file=paste(dir.db.base,"external",org,"RNA/miRNA/miRNAdata.RData",sep="/"))
hairpin.id <- sapply(strsplit(names(hairpin)," "),"[[",1)
if (length(which(pre.miRNA$ID != hairpin.id)) > 0) print(paste(org,"check hairpin an pre.miRNA, id discrepancy!!",sep=" "))

#rat
org <- "rno"
pre.miRNA <- pre.miRNA.tmp[grep(org,pre.miRNA.tmp$ID),]
miRNA <- miRNA.tmp[grep(org,miRNA.tmp$ID),]
save(pre.miRNA,file="obj_premiRNArno.out")
save(miRNA,file="obj_miRNArno.out")
mature <- mature.tmp[grep(org,mature.ID)]
hairpin <- hairpin.tmp[ grep(org,hairpin.ID)]
writeXStringSet(DNAStringSet(mature), file=paste(dir.db.base, "external", org, "RNA/miRNA/mature20.fa",sep="/"), format="fasta")
writeXStringSet(DNAStringSet(hairpin), file=paste(dir.db.base, "external", org, "RNA/miRNA/hairpin20.fa",sep="/"), format="fasta")
save(hairpin,mature,pre.miRNA,miRNA, file=paste(dir.db.base,"external",org,"RNA/miRNA/miRNAdata.RData",sep="/"))
hairpin.id <- sapply(strsplit(names(hairpin)," "),"[[",1)
if (length(which(pre.miRNA$ID != hairpin.id)) > 0) print(paste(org,"check hairpin an pre.miRNA, id discrepancy!!",sep=" "))

#human
org <- "hsa"
pre.miRNA <- pre.miRNA.tmp[grep(org,pre.miRNA.tmp$ID),]
miRNA <- miRNA.tmp[grep(org,miRNA.tmp$ID),]
save(pre.miRNA,file="obj_premiRNAhsa.out")
save(miRNA,file="obj_miRNAhsa.out")
mature <- mature.tmp[grep(org,mature.ID)]
hairpin <- hairpin.tmp[ grep(org,hairpin.ID)]
writeXStringSet(DNAStringSet(mature), file=paste(dir.db.base, "external", org, "RNA/miRNA/mature20.fa",sep="/"), format="fasta")
writeXStringSet(DNAStringSet(hairpin), file=paste(dir.db.base, "external", org, "RNA/miRNA/hairpin20.fa",sep="/"), format="fasta")
save(hairpin,mature,pre.miRNA,miRNA, file=paste(dir.db.base,"external",org,"RNA/miRNA/miRNAdata.RData",sep="/"))
hairpin.id <- sapply(strsplit(names(hairpin)," "),"[[",1)
if (length(which(pre.miRNA$ID != hairpin.id)) > 0) print(paste(org,"check hairpin an pre.miRNA, id discrepancy!!",sep=" "))

#fly
org <- "dme"
pre.miRNA <- pre.miRNA.tmp[grep(org,pre.miRNA.tmp$ID),]
miRNA <- miRNA.tmp[grep(org,miRNA.tmp$ID),]
save(pre.miRNA,file="obj_premiRNAdme.out")
save(miRNA,file="obj_miRNAdme.out")
mature <- mature.tmp[grep(org,mature.ID)]
hairpin <- hairpin.tmp[ grep(org,hairpin.ID)]
writeXStringSet(DNAStringSet(mature), file=paste(dir.db.base, "external", org, "RNA/miRNA/mature20.fa",sep="/"), format="fasta")
writeXStringSet(DNAStringSet(hairpin), file=paste(dir.db.base, "external", org, "RNA/miRNA/hairpin20.fa",sep="/"), format="fasta")
save(hairpin,mature,pre.miRNA,miRNA, file=paste(dir.db.base,"external",org,"RNA/miRNA/miRNAdata.RData",sep="/"))
hairpin.id <- sapply(strsplit(names(hairpin)," "),"[[",1)
if (length(which(pre.miRNA$ID != hairpin.id)) > 0) print(paste(org,"check hairpin an pre.miRNA, id discrepancy!!",sep=" "))

#worm
org <- "cel"
pre.miRNA <- pre.miRNA.tmp[grep(org,pre.miRNA.tmp$ID),]
miRNA <- miRNA.tmp[grep(org,miRNA.tmp$ID),]
save(pre.miRNA,file="obj_premiRNAcel.out")
save(miRNA,file="obj_miRNAcel.out")
mature <- mature.tmp[grep(org,mature.ID)]
hairpin <- hairpin.tmp[ grep(org,hairpin.ID)]
writeXStringSet(DNAStringSet(mature), file=paste(dir.db.base, "external", org, "RNA/miRNA/mature20.fa",sep="/"), format="fasta")
writeXStringSet(DNAStringSet(hairpin), file=paste(dir.db.base, "external", org, "RNA/miRNA/hairpin20.fa",sep="/"), format="fasta")
save(hairpin,mature,pre.miRNA,miRNA, file=paste(dir.db.base,"external",org,"RNA/miRNA/miRNAdata.RData",sep="/"))
hairpin.id <- sapply(strsplit(names(hairpin)," "),"[[",1)
if (length(which(pre.miRNA$ID != hairpin.id)) > 0) print(paste(org,"check hairpin an pre.miRNA, id discrepancy!!",sep=" "))

#arabidopsis
org <- "ath"
pre.miRNA <- pre.miRNA.tmp[grep(org,pre.miRNA.tmp$ID),]
miRNA <- miRNA.tmp[grep(org,miRNA.tmp$ID),]
save(pre.miRNA,file="obj_premiRNAath.out")
save(miRNA,file="obj_miRNAath.out")
mature <- mature.tmp[grep(org,mature.ID)]
hairpin <- hairpin.tmp[ grep(org,hairpin.ID)]
writeXStringSet(DNAStringSet(mature), file=paste(dir.db.base, "external", org, "RNA/miRNA/mature20.fa",sep="/"), format="fasta")
writeXStringSet(DNAStringSet(hairpin), file=paste(dir.db.base, "external", org, "RNA/miRNA/hairpin20.fa",sep="/"), format="fasta")
save(hairpin,mature,pre.miRNA,miRNA, file=paste(dir.db.base,"external",org,"RNA/miRNA/miRNAdata.RData",sep="/"))
hairpin.id <- sapply(strsplit(names(hairpin)," "),"[[",1)
if (length(which(pre.miRNA$ID != hairpin.id)) > 0) print(paste(org,"check hairpin an pre.miRNA, id discrepancy!!",sep=" "))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#check occurrences of bases in hairpins sequences
colSums(alphabetFrequency(hairpin))
A    C    G    U    M    R    W    S    Y    K    V    H    D    B    N    -    + 
  7548 6757 8053 9840    0    0    0    0    0    0    0    0    0    0    0    0    0 


