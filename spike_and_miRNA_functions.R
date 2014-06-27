#### Functions for processing aligmnents
#### 2014-05-12 Inez Terpstra
#### getBam: to read in bam file and convert it to data.frame
#### correctBamNM: to get gap extension of InDels
#### addBamAln: to set filters for 10% mismatch and 3p overhang
#### hitsNM: to split read counts on number of mismatches
#### assignMature: to assign aligned reads to mature or hairpin
#### count.isomiRs: to count number of aligend reads per isoform
#### getMappingData: to wrap functions to get valid alignments per run and save read counts
#==========================================================================================
#read in bam file
#-----------------------------------------------------------------------------------------
getBam <- function(bam.file){
print("reading")
  what <- scanBamWhat()
  p1 <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery=FALSE), simpleCigar = FALSE,tag = c("AS","NM"),what = what)
  bam <- scanBam(bam.file,type="BAM",param=p1)[[1]]
  bam.df <- data.frame(qname=bam$qname,flag=bam$flag,rname=bam$rname,
                       strand=bam$strand,pos=bam$pos,qwidth=bam$qwidth,
                       mapq=bam$mapq,cigar=bam$cigar,mrnm=bam$mrnm,  
                       mpos=bam$mpos,isize=bam$isize,seq=as.character(bam$seq),
                       qual=as.character(bam$qual),AS=bam$tag[[1]],NM=bam$tag[[2]])
  bam.df
}

#==========================================================================================
#split cigar string on InDel to obtain gap extensions
#-----------------------------------------------------------------------------------------
correctBamNM <- function(cigar){
  NMcorrection <- sapply(
    lapply(
      strsplit(gsub("[0-9]+M","",cigar),"I|D"),
      as.numeric
    ),
    function(x) {sum(x)-length(x)}
  )
  NMcorrection
}

#==========================================================================================
#get alignment length local mapping -- extract nr of 5'and 3'softclipped nucleotides from cigar string
#-----------------------------------------------------------------------------------------
addBamAln <- function(bam,aln.mode,perc.mm,un5p,un3p,fgp){
  
  print("clipping")

  aln.txt <- ifelse (aln.mode=="local", "S","I")
  aln.txt.1 <- paste("^[0-9]+",aln.txt,"[0-9]",sep="")
  aln.txt.2 <- paste(aln.txt,".*",sep="")
  aln.txt.3 <- paste(".+",aln.txt,"$",sep="")
  aln.txt.4 <- paste(aln.txt,"$",sep="")
  idx.S5p <- grep(aln.txt.1,bam$cigar)
  bam$S5p <- NA
  bam$S5p[idx.S5p] <- as.numeric(sub(aln.txt.2,"",bam$cigar[idx.S5p]))
  idx.S3p <- grep(aln.txt.3,bam$cigar)
  bam$S3p <- NA
  bam$S3p[idx.S3p] <-  as.numeric(sub(".*[^0-9]","", sub(aln.txt.4,"",bam$cigar[idx.S3p])))
  bam$softclipped <- rowSums(cbind(bam$S5p,bam$S3p),na.rm=TRUE)
  bam$aln <- bam$qwidth - bam$softclipped
  #edit distance has to be corrected for gap extensions
  #extract indel information from cigar string
  #strip cigar string first from 5p softclipped bases
  bam$cigar5p <- NA
  bam$cigar5p[is.na(bam$S5p)] <- as.character(bam$cigar[is.na(bam$S5p)])
  bam$cigar5p[is.na(bam$S5p)==FALSE] <- substr(as.character(bam$cigar[is.na(bam$S5p)==FALSE]),nchar(as.character(bam$S5p[is.na(bam$S5p)==FALSE]))+2,nchar(as.character(bam$cigar[is.na(bam$S5p)==FALSE])))
  #strip 5p stripped cigar string from 3p softclipped bases
  bam$cigar3p <- NA
  bam$cigar3p[is.na(bam$S3p)] <- as.character(bam$cigar5p[is.na(bam$S3p)])
  bam$cigar3p[is.na(bam$S3p)==FALSE] <- substr(as.character(bam$cigar5p[is.na(bam$S3p)==FALSE]),1,nchar(as.character(bam$cigar5p[is.na(bam$S3p)==FALSE]))-nchar(as.character(bam$S3p[is.na(bam$S3p)==FALSE]))-1)
  if (fgp == FALSE) bam$NM.corr <- bam$NM
  #-----------------------------------------------------------------------------------------
  #set filters for alignments
  #-----------------------------------------------------------------------------------------
  #filter 1 - number of mismatches and indels < floor percmm of alignment length
  if (fgp) bam$filter1 <- bam$NM <= floor(bam$aln/perc.mm)
  else bam$filter1 <- bam$NM.corr <= floor(bam$aln/perc.mm)
  #filter 2 - number of unmatched nucleotides on 3p side is <= 3
  bam$filter3 <- TRUE
  bam$filter3[is.na(bam$S3p)==FALSE] <- bam$S3p[is.na(bam$S3p)==FALSE] <= un3p
  #filter 3 - nonaligned stretches more than 3 ns --> discard
  bam$filter5 <- TRUE
  bam$filter5[is.na(bam$S5p)==FALSE] <- bam$S5p[is.na(bam$S5p)==FALSE] <= un5p
  #all filters true for valid alignment
  bam$mapped <-  rowSums(cbind(bam$filter1,bam$filter3,bam$filter5),na.rm=TRUE)>=3
  bam
}

#==========================================================================================
# list number of hits per corrected edit distance (zero or more)
#-----------------------------------------------------------------------------------------
hitsNM <- function(rname,NM.corr,pre.miRNA){
  NM <- table(NM.corr)
  mapped.tmp <- matrix(0,dim(pre.miRNA)[1],length(NM))
  colnames(mapped.tmp) <- paste("NM",names(NM),sep="")
  rownames(mapped.tmp) <- pre.miRNA$ID
  for (j in names(NM)){
    mapped.tmp[names(table(rname[NM.corr==j])),paste("NM",j,sep="")] <- table(rname[NM.corr==j])
  }
  mapped.tmp
}

#==========================================================================================
# assign valid aligned reads to mature or pre-miRNA
#-----------------------------------------------------------------------------------------
assignMature <- function(mapped.hp,miRNA,ovl){
  mapped.hp$part <- "hairpin"
  for (i in 1:dim(mapped.hp)[1]){ 
    idx <- which(miRNA$ID==mapped.hp$ID[i])
    for (j in 1:length(idx)){
      if (mapped.hp$part[i]=="hairpin"){
        if (mapped.hp$start[i] >= (miRNA$start[idx][j]-ovl) & mapped.hp$stop[i] <= (miRNA$stop[idx][j]+ovl)) mapped.hp$part[i] <- miRNA$mature[idx][j]
      }
    }
  }
  mapped.hp
}

#==========================================================================================
# count number of aligned reads per isoform
#-----------------------------------------------------------------------------------------
count.isomiRs <- function(rname,start.pos,aln,miRNA,ovl){
  isoforms <- table(paste(rname,start.pos,aln))
  premiR <- sapply(strsplit(names(isoforms)," "),"[[",1)
  str.pos <- as.double(sapply(strsplit(names(isoforms)," "),"[[",2))
  aln.len <- as.double(sapply(strsplit(names(isoforms)," "),"[[",3))
  end.pos <- str.pos+aln.len-1
  mapped.hairpin.parts <- data.frame(ID=premiR,part="hairpin", aln.reads=as.vector(isoforms),start=str.pos,stop=end.pos,len=aln.len,stringsAsFactors=FALSE)
  for (i in 1:dim(mapped.hairpin.parts)[1]){ 
    idx <- which(miRNA$ID==mapped.hairpin.parts$ID[i])
    for (j in 1:length(idx)){
      if (mapped.hairpin.parts$part[i]=="hairpin"){
        if (mapped.hairpin.parts$start[i] >= (miRNA$start[idx][j]-ovl) 
            & mapped.hairpin.parts$stop[i] <= (miRNA$stop[idx][j]+ovl)) {
          mapped.hairpin.parts$part[i] <- miRNA$mature[idx][j]
        }
      }
    }
  }
  mapped.hairpin.parts
}

#==========================================================================================
# process mapping results for each run
#-----------------------------------------------------------------------------------------
getMappingData <- function(run.name, bam.file.local, bam.file.global, percent.mismatches, unmatched5p, unmatched3p, allowMatureExtension, pre.miRNA, miRNA, fgp){
  
  print("getMappingData")

  bam.glob <- getBam(bam.file.global)
  bam.glob <- addBamAln(bam.glob,"global",percent.mismatches,unmatched5p,unmatched3p, fgp)
  
  mapped.glob <- bam.glob$qname[bam.glob$mapped]
  bam.nm <- bam.glob[bam.glob$mapped,]
  rm(bam.glob)
  
  bam.loc <- getBam(bam.file.local)
  bam.loc <- addBamAln(bam.loc,"local",percent.mismatches,unmatched5p,unmatched3p, fgp)

  # combine local and global mapped reads into bam.NM -- global mapping results are preferred
  mapped.loc <- bam.loc$qname[bam.loc$mapped]
  mapped.loc <- setdiff(mapped.loc,mapped.glob)
  
  bam.nm <- rbind(bam.loc[bam.loc$qname %in% mapped.loc, ],bam.nm)
  rm(bam.loc)
  
  #count number of aligned reads per isoform
  mapped.hairpin.parts <- count.isomiRs(bam.nm$rname,bam.nm$pos,bam.nm$aln,miRNA,allowMatureExtension)
  
  #create result matrix nr of mapped (pre-)miRNAs
  mapped.mature <- matrix(0,length(unique(miRNA$mature)),1)
  rownames(mapped.mature) <- unique(miRNA$mature)
  colnames(mapped.mature) <- run.name
  mapped.hairpin <- matrix(0,length(pre.miRNA$ID),1)
  rownames(mapped.hairpin) <- pre.miRNA$ID
  colnames(mapped.hairpin) <- run.name
  
  #count number of hits per pre-miRNA
  mapped.hairpin[names(table(bam.nm$rname)),1] <- table(bam.nm$rname)

  #for each mapped mature part of a pre-miRNA return mapped mature form
  for (m in 1:dim(mapped.mature)[1]){
    mapped.mature[m,1] <- sum(as.double(mapped.hairpin.parts$aln.reads[mapped.hairpin.parts$part %in% rownames(mapped.mature)[m]]))
  }
  
  #count number of hits per pre-miRNA and per edit distance
  if (fgp) mapped.hairpin.NM <- hitsNM(bam.nm$rname,bam.nm$NM,pre.miRNA)
  else mapped.hairpin.NM <- hitsNM(bam.nm$rname,bam.nm$NM.corr,pre.miRNA)
 
  #save results
  save(mapped.mature,file=paste("obj",run.name,"mapped.mature.out",sep="_"))
  save(mapped.hairpin,file=paste("obj",run.name,"mapped.hairpin.out",sep="_"))
  save(mapped.hairpin.NM,file=paste("obj",run.name,"mapped.hairpin.NM.out",sep="_")) 
  save(mapped.hairpin.parts, file=paste("obj",run.name,"mapped.hairpin.part.out",sep="_"))
  
  bam.nm
}

