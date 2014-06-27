#### Functions for processing aligmnents
#### 2014-05-19 Inez Terpstra
#### getBam: to read in bam file and convert it to data.frame
#### correctBamNM: to get gap extension of InDels
#### addBamAln: to set filters for 10% mismatch and 3p overhang
#### hitsNM: to split read counts on number of mismatches
#### getMappingData: to warp functions to get valid alignments per run and save read counts

#==========================================================================================
# list number of hits per corrected edit distance (zero or more)
#-----------------------------------------------------------------------------------------
hitsNM.smallRNA <- function(rname,NM.corr,smallRNA){
  NM <- table(NM.corr)
  mapped.tmp <- matrix(0,length(smallRNA)[1],length(NM))
  colnames(mapped.tmp) <- paste("NM",names(NM),sep="")
  rownames(mapped.tmp) <- sapply(strsplit(names(smallRNA)," "),"[[",1)
  for (j in names(NM)){
    mapped.tmp[names(table(rname[NM.corr==j])),paste("NM",j,sep="")] <- table(rname[NM.corr==j])
  }
  mapped.tmp
}

#==========================================================================================
# count number of aligned reads per isoform
#-----------------------------------------------------------------------------------------
count.isoforms <- function(rname,start.pos,aln,small){
  isoforms <- table(paste(rname,start.pos,aln))
  smallR <- sapply(strsplit(names(isoforms)," "),"[[",1)
  str.pos <- as.double(sapply(strsplit(names(isoforms)," "),"[[",2))
  aln.len <- as.double(sapply(strsplit(names(isoforms)," "),"[[",3))
  end.pos <- str.pos+aln.len-1
  mapped.smallRNA.parts <- data.frame(ID=smallR,part=small, aln.reads=as.vector(isoforms),start=str.pos,stop=end.pos,len=aln.len,stringsAsFactors=FALSE)
  mapped.smallRNA.parts
}


#==========================================================================================
# process mapping results for each run
#-----------------------------------------------------------------------------------------
getMappingData.smallRNA <- function(run.name, bam.file, percent.mismatches, unmatched5p, unmatched3p, small, smallRNA, fgp){
print("getMappingData")
  bam <- getBam(bam.file)
  bam <- addBamAln(bam,"global",percent.mismatches,unmatched5p,unmatched3p, fgp)
  
  mapped <- bam$qname[bam$mapped]
  bam.nm <- bam[bam$mapped,]
  rm(bam)
 
  #count number of aligned reads per isoform
  mapped.smallRNA.parts <- count.isoforms(bam.nm$rname,bam.nm$pos,bam.nm$aln,small)
   
  #create result matrix nr of mapped smallRNAs
  mapped.smallRNA <- matrix(0,length(smallRNA),1)
  rownames(mapped.smallRNA) <- sapply(strsplit(names(smallRNA)," "),"[[",1)
  colnames(mapped.smallRNA) <- run.name
   
  #count number of hits per smallRNA
  mapped.smallRNA[names(table(bam.nm$rname)),1] <- table(bam.nm$rname)

  #count number of hits per smallRNA and per edit distance
  if (fgp) mapped.smallRNA.NM <- hitsNM.smallRNA(bam.nm$rname,bam.nm$NM,smallRNA)
  else mapped.smallRNA.NM <- hitsNM.smallRNA(bam.nm$rname,bam.nm$NM.corr,smallRNA)
  
  #save results
  save(mapped.smallRNA,file=paste("obj_",run.name,"_mapped.",small,".out",sep=""))
  save(mapped.smallRNA.NM,file=paste("obj_",run.name,"_mapped.",small,".NM.out",sep="")) 
  save(mapped.smallRNA.parts, file=paste("obj_",run.name,"_mapped.",small,".part.out",sep=""))
  
  bam.nm
}

