#### Functions for normalization and generations of Box and MA plots
#### 2014-05-28 Inez Terpstra
#### sizeFactors.mad: to calculate scaling factors of a dataset - DESeq method
#### MA.mad: to generate MA plots of count data
#### Box.mad: to generate Box plots of count data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  Calculate DESeq size factors (Martijs)
#-------------------------------------------------------------------------------
sizeFactors.mad <- function (counts, locfunc = median){
    loggeomeans <- rowMeans(log(counts))
    apply(counts, 2, function(cnts) exp(locfunc((log(cnts) -
        loggeomeans)[is.finite(loggeomeans)])))
}

#-------------------------------------------------------------------------------
#  Generate MA plot (Martijs)
#-------------------------------------------------------------------------------
MA.mad <- function (small,plot.data,normtype){
  ## remove all zero
  rm <- which(rowSums(plot.data==0)==dim(plot.data)[2])
  if (length(rm) > 0) plot.data <- plot.data[-rm,]
  runs <- colnames(plot.data) 
  spike.data <- plot.data[grep("NSPK",rownames(plot.data)),]
  ref <- apply(plot.data,1,median)
  ref <- log2(ref+1)
  logdat <- log2(plot.data+1)
  ids <- rownames(plot.data)
  for(i in 1:ncol(logdat)){
    M <- (logdat[,i] - ref)
    A <- (logdat[,i] + ref)/2
    titel <- paste(small,"_",runs[i],"_",normtype,sep="")
    png(paste(titel,".png",sep=""),width=384,heigh=384)
      plot(A,M,main=titel)
      points(A[grep("NSPK",ids)],M[grep("NSPK",ids)],pch=20,col="red")
      points(A[grep("SSPK",ids)],M[grep("SSPK",ids)],pch=20,col="green")
      abline(h=0)
    graphics.off()
  }
}

#-------------------------------------------------------------------------------
#  Generate boxplot
#-------------------------------------------------------------------------------
Box.mad <- function (small,plot.data,normtype,op){
  ## remove all zero
  rm <- which(rowSums(plot.data==0)==dim(plot.data)[2])
  if (length(rm) > 0) plot.data <- plot.data[-rm,]
  data <- log2(plot.data+1)
  titel <- paste(small,normtype,sep="_")    
  png(paste(titel,".png",sep=""),width=480,heigh=480)
    nam <- colnames(data)
    boxplot(data,names=nam,ylab="log2(counts)",main=titel,outpch=op,las=2,cex.axis=0.8)
  dev.off()
}  
 


