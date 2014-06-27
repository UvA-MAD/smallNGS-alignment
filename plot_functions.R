#### Functions for generation of Box and MA plots
#### 2014-05-28 Inez Terpstra
#### MA.mad: to generate MA plots of count data
#### Box.mad: to generate Box plots of count data
#### 2014-06-04 Added PCA plots
#### PCA.mad: to generate PCA plots based on count tables and design file 
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
 
#-------------------------------------------------------------------------------
#  Generate PCA
#-------------------------------------------------------------------------------
PCA.mad <- function(small, plot.data, normtype, design){

  ## remove all zero
  rm <- which(rowSums(plot.data==0)==dim(plot.data)[2])
  if (length(rm) > 0) plot.data <- plot.data[-rm,]
  data <- log2(plot.data+1)

  #color samples based on the values of each condition
  nr.conditions <- dim(design[,-1])[2]
  pr <- prcomp(data)
  if (sum(pr$sdev)>0){
    titel <- paste(small,normtype,sep="_")
    png(filename = paste(titel,".png",sep=""),width=480,height=480,units="px",bg="white")
      barplot(pr$sdev^2/sum(pr$sdev^2),ylab="proportion of variance",main=titel)
    dev.off()

    for (i in 1:nr.conditions){
      k <- i+1
      condition <- colnames(design)[k]
      factors <-  factor(design[,k])
      legend.text <- levels(factors)
      if (typeof(design[,k])=="double") kleur <- rainbow(nlevels(factors),start=0,end=1/6) else kleur <- rainbow(nlevels(factors)) 
      levels(factors) <- kleur
      kleuren <- as.character(factors)
      titel <- paste(small,normtype,condition,sep="_")
      png(filename = paste(titel,".png",sep=""),width=480,height=480,units="px",bg="white")
        par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
        plot(pr$r[,1],pr$r[,2],pch=19,col=kleuren,lwd=2,main=titel,xlab="PC1",ylab="PC2",bty="L") 
        textxy(pr$r[,1],pr$r[,2],design[,1],cex=0.7)
        legend("topright",inset=c(-0.2,0),legend.text, pch=19, col=kleur,xpd=TRUE)
      dev.off()
    }
  }
}  


