#!/usr/bin/env Rscript
#### Plot spike alignment results
#### 2014-05-15 Inez Terpstra
#----------------------------------------------------------------------------------------------------
print("Start spike visualization")
#----------------------------------------------------------------------------------------------------
args <- commandArgs(TRUE)

base.dir <- args[1]
experiment <- args[2]
dir.db <- args[3]
db <- args[4]

#----------------------------------------------------------
suppressPackageStartupMessages(library(ShortRead))
#----------------------------------------------------------------------------------------------------------
dir.fq <- paste(base.dir,"/Scratch/Fastq",sep="")
dir.res <- paste(base.dir,"/Results/",experiment,"/Output",sep="") 
dir.plot <- paste(base.dir,"/Images/",experiment,sep="") 

#----------------------------------------------------------------------------------------------------------------------------------------
# retrieve spike data
#-------------------------------------------------------------------------------------------------------------------------
spike.data <- readDNAStringSet(paste(dir.db,"spikes.fasta",sep="/"))
spikes <- names(spike.data)
NS.spikes <- spikes[grep("NSPK",spikes)]
SS.spikes <- spikes[grep("SSPK",spikes)]

#---------------------------------------------------------------------------------------------------------
### read in alignment data
setwd(dir.res)

molecule <- "spike"
load(file=paste(dir.res,"/obj_",experiment,"_mapped_",molecule,".out",sep=""))
load(file=paste(dir.res,"/obj_",experiment,"_total.reads_", molecule,".out",sep=""))

# scale results: add 1 to all counts to prevent zero's before taking the log
# raw data
mapped <- mapped.all+0.5
# log data
mapped.log2 <- log2(mapped)

#-------------------------------------------------------------------------------------
#normalization spikes
#-----------------------------------------------------------------------------------------
titel <- "normalization spikes"
plot.data <-  as.matrix(mapped[NS.spikes,])
nr.runs <- dim(plot.data)[2]
png(filename = paste(titel,".png",sep=""), width=635, height=460, units="px", bg="white")
  plot.data <- as.matrix(plot.data[order(apply(plot.data,1,mean)),])
  plot(apply(plot.data,1,mean), apply(plot.data,1,max), ylim=c(min(plot.data), max(plot.data)), pch=19, log="xy", main=titel, xlab="mean hits", ylab="hits", type="n")
  for (i in 1:nr.runs) {
    points(apply(plot.data,1,mean),plot.data[,i],pch=20,type="b",col=rainbow(nr.runs)[i])
  }
dev.off()

data_raw <- as.matrix(mapped[SS.spikes,])
titel <- "size spikes"
plot.data <-  data_raw
png(filename = paste(titel,".png",sep=""), width=635, height=460, units="px", bg="white")
  x.data <- as.double(sapply(strsplit(SS.spikes,"_"),"[[",3))
  plot(x.data, apply(plot.data,1,max)+1, ylim=c(min(plot.data)+1, max(plot.data)), pch=19, log="y", main=titel, xlab="nt", ylab="hits", type="n", axes=FALSE)
  axis(1,at=x.data,labels=x.data)
  axis(2)
  for (i in 1:nr.runs) {
    points(x.data,plot.data[,i],pch=20,type="b",col=rainbow(nr.runs)[i])
  }
dev.off()


data_norm <- data_raw/rowMeans(data_raw)

titel <- "size spikes relative to experiment mean"
plot.data <-  data_norm
png(filename = paste(titel,".png",sep=""), width=635, height=460, units="px", bg="white")
  plot(x.data, apply(plot.data,1,max)+1, ylim=c(min(plot.data), max(plot.data)), pch=19, main=titel, xlab="nt", ylab="hits", type="n", axes=FALSE, log="y")
  axis(1,at=x.data,labels=x.data)
  axis(2)
  for (i in 1:nr.runs) {
    if (max(plot.data[,i]) > 2 & min(plot.data[,i]) >= 0.5){
      lty <- 1
      lwd <- 2
    } else if (min(plot.data[,i]) < 0.5 & max(plot.data[,i]) <= 2) { 
      lty <- 4
      lwd <- 3
    } else if (max(plot.data[,i]) > 2 & min(plot.data[,i]) < 0.5) { 
      lty <- 3
      lwd <- 4
    } else {
      lty <- 3
      lwd <- 1
    }
    points(x.data,plot.data[,i],pch=20,type="l",lwd=lwd,lty=lty,col=rainbow(nr.runs)[i])
  }
  abline(h=1)
dev.off()


data_norm_intern <- t(t(data_raw)/colSums(data_raw))

titel <- "size spikes relative to run mean"
plot.data <-  data_norm_intern
png(filename = paste(titel,".png",sep=""),width=635,height=460,units="px",bg="white")
  plot(x.data,apply(plot.data,1,max)+1,ylim=c(min(plot.data),max(plot.data)),pch=19,log="y",main=titel,xlab="nt",ylab="hits",type="n",axes=FALSE)
  axis(1,at=x.data,labels=x.data)
  axis(2)
  for (i in 1:nr.runs) {
    points(x.data,plot.data[,i],pch=20,type="l",col=rainbow(nr.runs)[i])
  }
  abline(h=1)
dev.off()


data_norm_intern_double <- data_norm_intern/rowMeans(data_norm_intern)

titel <- "size spikes relative to run mean and experiment mean"
plot.data <-  data_norm_intern_double
png(filename = paste(titel,".png",sep=""),width=635,height=460,units="px",bg="white")
  plot(x.data,apply(plot.data,1,max)+1,ylim=c(min(plot.data),max(plot.data)),pch=19,main=titel,xlab="nt",ylab="hits",type="n",axes=FALSE,log="y")
  axis(1,at=x.data,labels=x.data)
  axis(2)
  for (i in 1:nr.runs) {
    if (max(plot.data[,i]) > 2 & min(plot.data[,i]) >= 0.5){
      lty <- 1
      lwd <- 2
    } else if (min(plot.data[,i]) < 0.5 & max(plot.data[,i]) <= 2) { 
      lty <- 4
      lwd <- 3
    } else if (max(plot.data[,i]) > 2 & min(plot.data[,i]) < 0.5) { 
      lty <- 3
      lwd <- 4
    } else {
      lty <- 3
      lwd <- 1
    }
    points(x.data,plot.data[,i],pch=20,type="l",lwd=lwd,lty=lty,col=rainbow(nr.runs)[i])
  }
  abline(h=1)
dev.off()

#----------------------------------------------------------------------------------------
#hit densities
#----------------------------------------------------------------------------------------
titel <- paste("density",molecule,"hits",sep=" ")
plot.data <- mapped.log2

png(filename = paste(titel,".png",sep=""), width=635, height=460, units="px", bg="white")
  d.max <- 0
  for (i in 1:nr.runs){
    d <- density(plot.data[,i])
    if (max(d$y) > d.max) d.max <- max(d$y)
  }
  d <- density(plot.data[,1])
  plot(d,col=rainbow(15)[1],main=titel,ylim=c(0,round(d.max,2)))
  for (i in 2:nr.runs){
    d <- density(plot.data[,i])
    points(d,col=rainbow(15)[i],type="l")
  }
dev.off()

#-------------------------------------------------------------------------------
# plot nr of reads and hits per run
# plot percentage of hits per run
#------------------------------------------------------------------------------------------------------------------------------------------ 
for (logged.scale in c("y","")){
  titel <- paste("total number of reads and hits ", molecule, if(logged.scale=="y") "-Log scale",sep="")
  png(filename = paste(titel,".png",sep=""), width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
    barplot( t(total.reads.all[,"total"]), names.arg=colnames(mapped), col="lightgreen", border="black", main=titel, xlab="Run", ylab="Nr of reads", log=logged.scale, axes=TRUE, axisnames=TRUE, beside=TRUE, ylim=c(min(total.reads.all[,molecule])/2,max(total.reads.all[,"total"])))
    barplot( t(total.reads.all[,molecule]), col="blue",  border="black", add=TRUE, beside=TRUE, log=logged.scale, ylim=c(min(total.reads.all[,molecule])/2, max(total.reads.all[,"total"])))
    legend.text <- c("all",molecule)
    legend("topright",legend.text,fill=c("lightgreen","blue"))
  dev.off()
}

titel <- paste("precentage of hits ",molecule,sep="")
png(filename = paste(titel,".png",sep=""), width = 480, height = 480, units = "px", pointsize = 12, bg = "white")
  barplot( 100*t(total.reads.all[,molecule]/total.reads.all[,"total"]), names.arg=colnames(mapped),  col="lightblue", border="black", main=titel, xlab="Run", ylab="% of reads", axes=TRUE, axisnames=TRUE, beside=TRUE, ylim=c(0, max(100*total.reads.all[,molecule]/total.reads.all[,"total"])))
  legend.text <- "% of all"
  legend("topright",legend.text,fill="lightblue")
dev.off()





