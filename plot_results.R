#!/usr/bin/env Rscript
#### Plot alignment results
#### 2014-05-15 Inez Terpstra
#----------------------------------------------------------------------------------------------------
args <- commandArgs(TRUE)

base.dir <- args[1]
experiment <- args[2]
small <- args[3]
spikes <- args[4]

#----------------------------------------------------------
print(paste("Start",small,"visualization",sep=" "))
#----------------------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(ShortRead))
#----------------------------------------------------------------------------------------------------------
dir.res <- paste(base.dir,"/Results/",experiment,"/Output",sep="") 
dir.plot <- paste(base.dir,"/Images/",experiment,sep="") 

#---------------------------------------------------------------------------------------------------------
### read in alignment data
setwd(dir.res)

load(file=paste(dir.res,"/obj_",experiment,"_total.reads_",small,".out",sep=""))
load(file=paste(dir.res,"/obj_",experiment,"_mapped_",small,".out",sep=""))
load(file=paste(dir.res,"/obj_",experiment,"_read.lengths_",small,".out",sep=""))
load(file=paste(dir.res, "/obj_", experiment, "_aln.read.lengths_", small, ".out",sep=""))
# if spikes are mapped get spike alignment data else overwrite with small RNA aligned data
if (spikes=="y"){
  load(file=paste(dir.res,"/obj_",experiment,"_total.reads_spike",".out",sep=""))
  total.reads <- total.reads.all
  load(file=paste(dir.res,"/obj_",experiment,"_read.lengths_spike.out",sep=""))
  read.lengths <- read.lengths.all
} else {
  total.reads <- total.reads.all
  read.lengths <- read.lengths.all
}

# scale results: add 1 to all counts to prevent zero's before taking the log
# raw data
mapped <- mapped.all+0.5
# log data
mapped.log2 <- log2(mapped)

nr.runs <- dim(mapped)[2]
run.names <- colnames(mapped)

setwd(dir.plot)

#----------------------------------------------------------------------------------------
#hit densities
#----------------------------------------------------------------------------------------
titel <- paste("density",small,"hits",sep=" ")
plot.data <- mapped.log2
png(filename = paste(titel,".png",sep=""), width=635, height=460, units="px", bg="white")

if (dim(plot.data)[1] > 1){
    d.max <- 0
    for (i in 1:nr.runs){
      d <- density(plot.data[,i])
      if (max(d$y) > d.max) d.max <- max(d$y)
    }
    d <- density(plot.data[,1])
    plot(d,col=rainbow(15)[1],main=titel,ylim=c(0,round(d.max,2)))
    if (nr.runs>1){
      for (i in 2:nr.runs){
        d <- density(plot.data[,i])
        points(d,col=rainbow(15)[i],type="l")
      }
    }
} else if (length(plot.data) > 1){
  d <- density(plot.data)
  plot(d,main=titel)
}
dev.off()

#-------------------------------------------------------------------------------
#plot distributions
#------------------------------------------------------------------------------------------------------------------------------------------ 
#get axis limits 
max.counts <- matrix(NA,nr.runs,1)
max.lengths <- as.matrix(max.counts)
min.read.length <- 8
for (i in 1:nr.runs) {
  nrs <- read.lengths[[i]]
  tmp <- hist(nrs, breaks=c(min.read.length:max(nrs)), plot=FALSE)
  max.counts[i] <- max(tmp$counts)
  max.lengths[i] <- max(nrs)
} 
max.x <- max(max.lengths)
max.y <- max(max.counts)

xlabel <- "read length"
ylabel <- "freq"

for (hist.scale in c("ymax-scaled","x100.ymax-scaled")){
  titel.ext <- hist.scale
  if (substring(hist.scale,1,4)=="x100") breaks <- c(8:100) else breaks <- c(8:max.x)
  
  #plot
  for (i in 1:nr.runs){
    titel.base <- run.names[i]
    #create plotting device 
    png(filename = paste("Read length distr ",small," ",titel.ext," ",titel.base,".png",sep=""), width=29, height=21.7, units="cm", res=300, pointsize=12, bg="white")
      #set ylabels, xlabels depending on position of graph
      titel <- paste(small,titel.base,"-",prettyNum(total.reads[i,1],big.mark="."),sep=" ")
      nrs <- read.lengths[[i]]
      unm.nrs <- read.lengths.all[[i]]
      aln.nrs <- aln.read.lengths.all[[i]]
      plot.data <- nrs
      plot.data.unm <- unm.nrs
      plot.data.aln <- aln.nrs
      if (substring(hist.scale,1,4) == "x100"){
        plot.data <- nrs[nrs<=100]
        plot.data.unm <- unm.nrs[unm.nrs<=100]
        plot.data.aln <- aln.nrs[aln.nrs<=100]
      }
      if (hist.scale == "ymax-scaled") breaks <- c(8:max(nrs))
      ymax.lim <- ifelse (hist.scale == "x100-scaled", max(hist(plot.data, breaks=breaks, plot=FALSE)$counts), max.y)
      hist(plot.data, breaks=breaks, main=titel, xlab=xlabel, ylab=ylabel, col="lightgreen", border="black", ylim=c(0,ymax.lim), xlim=range(breaks))  
      hist(plot.data.unm, breaks=breaks, col="lightblue", add=TRUE)
      hist(plot.data.aln, breaks=breaks, col="blue", add=TRUE)
      legend.text <- c("all","unmapped",small)
      legend("topright", legend.text, fill=c("lightgreen","lightblue","blue"))
    #close current device
    dev.off()
  }
}

#-------------------------------------------------------------------------------
# plot nr of reads and hits per run
# plot percentage of hits per run
#------------------------------------------------------------------------------------------------------------------------------------------ 
if (small=="miRNA") small <- "hairpin"

for (logged.scale in c("y","")){
  titel <- paste("total number of reads and hits ", small, if(logged.scale=="y") "-Log scale",sep="")
  png(filename = paste(titel,".png",sep=""), width = 480, height = 480, units = "px", pointsize = 12,bg = "white")
    ymin <- min(total.reads.all[,small])/2
    ymax <- max(total.reads[,1])
    midpts <- barplot( t(total.reads[,1]), names.arg="", col="lightgreen", border="black", main=titel, ylab="Nr of reads", log=logged.scale, axes=FALSE, axisnames=FALSE, ylim=c(ymin,ymax))
    axis(2)
    axis(1, at = midpts, colnames(mapped), cex.axis = 0.8, las=2)
    barplot( t(total.reads.all[,"unmapped"]), col="lightblue", border="black", add=TRUE, log=logged.scale, ylim=c(ymin, ymax), axes=FALSE, axisnames=FALSE)
    barplot( t(total.reads.all[,small]), col="blue",  border="black", add=TRUE, log=logged.scale, ylim=c(ymin, ymax), axes=FALSE, axisnames=FALSE)
    legend.text <- c("all","unmapped",small)
    legend("topright",legend.text,fill=c("lightgreen","lightblue","blue"))
  dev.off()
}

titel <- paste("precentage of hits ",small,sep="")
png(filename = paste(titel,".png",sep=""), width = 480, height = 480, units = "px", pointsize = 12, bg = "white")
    ymin <- 0
    ymax <- max(100*total.reads.all[,small]/total.reads.all[,"unmapped"])+5
  barplot( 100*t(total.reads.all[,small]/total.reads.all[,"unmapped"]), col="lightblue", border="black", main=titel, ylab="% of reads", axes=FALSE, axisnames=FALSE, ylim=c(0, ymax))
    axis(2)
    axis(1, at = midpts, colnames(mapped), cex.axis = 0.8, las=2)
  barplot( 100*t(total.reads.all[,small]/total.reads[,1]), col="blue", border="black", ylim=c(ymin,ymax),add=TRUE, axes=FALSE, axisnames=FALSE)
  legend.text <- c("% of unmapped","% of all")
  legend("topright",legend.text,fill=c("lightblue","blue"))
dev.off()


