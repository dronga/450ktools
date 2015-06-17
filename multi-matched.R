require(RColorBrewer)
require(foreach)

sciNotation <- function(x, digits = 1) {
      if (length(x) > 1) {
                return(append(sciNotation(x[1]), sciNotation(x[-1])))
              }
          if (!x) return(0)
          exponent <- floor(log10(x))
          ##base <- sprintf("%.1f",x / 10^exponent)
      base <- formatC( round( x / 10^exponent, digits ), format='f', digits=digits )
      as.expression(substitute(base %*% 10^exponent,
              list(base = base, exponent = exponent)))
    }

require(doMC)

registerDoMC(cores=8)

##make priceanno
load("../anno.RData")

hits <- as.character(read.table("../bmi.hits.187.txt")[,1])

##load("bmi_matched_cpgs_160615.RData")
##matched.logic <- apply(match[,2:10001], 2, function(x){anno$IlmnID %in% x})
##save(matched.logic, file="matched.logic.RData")

load("matched.logic.RData")

## merge annotations
require(GenomicRanges)

files <- dir("IN")

hmark <- "DNase.macs2.narrowPeak"
hmark <- "H3K4me1.narrowPeak"
hmark <- "H3K4me3.narrowPeak"
hmark <- "H3K27ac.narrowPeak"

files <- files[grep(hmark, files)]


anno.granges <- GRanges(Rle(anno$CHR), IRanges(anno$MAPINFO, anno$MAPINFO), strand="*", Rle(anno$IlmnID))


peak.granges <- foreach(file=files, .combine=rbind) %dopar% {

print(file)
dataf <- read.table(paste("IN/",file,sep=""), as.is=T)

  dataf[,1] <- gsub("chr","",dataf[,1])
dataf[dataf[,1]==23, 1] <- "X"
dataf[dataf[,1]==24, 1] <- "Y"

data.granges <- GRanges(Rle(dataf[, 1]), IRanges(dataf[, 2], dataf[, 3]), strand="*")

ov <- as.matrix(findOverlaps(anno.granges,data.granges))
x <-1:nrow(anno)%in% ov[,1]
names(x) <- rownames(anno)

##p.vals <- foreach(x=features, .combine=c) %do% {
  ##x <- factor(x)
  ##p.vals <- foreach(y=levels(x), .combine=c) %do% {
  ##  p <- fisher.test(factor(x), anno$IlmnID %in% hits)$p.value

  obs <- sum(x & (anno$IlmnID %in% hits))

  #system.time(null <- apply(match,2,function(y){sum(x & (anno$IlmnID %in% y))}))

  system.time(null <- colSums(x&matched.logic))
  
  p.vals <- sum(null>=obs)/10001
  
  ##p.vals <- signif(p,3)
    #names(p.vals) <- levels(x)
  ##m1<- sum(x)/nrow(anno)
  ##m2<- sum(x[hits])/length(hits)
 c(p.vals,obs,mean(null))

}

peak.granges[peak.granges[,1]==0,1] <- 1e-4 

rownames(peak.granges) <- substr(files,0,4)

save(hmark,files,peak.granges, file=paste("results-",hmark,".RData",sep=""))


samples <- read.table("roadmap_samples.txt",as.is=T, sep="\t",comment.char="")
rownames(samples) <- samples[,1]

samples.id <- intersect(samples[,1], rownames(peak.granges))

peak.granges <- peak.granges[samples.id,]
samples <- samples[samples.id,]

legend <- samples[!duplicated(samples[,2:3]),2:3]

## pdf(paste("~/temp/",hmark,"-p-matched.pdf",sep=""), width=20, height=10)
## layout(rbind(1,2), heights=c(7,1))
## par(mar=c(20, 4, 4, 2))
## plot(1:nrow(peak.granges),-log10(peak.granges[,1]), pch=19, col=samples[,3],xaxt="n", xlab="", ylab="-logp", main=hmark)
## #barplot(-log(peak.granges[,1]), width=1, col=samples[,3],xaxt="n", xlab="", ylab="-logp", main=hmark)
## grid()
## axis(side=1,at=1:nrow(peak.granges),labels=F, las=2)
## mtext(text=samples[,4], side=1,at=1:nrow(samples),col=samples[,3],las=2,line=1)
## abline(h=0.05/nrow(peak.granges))
## par(mar=c(0,0,0,0))
## legend("topleft", legend=legend[,1], col=legend[,2], ncol=3, pch=19, cex=0.8)
## dev.off()

bonf <- 0.05 / nrow(peak.granges)

pdf(paste("",hmark,"-fold-matched_160515.pdf",sep=""), width=20, height=10)
layout(rbind(1,2), heights=c(7,1))
par(mar=c(20, 4, 4, 2))
y <- (peak.granges[,2]/peak.granges[,3])
plot(1:nrow(peak.granges),y, col=samples[,3],xaxt="n", xlab="", ylab="Observed/Expected", main=hmark, cex.axis=1.2, cex.lab=1.2,ylim=c(min(y),max(y)+1), pch=c(1,19)[as.numeric(peak.granges[,1]<bonf)+1])
#barplot(-log(peak.granges[,1]), width=1, col=samples[,3],xaxt="n", xlab="", ylab="-logp", main=hmark)
grid(nx=NA,ny=NULL)
axis(side=1,at=1:nrow(peak.granges),labels=F, las=2)
mtext(text=samples[,4], side=1,at=1:nrow(samples),col=samples[,3],las=2,line=1)
abline(h=bonf)
vlines <- cumsum(table(factor(samples[,2],levels=unique(samples[,2])))) + .5
abline(v=vlines,lty=2, col="grey")
text(x=vlines+0.0,y=max(y)+1,labels=legend[,1], srt=90, col=legend[,2], pos=2)
par(mar=c(0,0,0,0))
#legend("topleft", legend=legend[,1], col=legend[,2], ncol=3, pch=19, cex=0.8)
dev.off()

pdf(paste("",hmark,"-p-matched_160515.pdf",sep=""), width=20, height=10)
layout(rbind(1,2), heights=c(7,1))
par(mar=c(20, 4, 4, 2))
y <- -log10(peak.granges[,1])
plot(1:nrow(peak.granges),y, col=samples[,3],xaxt="n", xlab="", ylab="-log10(p)", main=hmark, cex.axis=1.2, cex.lab=1.2,ylim=c(min(y),max(y)+2), pch=c(1,19)[as.numeric(peak.granges[,1]<bonf)+1])
#barplot(-log(peak.granges[,1]), width=1, col=samples[,3],xaxt="n", xlab="", ylab="-logp", main=hmark)
grid(nx=NA,ny=NULL)
axis(side=1,at=1:nrow(peak.granges),labels=F, las=2)
mtext(text=samples[,4], side=1,at=1:nrow(samples),col=samples[,3],las=2,line=1)
abline(h=1)
vlines <- cumsum(table(factor(samples[,2],levels=unique(samples[,2])))) + .5
abline(v=vlines,lty=2,col="grey")
text(x=vlines+0.0,y=max(y)+2,labels=legend[,1], srt=90, col=legend[,2], pos=2)
par(mar=c(0,0,0,0))
#legend("topleft", legend=legend[,1], col=legend[,2], ncol=3, pch=19, cex=0.8)
dev.off()
