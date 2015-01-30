require(minfi)
require(IlluminaHumanMethylation450kmanifest)
require(RColorBrewer)
require(foreach)
require(ggplot2)

source("cluster.R")
source("qnorm.R")


idat.folder <- "/well/zondervan/ENDOX/Illumina450K/idat" #no trailing /
pheno.file <- "/well/zondervan/ENDOX/Illumina450K/pheno.txt"
id.name <- "Array.ID"
group.id <- "Sample.ID"

##this tells minfi which idats to read, don't put any additional files into this folder
filenames <- paste(idat.folder,"/",unique(substr(dir(idat.folder),1,17)), sep="")

##read in phenos
phe <- read.table(pheno.file, header=T, sep="\t")
rownames(phe) <- phe[,id.name]

##read in idats, reorder to pheno and add names
RGset <- read.450k(filenames, verbose=TRUE)
RGset <- RGset[, rownames(phe)]

save(RGset, file="out/RGset.RData")

##EpiMigrant Pipeline - includes detection p-value (1E-16) and call rate QC: cr>98%
beta <- QNorm(RGset)

##minfi
##MSet.norm <- preprocessIllumina(RGset, bg.correct = TRUE, normalize = "controls", reference = 2)
##beta <- getBeta(MSet.norm)


colnames(beta) <- phe[colnames(beta), group.id]

## RSid Probes
snpProbesI <- getProbeInfo(RGset, type = "SnpI")
snpProbesII <- getProbeInfo(RGset, type = "SnpII")

IA <- rbind(getRed(RGset)[snpProbesI[snpProbesI$Color=="Red","AddressA"],],getGreen(RGset)[snpProbesI[snpProbesI$Color=="Grn","AddressA"],])
IB <- rbind(getRed(RGset)[snpProbesI[snpProbesI$Color=="Red","AddressB"],],getGreen(RGset)[snpProbesI[snpProbesI$Color=="Grn","AddressB"],])

IIA <- getRed(RGset)[snpProbesII[,"AddressA"],]
IIB <- getGreen(RGset)[snpProbesII[,"AddressA"],]

SNP.A <- rbind(IA, IIA)
SNP.B <- rbind(IB, IIB)

SNP <- SNP.A/(SNP.A + SNP.B)

png("snp.png",height=1000,width=1000)
par(mfrow=c(5,5))
for(i in 1:25){plot(IA[,i],IB[,i])}
dev.off()

SNP <- SNP[,rownames(phe) ]
colnames(SNP) <- phe[colnames(SNP), "Sample.ID"]

pdf("plots/snp.clust.pdf", width=18)
plot.cluster.colors(SNP, colors=as.numeric(phe$Subject))
##legend("topright",legend=levels(phe$Tissue5),col=brewer.pal(10, "Set3"), pch=15)
dev.off()

## missingness (per sample QC) mind>2% :
mind <- apply(beta,2,function(x){sum(is.na(x))/length(x)})
remove <- names(which(mind>0.02))
if(length(remove)>0){
  beta <- beta[,-which(colnames(beta) %in% remove)]
  phe <- phe[-which(rownames(phe) %in% remove),]
}

save(beta, phe, file="out/beta.RData")

detP <- detectionP(RGset)


p.mat <- foreach(i=-1:-16 ,.combine=c) %do% sum(detP<10^(i))/prod(dim(detP))



foreach(i=-1:-16 ,.combine=cbind) %do% {

  sample.cr <- apply(detP, 2, function(x){sum(x<10^(i))/nrow(detP)})
  c(sum(sample.cr<0.99),sum(sample.cr<.98),sum(sample.cr<.95),sum(sample.cr<0.90),sum(sample.cr<0.50))
}
