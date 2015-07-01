require(limma)

`%notin%` <- function(x,y) !(x %in% y)


QNorm <- function(RGset, detP = 0.01, MinD = 0.02, CallR = 0) {

##get detection p-values:

dp <<- detectionP(RGset, type = "m+u")

###### our pipeline – QN in 6 categories

## Type II probes
TypeII.Name <- getProbeInfo(RGset, type = "II")$Name
TypeII.Green <- getGreen(RGset)[getProbeInfo(RGset, type = "II")$Address,]
TypeII.Red <- getRed(RGset)[getProbeInfo(RGset, type = "II")$Address,]
rownames(TypeII.Red) <- TypeII.Name
colnames(TypeII.Red) <- sampleNames(RGset)
rownames(TypeII.Green) <- TypeII.Name
colnames(TypeII.Green) <- sampleNames(RGset)



## Type I probes, split into green and red channels
TypeI.Green.Name <- getProbeInfo(RGset, type = "I-Green")$Name
TypeI.Green.M <- getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressB,]
rownames(TypeI.Green.M) <- TypeI.Green.Name
colnames(TypeI.Green.M) <- sampleNames(RGset)
TypeI.Green.U <- getGreen(RGset)[getProbeInfo(RGset, type = "I-Green")$AddressA,]
rownames(TypeI.Green.U) <- TypeI.Green.Name
colnames(TypeI.Green.U) <- sampleNames(RGset)

TypeI.Red.Name <- getProbeInfo(RGset, type = "I-Red")$Name
TypeI.Red.M <- getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressB,]
rownames(TypeI.Red.M) <- TypeI.Red.Name
colnames(TypeI.Red.M) <- sampleNames(RGset)
TypeI.Red.U <- getRed(RGset)[getProbeInfo(RGset, type = "I-Red")$AddressA,]
rownames(TypeI.Red.U) <- TypeI.Red.Name
colnames(TypeI.Red.U) <- sampleNames(RGset)

##remove high missingness samples
d = ifelse(dp<detP,1,NA)
cr = data.frame(rowSums(is.na(d))/length(d[1,]))
exclude.badcalls = rownames(cbind(cr,rownames(cr))[cbind(cr,rownames(cr))[,1]> CallR,])

exclude.sites = exclude.badcalls
##exclude.sites = unique(rbind(as.matrix(exclude.chrX), as.matrix(exclude.chrY),as.matrix(exclude.cas),as.matrix(exclude.snps),as.matrix(crossmap),as.matrix(exclude.badcalls), as.matrix(exclude.mhc)))

mind = data.frame(colSums(is.na(d))/length(d[,1]))
remove.mind = rownames(cbind(mind,rownames(mind))[cbind(mind,rownames(mind))[,1]> MinD,])
samples = rownames(cbind(mind,rownames(mind))[cbind(mind,rownames(mind))[,1]<MinD,])

TypeII.Green =subset(TypeII.Green, select=samples)
TypeII.Red = subset(TypeII.Red, select=samples)
TypeI.Green.M = subset(TypeI.Green.M, select=samples)
TypeI.Green.U = subset(TypeI.Green.U, select=samples)
TypeI.Red.M = subset(TypeI.Red.M, select=samples)
TypeI.Red.U = subset(TypeI.Red.U, select=samples)

##set NAs
d = subset(dp, select = samples)
TypeII.Green = TypeII.Green * ifelse(d[rownames(TypeII.Green),]==0,1,NA)
TypeII.Red = TypeII.Red * ifelse(d[rownames(TypeII.Red),]==0,1,NA)
TypeI.Green.M = TypeI.Green.M * ifelse(d[rownames(TypeI.Green.M),]==0,1,NA)
TypeI.Green.U = TypeI.Green.U * ifelse(d[rownames(TypeI.Green.U),]==0,1,NA)
TypeI.Red.M = TypeI.Red.M * ifelse(d[rownames(TypeI.Red.M),]==0,1,NA)
TypeI.Red.U = TypeI.Red.U * ifelse(d[rownames(TypeI.Red.U),]==0,1,NA)

#--------------------------------------------------------------------------------------------------------------------------------
##calculate betas - no QN
## TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
## TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
## TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
## beta.noQN=rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)
## colnames(beta.noQN)<-gsub("^X","",colnames(beta.noQN))
## beta.noQN=as.matrix(beta.noQN)
##save(beta.noQN, file="beta_noQN.RData")
##rm(beta.noQN, TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)

##QN

#exclude sites
TypeII.Green = TypeII.Green[rownames(TypeII.Green) %notin% as.matrix(exclude.sites),]
TypeII.Red = TypeII.Red[rownames(TypeII.Red) %notin% as.matrix(exclude.sites),]
TypeI.Green.M = TypeI.Green.M[rownames(TypeI.Green.M) %notin% as.matrix(exclude.sites),]
TypeI.Green.U = TypeI.Green.U[rownames(TypeI.Green.U) %notin% as.matrix(exclude.sites),]
TypeI.Red.M = TypeI.Red.M[rownames(TypeI.Red.M) %notin% as.matrix(exclude.sites),]
TypeI.Red.U = TypeI.Red.U[rownames(TypeI.Red.U) %notin% as.matrix(exclude.sites),]


TypeII.Red.norm=normalizeQuantiles(TypeII.Red)
TypeII.Green.norm=normalizeQuantiles(TypeII.Green)
TypeI.Green.M.norm=normalizeQuantiles(TypeI.Green.M)
TypeI.Green.U.norm=normalizeQuantiles(TypeI.Green.U)
TypeI.Red.M.norm=normalizeQuantiles(TypeI.Red.M)
TypeI.Red.U.norm=normalizeQuantiles(TypeI.Red.U)
rm(TypeII.Red,TypeII.Green,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U)

#calculate betas
TypeII.betas = TypeII.Green.norm/(TypeII.Red.norm+TypeII.Green.norm+100)
TypeI.Green.betas = TypeI.Green.M.norm/(TypeI.Green.M.norm+TypeI.Green.U.norm+100)
TypeI.Red.betas = TypeI.Red.M.norm/(TypeI.Red.M.norm+TypeI.Red.U.norm+100)
beta=rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)
colnames(beta)<-gsub("^X","",colnames(beta))
beta <- as.matrix(beta)
return(beta)
##save(beta, beta_excluded, file="beta_QN2.RData")

}
