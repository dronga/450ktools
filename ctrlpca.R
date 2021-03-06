require(minfi)

CtrlProbePC <- function(RGset) {
  
  #BSC1 control probes
  BSCI.Green.Name =getProbeInfo(RGset, type = "Control")[16:18,]$ExtendedType
  BSCI.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(BSCI.Green.Name), dimnames = list(BSCI.Green.Name, sampleNames(RGset)))
  BSCI.Green[BSCI.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[16:18,]$Address,]
  BSCI.Red.Name =getProbeInfo(RGset, type = "Control")[22:24,]$ExtendedType
  BSCI.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(BSCI.Red.Name), dimnames = list(BSCI.Red.Name, sampleNames(RGset)))
  BSCI.Red[BSCI.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[22:24,]$Address,]
  
  #BSC2 control probes
  BSCII.Red.Name =getProbeInfo(RGset, type = "Control")[28:31,]$ExtendedType
  BSCII.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(BSCII.Red.Name), dimnames = list(BSCII.Red.Name, sampleNames(RGset)))
  BSCII.Red[BSCII.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[28:31,]$Address,]
  
  #STAINING
  stain.Red.Name =getProbeInfo(RGset, type = "Control")[2,]$ExtendedType
  stain.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(stain.Red.Name), dimnames = list(stain.Red.Name, sampleNames(RGset)))
  stain.Red[stain.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[2,]$Address,]
  stain.Green.Name =getProbeInfo(RGset, type = "Control")[4,]$ExtendedType
  stain.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(stain.Green.Name), dimnames = list(stain.Green.Name, sampleNames(RGset)))
  stain.Green[stain.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[4,]$Address,]
  
  #EXTENSION
  extensionA.Red.Name =getProbeInfo(RGset, type = "Control")[7,]$ExtendedType
  extensionA.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(extensionA.Red.Name), dimnames = list(extensionA.Red.Name, sampleNames(RGset)))
  extensionA.Red[extensionA.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[7,]$Address,]
  extensionT.Red.Name =getProbeInfo(RGset, type = "Control")[8,]$ExtendedType
  extensionT.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(extensionT.Red.Name), dimnames = list(extensionT.Red.Name, sampleNames(RGset)))
  extensionT.Red[extensionT.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[8,]$Address,]
  extensionC.Green.Name =getProbeInfo(RGset, type = "Control")[9,]$ExtendedType
  extensionC.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(extensionC.Green.Name), dimnames = list(extensionC.Green.Name, sampleNames(RGset)))
  extensionC.Green[extensionC.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[9,]$Address,]
  extensionG.Green.Name =getProbeInfo(RGset, type = "Control")[10,]$ExtendedType
  extensionG.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(extensionG.Green.Name), dimnames = list(extensionG.Green.Name, sampleNames(RGset)))
  extensionG.Green[extensionG.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[10,]$Address,]
  
  #HYBRIDISATION
  hybridH.Green.Name =getProbeInfo(RGset, type = "Control")[11,]$ExtendedType
  hybridH.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(hybridH.Green.Name), dimnames = list(hybridH.Green.Name, sampleNames(RGset)))
  hybridH.Green[hybridH.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[11,]$Address,]
  hybridM.Green.Name =getProbeInfo(RGset, type = "Control")[12,]$ExtendedType
  hybridM.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(hybridM.Green.Name), dimnames = list(hybridM.Green.Name, sampleNames(RGset)))
  hybridM.Green[hybridM.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[12,]$Address,]
  hybridL.Green.Name =getProbeInfo(RGset, type = "Control")[13,]$ExtendedType
  hybridL.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(hybridL.Green.Name), dimnames = list(hybridL.Green.Name, sampleNames(RGset)))
  hybridL.Green[hybridL.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[13,]$Address,]
  
  #TARGET REMOVAL
  target.Green.Name =getProbeInfo(RGset, type = "Control")[14:15,]$ExtendedType
  target.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(target.Green.Name), dimnames = list(target.Green.Name, sampleNames(RGset)))
  target.Green[target.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[14:15,]$Address,]
  
  #Specificity I
  specI.Green.Name =getProbeInfo(RGset, type = "Control")[32:34,]$ExtendedType
  specI.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(specI.Green.Name), dimnames = list(specI.Green.Name, sampleNames(RGset)))
  specI.Green[specI.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[32:34,]$Address,]
  specI.Red.Name =getProbeInfo(RGset, type = "Control")[38:40,]$ExtendedType
  specI.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(specI.Red.Name), dimnames = list(specI.Red.Name, sampleNames(RGset)))
  specI.Red[specI.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[38:40,]$Address,]
  
  #Specificity II
  specII.Red.Name =getProbeInfo(RGset, type = "Control")[44:46,]$ExtendedType
  specII.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(specII.Red.Name), dimnames = list(specII.Red.Name, sampleNames(RGset)))
  specII.Red[specII.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[44:46,]$Address,]
  
  #NON POLYMORPHIC
  np.Red.Name =getProbeInfo(RGset, type = "Control")[47:48,]$ExtendedType
  np.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(np.Red.Name), dimnames = list(np.Red.Name, sampleNames(RGset)))
  np.Red[np.Red.Name,] <- getRed(RGset)[getProbeInfo(RGset, type = "Control")[47:48,]$Address,]
  np.Green.Name =getProbeInfo(RGset, type = "Control")[49:50,]$ExtendedType
  np.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(np.Green.Name), dimnames = list(np.Green.Name, sampleNames(RGset)))
  np.Green[np.Green.Name,] <- getGreen(RGset)[getProbeInfo(RGset, type = "Control")[49:50,]$Address,]
  
  #Normalisation
  control=getProbeInfo(RGset, type = "Control")
  normC.Green.Name=control[control[,2]=='NORM_C',4]
  normC.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normC.Green.Name), dimnames = list(normC.Green.Name, sampleNames(RGset)))
  normC.Green[normC.Green.Name,] <- getGreen(RGset)[control[control[,2]=='NORM_C',1],]
  normG.Green.Name=control[control[,2]=='NORM_G',4]
  normG.Green <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normG.Green.Name), dimnames = list(normG.Green.Name, sampleNames(RGset)))
  normG.Green[normG.Green.Name,] <- getGreen(RGset)[control[control[,2]=='NORM_G',1],]
  normA.Red.Name=control[control[,2]=='NORM_A',4]
  normA.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normA.Red.Name), dimnames = list(normA.Red.Name, sampleNames(RGset)))
  normA.Red[normA.Red.Name,] <- getRed(RGset)[control[control[,2]=='NORM_A',1],]
  normT.Red.Name=control[control[,2]=='NORM_T',4]
  normT.Red <- matrix(NA_real_, ncol = ncol(RGset), nrow = length(normT.Red.Name), dimnames = list(normT.Red.Name, sampleNames(RGset)))
  normT.Red[normT.Red.Name,] <- getRed(RGset)[control[control[,2]=='NORM_T',1],]
  
  #combine ctrl probe intensities
  ctrl = cbind(as.matrix(colMeans(BSCI.Green)), as.matrix(colMeans(BSCI.Red)), as.matrix(colMeans(BSCII.Red)), t(stain.Red), t(stain.Green), t(extensionA.Red), t(extensionT.Red), t(extensionC.Green), t(extensionG.Green), t(hybridH.Green), t(hybridM.Green), t(hybridL.Green),as.matrix(colMeans(target.Green)),as.matrix(colMeans(specI.Green)),as.matrix(colMeans(specI.Red)), as.matrix(colMeans(specII.Red)),as.matrix(np.Red[1,]),as.matrix(np.Red[2,]),as.matrix(np.Green[1,]),as.matrix(np.Green[2,]),as.matrix(colMeans(normC.Green)),as.matrix(colMeans(normG.Green)), as.matrix(colMeans(normA.Red)),as.matrix(colMeans(normT.Red)))
  
  #perform PCA
  pca <- prcomp(ctrl)
  pc <- predict(pca)
  colnames(pc) <- paste(colnames(pc),'_cp', sep='')
  
  return(pc)
  
}
