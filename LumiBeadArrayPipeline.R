###Load the Lumi package for Illumina BeadArray analyses

library("lumi", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("limma", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("beadarray", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")

###Read Illumina BeadArray data

x.lumi <- 'CB-blood-HJ_JingWen_7chip__011113_Allsamples_Bkg Subtracted_No Norm_Sample_Probe_Profile.txt'
IFNARKO <- lumiR.batch(x.lumi, lib.mapping = 'lumiMouseIDMapping', convertNuID = T)

###Generate QC plots

plot(IFNARKO, what = 'density')
plotCDF(IFNARKO, reverse = T)
plot(IFNARKO, what = 'pair')
pairs(IFNARKO, smoothScatter = T)
MAplot(IFNARKO, smoothScatter = T)
plot(IFNARKO, what = 'boxplot')
plot(IFNARKO, what = 'sampleRelation')
plotSampleRelation(IFNARKO, method = 'mds')
plotSampleRelation(IFNARKO, method = 'cluster')

##Do default VST variance stabilizing transform

IFNARKO.T <- lumiT(IFNARKO, ifPlot = T, method = 'log2')

##Plot VST transformation

trans <- plotVST(IFNARKO.T)
matplot(log2(trans$untransformed), trans$transformed)

##Do quantile between microarray normaliazation
IFNARKO.TN <- lumiN(IFNARKO.T)

##Do quality control estimation after normalization
IFNARKO.TNQ <- lumiQ(IFNARKO.TN)

###Generate QC plots after normalisation

plot(IFNARKO.TNQ, what = 'density')
plotCDF(IFNARKO.TNQ, reverse = T)
plot(IFNARKO.TNQ, what = 'pair')
pairs(IFNARKO.TNQ, smoothScatter = T)
MAplot(IFNARKO.TNQ, smoothScatter = T)
plot(IFNARKO.TNQ, what = 'sampleRelation')
plotSampleRelation(IFNARKO.TNQ, method = 'mds')
plotSampleRelation(IFNARKO.TNQ, method = 'cluster')
plot(IFNARKO.TNQ, what = 'boxplot')

###Write expression results matrix

write.exprs(IFNARKO.TNQ, file = 'IFNARKO_lumiTNQ_matrix.txt')

###Differential expression analysis

##Remove unexpressed genes

presentCount <- detectionCall(IFNARKO)
selDataMatrix <- IFNARKO.TNQ[presentCount > 0,]
probeList <- rownames(selDataMatrix)

###Define sample groups

sampleType <- c('Naive', 'Infected', 'Infected', 'Naive', 'Infected', 'Infected', 'Naive', 'Infected', 'Naive', 'Infected', 'Naive')

###Create comparative design

design <- model.matrix(~ 0 + factor(sampleType))
colnames(design) <- c('BD0', 'SD0')
design

###Fit model

IFNARKOfit <- lmFit(BvS, design)
IFNARKOfit <- eBayes(IFNARKOfit)

##Make contrasts

contrast.matrix <- makeContrasts(BD0 - SD0, levels = design)

##Re-fit

IFNARKOfit2 <- contrasts.fit(IFNARKOfit, contrast.matrix)

IFNARKOfit2 <- eBayes(IFNARKOfit2)


topTable(IFNARKOfit2, adjust.method = 'fdr', number = 40, lfc = 0.3, sort.by = 'p')












