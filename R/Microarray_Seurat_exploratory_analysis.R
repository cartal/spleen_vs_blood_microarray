###Script for exploraroty analysis of RNA-Seq data with Seurat
###Not intended for full-scale analysis

##Load required packages

library("Seurat")
library("dplyr")
library("readr")
my_palette <- colorRampPalette(c("steelblue", "white", "red"))(n = 299)

###Read table

dataset <- read_delim("Blood_vs_Spleen_merged_datasets.csv", ",", escape_double = F, trim_ws = T)

###Conver gene identifiers to uppercase

dataset <- data.frame(lapply(dataset, function(v) {if (is.character(v)) return(toupper(v))else return(v)}))

###Make gene names unique

rownames(dataset) <- make.unique(as.character(dataset$genes), sep = '_')
dataset$genes <- NULL

##Create a Seurat object

microarray <- CreateSeuratObject(raw.data = dataset, min.cells = 2, min.genes = 500, do.scale = T)

##Normalise data

#microarray <- NormalizeData(microarray)

##Scale data - in this case I will try the 'negative binomial', but default is 'linear'. (Other options: 'poisson')

microarray <- ScaleData(microarray, model.use = 'linear', use.umi = F)

##Calculate gene disperssion and vairable genes across samples

microarray <- FindVariableGenes(microarray, num.bin = 10, do.plot = T)

##Perform linear dimensionality reduction a.k.a PCA clustering

microarray <- RunPCA(microarray, pc.genes = microarray@var.genes, do.print = T, genes.print = 5)

##Projet PCA analysis. This part scores each gene with its correlation with the calculated components. 

microarray <- ProjectPCA(microarray)

##Visualise PCA 

VizPCA(microarray, 1:4)

#tiff("/home/cartal/experiments/5-spleen_vs_blood_microarrays/Spleen_vs_Blood_merged_lumi_Seurat_Markers/Spleen_vs_Blood_normalised_PCA.tiff", width = 14, height = 6, units = 'in', res = 600)
PCAPlot(microarray, 1, 2, pt.size = 3.5, do.label = T, plot.title = 'PCA of gene expression values between Blood and Spleen')
#dev.off()

##Visualise componets in heatmaps

PCHeatmap(microarray, pc.use = 1:6, cells.use = 51, do.balanced = T, remove.key = T, label.columns = F, use.full = T, col.use = my_palette) #PCA panels

##JackStraw sampling procedure to determine statistically significant genes

microarray <- JackStraw(microarray, num.pc = 9, num.replicate = 5000, do.print = T)
JackStrawPlot(microarray, PCs = 1:9)

##Run t-SNE clustering for group identification

microarray <- RunTSNE(microarray, dims.use = c(1,2,3,6), do.fast = F, seed.use = 1712, add.iter = 90000)
TSNEPlot(microarray, pt.size = 2.5, do.label = T)

###Find population-specific markers 

##Student's t-test 

microarray.t_markers <- FindAllMarkers(microarray, only.pos = T, thresh.use = 0.35, return.thresh = 1e-5, test.use = 't', min.cells = 2, logfc.threshold = 1.5, random.seed = 1712, latent.vars = c('BD0', 'SD0'))
microarray.t_markers_sig <- subset(microarray.t_markers, p_val_adj < 1e-4)
t_top30 <- microarray.t_markers %>% group_by(cluster) %>% top_n(30, avg_logFC)
microarray@ident <- ordered(microarray@ident, levels = c("BD0", "BD2", "BD4", "BD6", "BD8", "BD10", "BD12", "SD0", "SD2", "SD4", "SD6", "SD8", "SD10", "SD12"))

tiff("/home/cartal/experiments/5-spleen_vs_blood_microarrays/Spleen_vs_Blood_merged_lumi_Seurat_Markers/Spleen_vs_Blood_normalised_180120_t_markers_heatmap.tiff", width = 14, height = 8, units = 'in', res = 600)
DoHeatmap(object = microarray, genes.use = t_top20$gene, slim.col.label = T, remove.key = T, col.low = 'steelblue', col.mid = 'white', col.high = 'red', cex.row = 4)
dev.off()
write.table(microarray.t_markers_sig, sep = ',', row.names = T, './Spleen_vs_Blood_merged_lumi_Seurat_Markers/Spleen_vs_Blood.all_markers_180120_linear.t_markers_sig.csv')

##AUC ROC classfier as implemented in Seurat 

microarray.roc_markers <- FindAllMarkers(microarray, only.pos = T, thresh.use = 0.35, return.thresh = 1e-5, test.use = 'roc', min.cells = 2, logfc.threshold = 1.5, random.seed = 1712, latent.vars = c('BD0', 'SD0'))
microarray.roc_markers_sig <- subset(microarray.roc_markers, power > 0.90)
roc_top20 <- microarray.roc_markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

tiff("/home/cartal/experiments/5-spleen_vs_blood_microarrays/Spleen_vs_Blood_merged_lumi_Seurat_Markers/Spleen_vs_Blood_normalised_180120_roc_markers_heatmap.tiff", width = 14, height = 8, units = 'in', res = 600)
DoHeatmap(object = microarray, genes.use = roc_top20$gene, slim.col.label = T, remove.key = T, col.low = 'steelblue', col.mid = 'white', col.high = 'red', cex.row = 4)
dev.off()
write.table(microarray.roc_markers_sig, sep = ',', row.names = T, './Spleen_vs_Blood_merged_lumi_Seurat_Markers/Spleen_vs_Blood.all_markers_180120_linear.roc_markers_sig.csv')

##Tobit-test for DE as implemented in Seurat 

microarray.tobit_markers <- FindAllMarkers(microarray, only.pos = T, thresh.use = 0.35, return.thresh = 1e-5, test.use = 'tobit', min.cells = 2, logfc.threshold = 1.25, random.seed = 1712, latent.vars = c('BD0', 'SD0'))
microarray.tobit_markers_sig <- subset(microarray.tobit_markers, p_val_adj < 1e-4)
tobit_top20 <- microarray.tobit_markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

tiff("/home/cartal/experiments/5-spleen_vs_blood_microarrays/Spleen_vs_Blood_merged_lumi_Seurat_Markers/Spleen_vs_Blood_normalised_180120_tobit_markers_heatmap.tiff", width = 14, height = 8, units = 'in', res = 600)
DoHeatmap(object = microarray, genes.use = tobit_top20$gene, slim.col.label = T, remove.key = T, col.low = 'steelblue', col.mid = 'white', col.high = 'red', cex.row = 4)
dev.off()
write.table(microarray.roc_markers_sig, sep = ',', row.names = T, './Spleen_vs_Blood_merged_lumi_Seurat_Markers/Spleen_vs_Blood.all_markers_180120_poisson.tobit_markers_sig.csv')

###Find condition-specific markers using the ROC AUC classifier 

SD2.roc_markers <- FindMarkers(microarray, ident.1 = c('SD2'), ident.2 = c('SD0', "SD4", "SD6", "SD8"), thresh.use = 0.35, test.use = 'roc', min.cells = 2, logfc.threshold = 0.75, random.seed = 1786)
SD2.roc_markers_sig <- subset(SD2.roc_markers, power > 0.90)
print(head(SD2.roc_markers_sig, 20))
DoHeatmap(object = microarray, genes.use = rownames(SD2.roc_markers_sig), slim.col.label = T, remove.key = T, col.low = 'steelblue', col.mid = 'white', col.high = 'red', cex.row = 5)
write.table(SD2.roc_markers_sig, sep = ',', row.names = T, 'SD2_blood_specific_merged_180118.roc_markers_sig.csv')



BD2.roc_markers <- FindMarkers(microarray, ident.1 = c('BD2'), ident.2 = c('BD0', "BD4", "BD6", "BD8", "BD10", "BD12"), thresh.use = 0.35, test.use = 'roc', min.cells = 2, logfc.threshold = 0.75, random.seed = 1786)
BD2.roc_markers_sig <- subset(BD2.roc_markers, power > 0.90)
print(head(BD2.roc_markers_sig, 20))
DoHeatmap(object = microarray, genes.use = rownames(BD2.roc_markers_sig), slim.col.label = T, remove.key = T, col.low = 'steelblue', col.mid = 'white', col.high = 'red', cex.row = 5)
write.table(BD2.roc_markers_sig, sep = ',', row.names = T, 'BD2_blood_specific_merged_180118.roc_markers_sig.csv')
  

