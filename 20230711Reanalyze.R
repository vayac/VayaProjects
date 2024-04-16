rm(list = ls()) # clear the environment 
setwd("~/Work2020/Projects/20201109SingleCell")

library(Seurat)
library(org.Mm.eg.db)
library(clusterProfiler)
library(openxlsx)
library(ggplot2)
library(enrichplot)
library(dplyr)
library(SingleR)
library(patchwork)

projectName <- "20230711Reanalyze"
toAnalyseSerObj<-readRDS("~/Work2020/Projects/20201109SingleCell/AT_filter50~6000mt20.rds")

colors <- c("red","green","aquamarine2","yellow","tan","purple","magenta","lightblue","khaki2","brown","wheat","lightsalmon","rosybrown4","deepskyblue4","darkolivegreen","gold")

#SingleR Notation
ref <- ImmGenData()
pro<-SingleR(toAnalyseSerObj@assays$RNA@data,ref,labels = ref$label.main)
toAnalyseSerObj$label <- pro$labels
DimPlot(toAnalyseSerObj, reduction = "umap", group.by = "label",cols = colors) + ggtitle("") 
ggsave(paste(projectName,"/AdiposeTissueUmap_final.pdf",sep = ""),width = 6,height = 5)
DimPlot(toAnalyseSerObj, reduction = "umap", split.by = "orig.ident",group.by = "label",ncol = 2,cols = colors)
ggsave(paste(projectName,"/AdiposeTissueUmap_origIdent.pdf",sep = ""),width = 11, height = 11)
 openxlsx::write.xlsx(table(predicted=toAnalyseSerObj$label,toAnalyseSerObj$orig.ident),paste(projectName,"/SingleRLabelNumAmongSamples.xlsx", sep = ""),rowNames = TRUE, colNames = TRUE)

macrophageGenes <- c("Mrc1","Cxcl2","Adgre1","Csf1r")
for (gene in macrophageGenes) {
  
  VlnPlot(toAnalyseSerObj,gene,group.by = "label",cols = c("red","green","aquamarine2","yellow","tan","purple","magenta","lightblue","khaki2","brown","wheat","lightsalmon","rosybrown4","deepskyblue4","darkolivegreen","gold"))
  ggsave(paste(projectName,"/",gene,"_Vln.pdf",sep = ""),height = 6,width = 12)
  FeaturePlot(toAnalyseSerObj,gene)
  ggsave(paste(projectName,"/",gene,"_FeaPlot.pdf",sep = ""),height = 8, width =7.8)
}

VlnPlot(toAnalyseSerObj,c("Mrc1","Cxcl2","Adgre1","Csf1r"),group.by = "lable", ncol = 2)
ggsave(paste(projectName,"/MacrophageMarker.pdf",sep = ""))
FeaturePlot(toAnalyseSerObj,c("Mrc1","Cxcl2","Adgre1","Csf1r"))
ggsave(paste(projectName,"/UmapMacrophageMarker.pdf",sep = ""))

#get subset Macrophage
#toAnalyseSerObj <- subset(toAnalyseSerObj, subset = lable =="Macrophages")

Macrophage <- readRDS("~/Work2020/Projects/20201109SingleCell/Macrophage_AT_filter50~6000mt20.rds")
DimPlot(Macrophage)
#macrophage <- subset(Macrophage, idents = c(0:8) )
macrophage <- subset(Macrophage, subset = Gata6>1 & Fn1 > 3 & F5 >1 & Selp>1 )
DimPlot(Macrophage, reduction = "umap",label = TRUE,pt.size = 0.01)#, split.by = "seurat_clusters",ncol =3) + ggtitle("") 
ggsave(paste(projectName,"/MacrophageUmap.pdf",sep = ""), width = 3, height=2.7)
DimPlot(Macrophage, reduction = "umap",split.by = "orig.ident",ncol = 2,pt.size = 0.01)
ggsave(paste(projectName,"/MacrophageSampleUmap.pdf",sep = ""), width = 3.7, height=4)






toAnalyseSerObj<- macrophage
toAnalyseSerObj <- NormalizeData(toAnalyseSerObj)
all.genes <- rownames(toAnalyseSerObj)
toAnalyseSerObj <- ScaleData(toAnalyseSerObj, features = all.genes)
#Perform PCA
toAnalyseSerObj <- RunPCA(toAnalyseSerObj, features = VariableFeatures(toAnalyseSerObj))
print(toAnalyseSerObj[["pca"]], dims = 1:5, nfeatures = 5)

toAnalyseSerObj <- JackStraw(toAnalyseSerObj, num.replicate = 100)
toAnalyseSerObj <- ScoreJackStraw(toAnalyseSerObj,dim = 1:20)
JackStrawPlot(toAnalyseSerObj, dims = 1:15)
ggsave(paste(projectName,"/JackStrawPlot.pdf",sep = ""))
ElbowPlot(toAnalyseSerObj)
ggsave(paste(projectName,"/ElbowPlot.pdf",sep = ""))

#Using different PCA dimensions to find neighbors and form cluster
toAnalyseSerObj <- FindNeighbors(toAnalyseSerObj, dims = 1:15)
toAnalyseSerObj <- FindClusters(toAnalyseSerObj, resolution = 0.25)

#Different dimension of PCA in use may form different Umap Even the cluster would be the same 
toAnalyseSerObj <- RunUMAP(toAnalyseSerObj, dims = 1:15)
DimPlot(toAnalyseSerObj, reduction = "umap", label = TRUE)
ggsave(paste(projectName,"/Umap_clusters.pdf",sep = ""))

FeaturePlot(toAnalyseSerObj,"Trem2")

VlnPlot(toAnalyseSerObj,c("Ccl3","Ccl4","Tnf"))
VlnPlot(toAnalyseSerObj,c("Lyve1","Cd209d","F13a1"))
