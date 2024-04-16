rm(list = ls()) # clear the environment 
setwd("~/Work2020/Projects/20201109SingleCell")

library(Seurat)
library(org.Mm.eg.db)
library(clusterProfiler)
library(openxlsx)
library(ggplot2)
library(enrichplot)
library(dplyr)

macrophages <- readRDS("Macrophage_AT_filter50~6000mt20.rds")

projectName<- "Macrophage_AT_filter50~6000mt20"
for (i in 0:3) {
  for (j in (i+1):4) {
    name <- paste(i,"vs",j,sep = "")
    clusterMarkers <- FindMarkers(macrophages, ident.1 = i, ident.2 = j, min.pct = 0.25)
    path_res<-PathwayAnalyse(clusterMarker = clusterMarkers , paste(projectName,"/PathwayAnalyse",i,"vs",j,".xlsx", sep = ""),p=0.25)
    dotplot(path_res, showCategory = 10, title = paste("Enriched Pathways of",name) , split=".sign") + facet_grid(.~.sign)
    ggsave(paste(projectName,"/PathwayDotPlot",i,"vs",j,".pdf", sep = ""))
    
  }
}


#Compare Cluster 3 to Cluster 01245678
name <- "3vs01245678"
clusterMarkers <- FindMarkers(macrophages, ident.1 = 3, ident.2 = c(0,1,2,4,5,6,7,8), min.pct = 0.25)
path_res<-PathwayAnalyse(clusterMarker = clusterMarkers,fileName=paste(projectName,"/PathwayANalyse",name,".xlsx",sep = ""),p=1)
#Draw Kegg Pathway Dot Plot according to the list given by my boss
accPathway <- c("Cholesterol metabolism","Glycolysis / Gluconeogenesis","PPAR signaling pathway","Synaptic vesicle cycle","Lysosome","Ribosome","Oxidative phosphorylation","Carbon metabolism","Phagosome","Thermogenesis","Metabolic pathways")
accDot<-dotplot(path_res, showCategory = accPathway, title = "Activated in Cluster 3 Compared to Cluster 01245678", split=".sign") + facet_grid(.~.sign)

decPathway <- c("Endocytosis","Cytokine-cytokine receptor interaction","Chemokine signaling pathway","MAPK signaling pathway","C-type lectin receptor signaling pathway","Measles","Herpes simplex virus 1 infection","Toll-like receptor signaling pathway","AGE-RAGE signaling pathway in diabetic complications","Signaling pathways regulating pluripotency of stem cells","NF-kappa B signaling pathway","Th17 cell differentiation","Inflammatory bowel disease","Th1 and Th2 cell differentiation","TNF signaling pathway","IL-17 signaling pathway","Viral protein interaction with cytokine and cytokine receptor")
decDot<-dotplot(path_res, showCategory = decPathway, title = "Suppressed in Cluster 3 Compared to Cluster 01245678", split=".sign") + facet_grid(.~.sign)

accDot+decDot
ggsave(paste(projectName,"/PathwayDotPlot",name,".pdf", sep = ""),width = 15,height = 10)

#GSEA Plot
gseaplot2(path_res,geneSetID = c("mmu04979","mmu03320","mmu04142","mmu04721"))
ggsave(paste(projectName,"/GSEA3vs01245678ActivatedPathway.pdf", sep = ""),width =8, height = 5)

gseaplot2(path_res,geneSetID = c("mmu04668","mmu04064","mmu04657","mmu04061"), pvalue_table = TRUE)
ggsave(paste(projectName,"/GSEA3vs01245678SuppressedPathway.pdf", sep = ""),width = 15, height = 7)

#Compare Cluster 3 to Cluster 01245678
name <- "3vs01256"
clusterMarkers <- FindMarkers(macrophages, ident.1 = 3, ident.2 = c(0,1,2,5,6), min.pct = 0.25)
path_res<-PathwayAnalyse(clusterMarker = clusterMarkers,fileName=paste(projectName,"/PathwayANalyse",name,".xlsx",sep = ""),p=1)
#Draw Kegg Pathway Dot Plot according to the list given by my boss
accPathway <- c("Cholesterol metabolism","Glycolysis / Gluconeogenesis","PPAR signaling pathway","Synaptic vesicle cycle","Lysosome","Ribosome","Oxidative phosphorylation","Carbon metabolism","Phagosome","Thermogenesis","Metabolic pathways")
accDot<-dotplot(path_res, showCategory = accPathway, title = "Activated in Cluster 3 Compared to Cluster 01245678", split=".sign") + facet_grid(.~.sign)

decPathway <- c("Endocytosis","Cytokine-cytokine receptor interaction","Chemokine signaling pathway","MAPK signaling pathway","C-type lectin receptor signaling pathway","Measles","Herpes simplex virus 1 infection","Toll-like receptor signaling pathway","AGE-RAGE signaling pathway in diabetic complications","Signaling pathways regulating pluripotency of stem cells","NF-kappa B signaling pathway","Th17 cell differentiation","Inflammatory bowel disease","Th1 and Th2 cell differentiation","TNF signaling pathway","IL-17 signaling pathway","Viral protein interaction with cytokine and cytokine receptor")
decDot<-dotplot(path_res, showCategory = decPathway, title = "Suppressed in Cluster 3 Compared to Cluster 01245678", split=".sign") + facet_grid(.~.sign)

accDot+decDot
ggsave(paste(projectName,"/PathwayDotPlot",name,".pdf", sep = ""))

#GSEA Plot
path_res@result["mmu04979",]$Description<-"Cholesterol metabolism"
path_res@result["mmu03320",]$Description<-"PPAR signaling pathway"
path_res@result["mmu04142",]$Description<-"Lysosome"
path_res@result["mmu04721",]$Description<-"Synaptic vesicle cycle"
gseaplot2(path_res,geneSetID = c("mmu04979","mmu03320","mmu04142","mmu04721"), pvalue_table = TRUE)
ggsave(paste(projectName,"/GSEA3vs01256ActivatedPathway.pdf", sep = ""))

gseaplot2(path_res,geneSetID = c("mmu04668","mmu04064","mmu04657","mmu04061"), pvalue_table = TRUE)
ggsave(paste(projectName,"/GSEA3vs01256SuppressedPathway.pdf", sep = ""),width = 15, height = 7)



VlnPlot(macrophages,c("Nr1h3","Nfe2l1","Pparg"))

DimPlot(macrophages, reduction = "umap", label = TRUE,split.by = 'orig.ident')
ggsave(paste(projectName,"/UmapSplitBySample.pdf", sep = ""),width = 10, height = 3)

VlnPlot(macrophages,c("Atp1a1","Pim1","Pcolce2","Plat","Serpine1"))


VlnPlot(macrophages,c("Ffar4","Gpr120"))
VlnPlot(macrophages,c("Cxcl2","Il1b","Tnf","Ccl2", "Ccl3","Ccl4","Ccl6","Lyve1","Mrc1","Cd209d",
                      "Ccr2","Retnla"))
ggsave(paste(projectName,"/marker genes1.pdf", sep = ""))
VlnPlot(macrophages,c("Trem2","Gpnmb","Fabp5","Ctsd","Lpl","Mki67","Stmn1"
                      ,"Top2a","Mcm3","Mcm6","Ear2","Cd226","Fn1"))
ggsave(paste(projectName,"/marker genes2.pdf", sep = ""))


VlnPlot(macrophages,c("F13a1","Lgals3","Cd9","Atp6v0d2","Tubb5"))


VlnPlot(macrophages,c("Lyz1"))


genelist <- c("F13a1","Lgals3","Cd9","Atp6v0d2","Tubb5")
genelist <- c("Cd74")
genelist <- c("Cxcl2","Il1b","Tnf","Ccl2", "Ccl3","Ccl4","Ccl6","Lyve1","Mrc1","Cd209d",
              "Ccr2","Retnla","Trem2","Gpnmb","Fabp5","Ctsd","Lpl","Mki67","Stmn1"
              ,"Top2a","Mcm3","Mcm6","Ear2","Cd226","Fn1","Mcm5")
genelist <- c("Cd36","Adgre1","Fcgr1","Csf1r")
genelist <- c("F5","Fn1","Selp","Hp","Garnl3","Serpinb2","Vmn2r26","Hdc","Cd62p")
#resident macrophage list
genelist <- c("Pf4","Cx3cr1","Ccl8")
genelist <- c("Cd163","Mertk","Hmox1","Slc40a1","Tfrc","Lrp6","Abca1","Stab1","Mafb")
genelist <- c("Lipa","Lpl","Ctsl","Fabp4","Lgals1")
for( gene in genelist)
{
  VlnPlot(macrophages,gene)
  ggsave(paste("violinPlot/",gene,".pdf",sep = ""),width = 4, height = 3)
}

toAnalyseSerObj <- readRDS("Macrophage_AT_filter50~6000mt20.rds")
ggsave(paste(projectName,"/top10GeneFeatures.pdf", sep = ""),width = 15, height =10)

genelist<-c("Trem2","Gpnmb","Fabp5","Cd9","Lgals3","Atp6v0d2")
genelist<-c("Pcolce2","Lyve1","Cd63","Cd81","Fcgr4")

for( gene in genelist)
{
  FeaturePlot(macrophages, features = gene)
  ggsave(paste("FeaturePlot/",gene,".pdf",sep = ""),width = 6, height = 5)
}
FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
                                          "CCL2", "PPBP"), min.cutoff = "q9")

clusterMarkers <- FindMarkers(toAnalyseSerObj, ident.1 = c(4,7,8), ident.2 = c(0,1,2,3,5,6))
write.xlsx(clusterMarkers, paste(projectName,"/ClusterMarkers478vs012356.xlsx", sep = ""),rowNames = TRUE)
path_res<-PathwayAnalyse(clusterMarker = clusterMarkers , paste(projectName,"/PathwayAnalyse478vs012356.xlsx", sep = ""),p=0.25)

clusterMarkers$label = NA

labels <- c("Lgals3", "Ftl1", "Fabp4", "Ighm", "Scd1", "Cavin1", "Cavin2", "Cav1","Trem2","Lars2")
for (labelname in labels)
{
  clusterMarkers[labelname,"label"]<-labelname
}
p<-ggplot(clusterMarkers, aes(x=avg_log2FC, y=-log10(p_val_adj),label = label,color='red',lineheight = 1)) + geom_point() + ggtitle("Volcano") + xlab("log2 FC") + ylab("-log10 adjusted p-value")
+geom_text(check_overlap = TRUE,color= "black",vjust = 0)
p+ geom_text_repel(data          = clusterMarkers[labels,],
                   colour        = "black",
                   nudge_y       = 5,
                   size          = 4,
                   box.padding   = 0.5,
                   point.padding = 0.5,
                   force         = 100,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "x")

EnhancedVolcano(clusterMarkers,
                lab = rownames(clusterMarkers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                labhjust = 0.25,
                labvjust = 0.5
)

ggsave(paste(projectName,"/volcano/478vs012356connected.pdf", sep = ""))
gene<- "Itgax"
gene<-"S100a8"
gene<-"Tyrobp"
VlnPlot(macrophages,gene)
ggsave(paste("violinPlot/",gene,".pdf",sep = ""),width = 4, height = 3)
VlnPlot(macrophages,gene,group.by = "orig.ident")
ggsave(paste("violinPlot/",gene,"orig.pdf",sep = ""),width = 3, height = 4)
gene<-c("Itgax","Fcgr1")
VlnPlot(macrophages,gene)
VlnPlot(macrophages,gene,group.by = "orig.ident")

#VAM genes
VlnPlot(macrophages,c("Mef2a","Foxo1","Foxp1","Irf3","Klf4","Klf2"))
VlnPlot(macrophages,c("Cebpd","Etv5","Maf","Jun","Rora"))

#Cavity Macropahge
VlnPlot(macrophages,c("Cd226","Itgax","Ccr2","Retnla"))
