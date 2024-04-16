generateUmap <- function (toRunSeurat,cellLables,projectName,lowFC = 50, highFC=3500, mtPer =20) 
{

  if(is.null(cellLables))
  {
    toAnalyseSerObj <-  toRunSeurat[[1]]
  }
  else
  {
    com_data <- merge(x = toRunSeurat[[1]], y = toRunSeurat[2:length(toRunSeurat)], 
                      add.cell.ids = cellLables,project =projectName)
    
    toAnalyseSerObj <- com_data
  }

  
  if(!dir.exists(projectName))
  {
    dir.create(projectName)
  }
  
  #Check genes.tsv and find mitochondrial related Genes are begin with lower case mt-
  toAnalyseSerObj[["percent.mt"]] <- PercentageFeatureSet(toAnalyseSerObj, pattern = "^mt-")
  
  # Visualize QC metrics as a violin plot
  VlnPlot(toAnalyseSerObj, features = c("nFeature_RNA","nCount_RNA","percent.mt"))
  ggsave(paste(projectName,"/1FeaturePlot.pdf",sep = ""))
  
  #FeatureScatter used to visualize feature to feature relationship
  #count_RNA has a nearly 0 Pearson Correlation with mt percent indicates not linearly related
  #count_RNA has a nearly 1 Pearson Correaltion with feature_rna indicates positively linearly related
  plot1 <- FeatureScatter(toAnalyseSerObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(toAnalyseSerObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  ggsave(paste(projectName,"/2ScatterPlot.pdf",sep = ""))
  
  #Filter,2576 cells to 2156 cells
  print(table(toAnalyseSerObj$orig.ident))
  #toAnalyseSerObj<- subset(toAnalyseSerObj, subset = nFeature_RNA > 50 & nFeature_RNA < 3500 & percent.mt<50)
  toAnalyseSerObj <- subset(toAnalyseSerObj, subset = nFeature_RNA > lowFC & nFeature_RNA < highFC & percent.mt<mtPer)
  
  #For Raw Data
  #toAnalyseSerObj <- subset(toAnalyseSerObj, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt<10)
  
  print(table(toAnalyseSerObj$orig.ident))
  
  
  #Normalize
  toAnalyseSerObj <- NormalizeData(toAnalyseSerObj)
  VlnPlot(toAnalyseSerObj, features = c("Hbb-bs","Hba-a1","Hba-a2","Hbb-bt"))
  ggsave(paste(projectName,"/2hbbExpression.pdf",sep = ""))
  #Filter "Hbb-bs","Hba-a1","Hba-a2","Hbb-bt">4
  toAnalyseSerObj <- subset(toAnalyseSerObj, subset = `Hbb-bs`< 4 &`Hba-a1`< 4 & `Hba-a2`< 4& `Hbb-bt`< 4)
  print(table(toAnalyseSerObj$orig.ident))
  
  #find 2000 variable count
  toAnalyseSerObj <- FindVariableFeatures(toAnalyseSerObj, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(toAnalyseSerObj), 10)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(toAnalyseSerObj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  # If error ocurrs,increse the size of viewport panel
  CombinePlots(plots = list(plot1,plot2))
  ggsave(paste(projectName,"/3VariableFeaturePlot.pdf",sep = ""), height = 5, width = 12)
  
  #Scaling the data for PCA
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
  toAnalyseSerObj <- FindNeighbors(toAnalyseSerObj, dims = 1:10)
  toAnalyseSerObj <- FindClusters(toAnalyseSerObj, resolution = 0.5)
  
  #Different dimension of PCA in use may form different Umap Even the cluster would be the same 
  toAnalyseSerObj <- RunUMAP(toAnalyseSerObj, dims = 1:10)
  DimPlot(toAnalyseSerObj, reduction = "umap", label = TRUE)
  ggsave(paste(projectName,"/Umap_clusters.pdf",sep = ""))
  
  DimPlot(toAnalyseSerObj, reduction = "umap", group.by = "orig.ident",cols = c("red","green","blue","grey","black","purple","brown","orange"))
  ggsave(paste(projectName,"/Umap_samples.pdf",sep = ""))
  saveRDS(toAnalyseSerObj, file = paste(projectName,".rds",sep = ""))
  
  #Write Excel about top20Markers each cluster
  markers <- FindAllMarkers(toAnalyseSerObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top20<-markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  openxlsx::write.xlsx(top20,paste(projectName,"/top20markers.xlsx", sep = ""))
  
markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  DoHeatmap(toAnalyseSerObj, features = top10$gene) 
  ggsave(paste(projectName,"/top10GeneFeatures.pdf", sep = ""), width = 40, height =30)
  
  #Cell Distribution Excel
  cellDistri<-table(toAnalyseSerObj@meta.data$seurat_clusters,toAnalyseSerObj@meta.data$orig.ident)
  openxlsx::write.xlsx(cellDistri,paste(projectName,"/cellDistribution.xlsx", sep = ""))
  
  
  # Maker Dir for markers
  markerDir<-paste(projectName,"/markers",sep = "")
  if(!dir.exists(markerDir))
  {
    dir.create(markerDir)
  }
  
  #Bcell Markers
  VlnPlot(toAnalyseSerObj, features = c("Ms4a1", "Cd79a","Mzb1","Cd79b","Ly6d"),pt.size = 0.05) 
  ggsave(paste(projectName,"/markers/BCellMarker.pdf", sep = ""))
  #TCell Markers
  VlnPlot(toAnalyseSerObj, features = c("Cxcr6","Icos","Cd3g","Il7r","Cd4","Cd8b1","Cd8a","Nkg7","Rag1","Lck","Cd3d"),pt.size = 0.05) 
  ggsave(paste(projectName,"/markers/TCellMarker.pdf", sep = ""))
  #Macrophage Markers
  VlnPlot(toAnalyseSerObj, features = c("Cxcl2","Cd14","Tnf","Adgre1","Fcgr1","Csf1r","Ccr2","Ly6c2","Il1b"),pt.size = 0.05)
  ggsave(paste(projectName,"/markers/MacrophageMarker.pdf", sep = ""))
  #Resident Macrophage
  VlnPlot(toAnalyseSerObj, features = c("Adgre1","Csf1r","Fcgr1","Cd68","F13a1","Lyve1","Gas6"),pt.size = 0.05) 
  ggsave(paste(projectName,"/markers/ResidentMacrophageMarker.pdf", sep = ""))
  #Master Cells Marker
  VlnPlot(toAnalyseSerObj, features = c("Furin","Il1rl1","Calca","Cd3d"),ncol = 2,pt.size = 0.05) 
  ggsave(paste(projectName,"/markers/MasterCellMarker.pdf", sep = ""))
  #Granulocytes Marker
  VlnPlot(toAnalyseSerObj, features = c("S100a8","S100a9","Ngp","Camp"), ncol = 2,pt.size = 0.05)
  ggsave(paste(projectName,"/markers/GranulocytesCellMarker.pdf", sep = ""))
  #NK cells
  #VlnPlot(toAnalyseSerObj, features = c("Klrb1c","Ncr1","Klra8","Klrc1","Klrb1"), ncol = 2)
  VlnPlot(toAnalyseSerObj,features = c("Klrb1c","Klra8","Klrc1","Cd3","Itgam","Cd27","Tnfrsf7","Il2rb","Cd161","Nkg2d","Nkp46","Klrk1","Ncr1"),ncol = 3,pt.size = 0.05)
  ggsave(paste(projectName,"/markers/NKCellMarker.pdf", sep = ""))
  #MoDCs/dendritic cells
  VlnPlot(toAnalyseSerObj, features = c("Cd209a","Cd74","Flt3","H2-Eb1"), ncol = 2,pt.size = 0.05) 
  ggsave(paste(projectName,"/markers/DCMarker.pdf", sep = ""))
  #Monocyte Marker
  VlnPlot(toAnalyseSerObj, features = c("Ly6c2","Ccr2","Csf1r","Adgre1","Itgam","Lilrb4a","Ly6c","Ly6c1","Nr4a1","Cx3cr1"), ncol = 3,pt.size = 0.05)
  ggsave(paste(projectName,"/markers/MonocyteMarker.pdf", sep = ""))
  
  
  VlnPlot(toAnalyseSerObj,features = c("Plat","Plau","Serpine1","Dach1","Lrp1","Atf6"),ncol = 3)
  ggsave(paste(projectName,"/markers/tPA-Related.pdf", sep = ""), width = 20, height = 10)
 
  VlnPlot(toAnalyseSerObj,features = c("Fech","Fth1","Cd163","Tfrc","Hmox1","Gpx4","Acsl4","Aifm2","Slc7a11","Txn1","Slc40a1"),ncol = 3,pt.size = 0.1)
  ggsave(paste(projectName,"/markers/ironMetabolismVlnPlot.pdf", sep = ""), width = 30, height = 20)
  
  VlnPlot(toAnalyseSerObj,features = c("Mrc1","Cd163","Retnla","Il10","Chil3","Tgfb1"),ncol = 3)
  ggsave(paste(projectName,"/markers/M2MarkersVlnPlot.pdf", sep = ""), width = 20, height = 10)
  
  VlnPlot(toAnalyseSerObj,features = c("Cx3cr1","Rac1","Mertk","Gas6","Lrp1","Stab1","Mfge8","S1pr1","Pros1"),ncol = 3)
  ggsave(paste(projectName,"/markers/EfferocytosisMarker.pdf", sep = ""), width = 25, height = 15)
  
  #VlnPlot(toAnalyseSerObj, features = c("Pcolce2","Scarb1","Atp1a1","Cd36","Pparg"), ncol = 3)
  VlnPlot(toAnalyseSerObj, features = c("Pcolce2","Scarb1","Atp1a1","Pparg","Cd36","Trem2","Pim1"), ncol = 3)
  ggsave(paste(projectName,"/markers/keyGeneVlnPlot.pdf", sep = ""), width = 20, height = 10)
  
  VlnPlot(toAnalyseSerObj,c("Acly"))
  ggsave(paste(projectName,"/markers/Acly.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Ccl2"))
  ggsave(paste(projectName,"/markers/Ccl2.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Fabp4","Fabp5"))
  ggsave(paste(projectName,"/markers/Fabp.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Cpt1a","Cpt1b","Cpt1c","Cpt2"))
  ggsave(paste(projectName,"/markers/Cpt.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Alox5","Alox12","Alox15"))
  ggsave(paste(projectName,"/markers/Alox.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Nrf1","Nrf2","Tfam","Ppargc1a","Tfb1m","Tfb2m","Surf1","Shmt2"))
  ggsave(paste(projectName,"/markers/MitochondriaSynthesi.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Pink1","Prkn","Atg3","Atg5","Atg7","Becn1"),pt.size = 0.1)
  ggsave(paste(projectName,"/markers/Mitophagy.pdf", sep = ""))
  
  # generateAortaList(c("Ccl8","Folr2","Pf4","F13a1","Cbr2","Lyve1","Ltc4s","C4b","Wfdc17"),"ResidentMac")
  # generateAortaList(c("Ccl2","Ccl3","Ccl4","Cxcl1","Cxcl2","Ccr12","Tnf","Il1b","Cd14"),"InflaMac")
  # generateAortaList(c("Trem2","Spp1","Gpnmb","Lgals3","Mmp12","Cd9","Ctsd","Ctsl","Lpl","Fabp5"),"Trem2FoamMac")
  # generateAortaList(c("Isg15","Irf7","Rsad2","Ifit3","Ifit1","Mnda","Ms4a4c","Ly6a","Ccl8","Ccl12"),"IfnicMac")
  # 
  VlnPlot(toAnalyseSerObj,c("Ccl2","Ccl3","Ccl4","Cxcl1","Cxcl2","Ccr12","Tnf","Il1b","Cd14"))
  ggsave(paste(projectName,"/markers/InflaMac.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Ccl8","Folr2","Pf4","F13a1","Cbr2","Lyve1","Ltc4s","C4b","Wfdc17"))
  ggsave(paste(projectName,"/markers/ResidentMac.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Trem2","Spp1","Gpnmb","Lgals3","Mmp12","Cd9","Ctsd","Ctsl","Lpl","Fabp5"))
  ggsave(paste(projectName,"/markers/Trem2FoamMac.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Isg15","Irf7","Rsad2","Ifit3","Ifit1","Mnda","Ms4a4c","Ly6a","Ccl8","Ccl12"))
  ggsave(paste(projectName,"/markers/IfnicMac.pdf", sep = ""))
  
  #Pcolce2, Plat, Cd36
  # FeaturePlot(toAnalyseSerObj, features = c("Pcolce2"))
  # ggsave(paste(projectName,"/markers/Pcolce2Distri.pdf", sep = ""))
  # 
  # FeaturePlot(toAnalyseSerObj, features = c("Plat"))
  # ggsave(paste(projectName,"/markers/PlatDistri.pdf", sep = ""))
  # 
  # FeaturePlot(toAnalyseSerObj, features = c("Cd36"))
  # ggsave(paste(projectName,"/markers/Cd36Distri.pdf", sep = ""))
  
  #SingleR Notation
  ref <- ImmGenData()
  pro<-SingleR(toAnalyseSerObj@assays$RNA@data,ref,labels = ref$label.fine,method= "cluster",clusters = toAnalyseSerObj@meta.data$seurat_clusters)
  openxlsx::write.xlsx(pro,paste(projectName,"/singleR.xlsx", sep = ""),rowNames = TRUE, colNames = TRUE)
  openxlsx::write.xlsx(pro$labels,paste(projectName,"/singleRJustLabel.xlsx", sep = ""),rowNames = TRUE, colNames = TRUE)

  VlnPlot(toAnalyseSerObj,c("Ccl8","Folr2","Pf4","F13a1","Cbr2","Lyve1","Ltc4s","C4b","Wfdc17"))
  ggsave(paste(projectName,"/markers/ResidentMac.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Ccl2","Ccl3","Ccl4","Cxcl1","Cxcl2","Ccr12","Tnf","Il1b","Cd14"))
  ggsave(paste(projectName,"/markers/InflaMac.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Trem2","Spp1","Gpnmb","Lgals3","Mmp12","Cd9","Ctsd","Ctsl","Lpl","Fabp5"))
  ggsave(paste(projectName,"/markers/Trem2FoamMac.pdf", sep = ""))
  
  VlnPlot(toAnalyseSerObj,c("Isg15","Irf7","Rsad2","Ifit3","Ifit1","Mnda","Ms4a4c","Ly6a","Ccl8","Ccl12"))
  ggsave(paste(projectName,"/markers/IfnicMac.pdf", sep = ""))
  
  
  clist <- c("Ear2","Fn1","Retnla","Clec4b1","Ccl6","Cripl","Ccr2","Lyz1")
  listName <- "ClassicalMonoDCgene"
  VlnPlot(toAnalyseSerObj,clist)
  ggsave(paste(projectName,"/markers/",listName,".pdf", sep = ""))
  
  clist <- c("Cd209a","Ifitm1","Napsa","Ifi30","Gm2")
  listName <- "MixedMonoDC"
  VlnPlot(toAnalyseSerObj,clist)
  ggsave(paste(projectName,"/markers/",listName,".pdf", sep = ""))
  
  clist <- c("Cd8a","Ccr9","Ldhb","Endou","Cd3d","Cd3g","Cd28","Cd8b1","Ly6d")
  listName <- "TCELL"
  VlnPlot(toAnalyseSerObj,clist)
  ggsave(paste(projectName,"/markers/",listName,".pdf", sep = ""))
  
  clist <- c("Cd79a","Mzb1","Ms4a1","Ebf1","Cd79b","Ccr7","Plac8","4930523C07Rik","Fcmr")
  listName <- "BCELL"
  VlnPlot(toAnalyseSerObj,clist)
  ggsave(paste(projectName,"/markers/",listName,".pdf", sep = ""))
}

library(Hmisc)

toCapitalize <- function(vec)
{
  vec<- tolower(vec)
  vec<- capitalize(vec)
  return(vec)
}
clist<-c("ACLY","CS","ACO1","ACO2","IDH1","IDH2","IDH3A","IDH3B","IDH3G","OGDH","DLST","DLD","SUCLG1","SUCLG2","SUCLA2","SDHA","SDHB","SDHC","SDHD","FH1","MDH1","MDH2")
library(data.table)

generateCSV <- function (toCSVSeurat,fileName) 
{
  
  toCSVSeurat[["percent.mt"]] <- PercentageFeatureSet(toCSVSeurat, pattern = "^mt-")
  toCSV<- subset(toCSV, subset = nFeature_RNA > 50 & nFeature_RNA < 3500 & percent.mt<50)
  data_to_write_out <- as.data.frame(as.matrix(toCSV@assays$RNA@data))
  fwrite(x = data_to_write_out, file = fileName,row.names = TRUE)
  
}


generateAortaList <- function (clist,listName) 
{
  
  VlnPlot(ApoE_Chow_Aorta,clist)
  ggsave(paste("tosent/ApoE_Chow_Aorta/",listName,".pdf", sep = ""))
  
  VlnPlot(ApoE_HFD_Aorta,clist)
  ggsave(paste("tosent/ApoE_HFD_Aorta/",listName,".pdf", sep = ""))
  
  VlnPlot(DKO_Chow_Aorta,clist)
  ggsave(paste("tosent/DKO_Chow_Aorta/",listName,".pdf", sep = ""))
  
  VlnPlot(DKO_HFD_Aorta,clist)
  ggsave(paste("tosent/DKO_HFD_Aorta/",listName,".pdf", sep = ""))
  
}

DrawMonocle <- function(toAnalyseSerObj)
{
  cds <- as.cell_data_set(toAnalyseSerObj)
  cds <- cluster_cells(cds)
  p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
  p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
  wrap_plots(p1, p2)
  
  integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
  cds <- as.cell_data_set(integrated.sub)
  cds <- learn_graph(cds)
  plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = TRUE, label_branch_points = TRUE)
  ggsave("p12Partition1Trajectory.pdf")
  
  #max.avp <- which(unlist(FetchData(integrated.sub, "Avp"))!=0)
  max.avp<- which(unlist(FetchData(integrated.sub, "Avp")) == max(unlist(FetchData(integrated.sub, "Avp"))))
  fmax.avp <- colnames(integrated.sub)[max.avp]
  cds <- order_cells(cds, root_cells = fmax.avp)
  #without specify the root_cells pops the graphic selection
  #cds <- order_cells(cds)
  plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
             label_branch_points = FALSE)
  ggsave("p12Partition1Pseudotime2.pdf")
}

PathwayAnalyse<- function(clusterMarker, fileName = "tosent/temp.xlsx",p = 0.05)
{
  #Extract the log fold changes for GSEA
  cluster_lfcs <- clusterMarker$avg_log2FC
  #names(cluster_lfcs) <- rownames(clusterMarker)
  names(cluster_lfcs) <- mapIds(org.Mm.eg.db, keys = rownames(clusterMarker), keytype = "SYMBOL", column = "ENTREZID")
  #names(cluster_lfcs) <- mapIds(org.Mm.eg.db, keys = names(cluster_lfcs), keytype = "SYMBOL", column = "ENTREZID")
  
  #Remove genes without ENTREZ IDs
  cluster_lfcs <- cluster_lfcs[!is.na(names(cluster_lfcs))]
  cluster_lfcs <- cluster_lfcs[!is.infinite(cluster_lfcs)]
  cluster_lfcs <- cluster_lfcs[!is.na(cluster_lfcs)]
  #Reorder the genes into descending order
  cluster_lfcs <- cluster_lfcs[order(cluster_lfcs, decreasing = TRUE)]
  
  
  #I've set the pvalue cutoff to 1 so we see all results (including non-significant)
  cluster2_vs_0_kegg <- gseKEGG(geneList = cluster_lfcs, organism = "mmu", pvalueCutoff = p)
  
  cluster2_vs_0_kegg <- setReadable(cluster2_vs_0_kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  
  View(cluster2_vs_0_kegg@result)
  write.xlsx(cluster2_vs_0_kegg@result,fileName)
  
  return(cluster2_vs_0_kegg)
  # dotplot(cluster2_vs_0_kegg, showCategory=30)
  # ggsave("goEnrichment.pdf")
}

