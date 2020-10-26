
#library(SingleR)
library(Seurat)
library(dplyr)
library(reshape2)
library(tidyr)
library(randomcoloR)
library(RColorBrewer)

setwd("/data/changc11/testtest/CT26_tumor_shBCL9_scRNAseq_results/")

#############load hsbcl9 and veh data#################

All_raw_counts<-read.table(file=paste0("All_Samples_Count.txt"),header=T)
rownames(All_raw_counts)<-make.unique(All_raw_counts[,2])
All_raw_counts<-All_raw_counts[,-(1:2)]
pbmc_col_HY_pbmc<- CreateSeuratObject(counts = All_raw_counts, min.cells = 3, min.features = 200, project = "col_HY")
col_HY_names<-data.frame(x=colnames(pbmc_col_HY_pbmc) ) %>% separate(x, c("Group","Sample"))
pbmc_col_HY_pbmc@meta.data[["Sample"]] <-col_HY_names$Sample
pbmc_col_HY_pbmc@meta.data[["Group"]] <-col_HY_names$Group

#############load NT and KD data#################

hsBCL9.pbmc.18 <- Read10X(data.dir = paste0("18"))
hsBCL9.pbmc.19 <- Read10X(data.dir = paste0("19"))
hsBCL9.pbmc.22 <- Read10X(data.dir = paste0("22"))
hsBCL9.pbmc.24 <- Read10X(data.dir = paste0("24"))
hsBCL9.pbmc.25 <- Read10X(data.dir = paste0("25"))
hsBCL9.pbmc.29 <- Read10X(data.dir = paste0("29"))

pbmc_18 <- CreateSeuratObject(counts = hsBCL9.pbmc.18, min.cells = 3, min.features = 200, project = "A18_KD")
pbmc_19 <- CreateSeuratObject(counts = hsBCL9.pbmc.19, min.cells = 3, min.features = 200, project = "B19_KD")
pbmc_22 <- CreateSeuratObject(counts = hsBCL9.pbmc.22, min.cells = 3, min.features = 200, project = "C22_KD")
pbmc_24 <- CreateSeuratObject(counts = hsBCL9.pbmc.24, min.cells = 3, min.features = 200, project = "D24_NT")
pbmc_25 <- CreateSeuratObject(counts = hsBCL9.pbmc.25, min.cells = 3, min.features = 200, project = "E25_NT")
pbmc_29 <- CreateSeuratObject(counts = hsBCL9.pbmc.29, min.cells = 3, min.features = 200, project = "F29_NT")

#############Merge data#################

pbmc_AB <- merge(pbmc_18, pbmc_19,add.cell.ids = c("A18_KD","B19_KD"),project = "hsBCL9",merge.data = TRUE)
pbmc_CD <- merge(pbmc_22, pbmc_24,add.cell.ids = c("C22_KD","D24_NT"),project = "hsBCL9",merge.data = TRUE)
pbmc_EF <- merge(pbmc_25, pbmc_29,add.cell.ids = c("E25_NT","F29_NT"),project = "hsBCL9",merge.data = TRUE)

pbmc_ABCD <- merge(pbmc_AB, pbmc_CD,project = "hsBCL9",merge.data = TRUE)
pbmc_ABCDEF <- merge(pbmc_ABCD, pbmc_EF,project = "hsBCL9",merge.data = TRUE)

sample_table<-data.frame(x=pbmc_ABCDEF$orig.ident ) %>% separate(x, c("Sample", "Group"))

pbmc_ABCDEF@meta.data[["Sample"]] <-sample_table$Sample
pbmc_ABCDEF@meta.data[["Group"]] <-sample_table$Group

################################################################
###############clustering of all cell population################
################################################################


pbmc_ABCDEF[["percent.mt"]] <- PercentageFeatureSet(pbmc_ABCDEF, pattern = "^mt-")
pbmc_ABCDEF <- NormalizeData(pbmc_ABCDEF, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_ABCDEF <- FindVariableFeatures(pbmc_ABCDEF, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_ABCDEF)
pbmc_ABCDEF <- ScaleData(pbmc_ABCDEF, features = rownames(pbmc_ABCDEF))

pbmc_col_HY_pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc_col_HY_pbmc, pattern = "^mt-")
pbmc_col_HY_pbmc <- NormalizeData(pbmc_col_HY_pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_col_HY_pbmc <- FindVariableFeatures(pbmc_col_HY_pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_col_HY_pbmc)

pbmc_col_HY_pbmc <- ScaleData(pbmc_col_HY_pbmc, features = rownames(pbmc_col_HY_pbmc))
#pbmc <- ScaleData(pbmc, vars.to.regress = c("nFeature_RNA","nCount_RNA", "percent.mt"), features = rownames(pbmc))

################################################################
#########Remove batch effects by FindIntegrationAnchors#########
################################################################

ob.list <- list(pbmc_col_HY_pbmc, pbmc_ABCDEF)
object.anchors <- FindIntegrationAnchors(object.list = ob.list, dims = 1:30)
pbmc <- IntegrateData(anchorset = object.anchors, dims = 1:30)


pbmc <- ScaleData(object = pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc)
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:5)
pbmc <- FindClusters(pbmc, resolution = 0.1)


################################################################
#########################Figure 1 A-C###########################
################################################################

palette_seurat_clusters <- randomColor(length(unique(pbmc$seurat_clusters)))
palette_Sample <- brewer.pal( length(unique(pbmc$Sample)),name ="Set3")
palette_Group <- brewer.pal( length(unique(pbmc$Group)),name ="Dark2")

pdf(paste0("Figure_1_C_hsBCL9_Seurat_individualgenes_pbmc_umap_cluster_",Sys.Date(),".pdf"))
print(DimPlot(pbmc, reduction = "umap",label = TRUE,cols = palette_seurat_clusters ))
dev.off()

pdf(paste0("hsBCL9_Seurat_individualgenes_pbmc_tsne_cluster_",Sys.Date(),".pdf"))
print(DimPlot(pbmc, reduction = "tsne",label = TRUE,cols = palette_seurat_clusters))
dev.off()

pdf(paste0("Figure_1_A_hsBCL9_Seurat_individualgenes_pbmc_umap_Sample_",Sys.Date(),".pdf"))
print(DimPlot(pbmc, reduction = "umap", group.by = "Sample",label = TRUE ,cols = palette_Sample))
dev.off()

pdf(paste0("Figure_1_B_hsBCL9_Seurat_individualgenes_pbmc_umap_Group_",Sys.Date(),".pdf"))
print(DimPlot(pbmc, reduction = "umap", group.by = "Group" ,cols = palette_Group))
dev.off()

pdf(paste0("hsBCL9_Seurat_individualgenes_pbmc_tsne_Sample_",Sys.Date(),".pdf"))
print(DimPlot(pbmc, reduction = "tsne", group.by = "Sample",label = TRUE ,cols = palette_Sample))
dev.off()

pdf(paste0("hsBCL9_Seurat_individualgenes_pbmc_tsne_Group_",Sys.Date(),".pdf"))
print(DimPlot(pbmc, reduction = "tsne", group.by = "Group" ,cols = palette_Group))

CAF_gene_list<-c("Mzb1","Ighg1","Ighg3","Fcrl5","Cd79a","Cd79b","Kdr","Clec4g","Ehd3","Plpp3","Ednrb","Vwf","Colec11","Col14a1","Col3a1","Cxcl12","Aebp1","Prss23","Dcn","Rarres2","Stab2","Csf1r","Cd3g","Ebf1","Irf8","Sox9","Apoc3","Top2A","Jchain","Dcn","Wnt9b","Rspo3","Cdh13","Wnt2","Ednrb","Jag1","Lrg1","Efnb1","Ltbp4","Adgrg6,Fcgr2b","Gpr182","Acta2","Aebp1","Ccdc80","Cfd","Col12a1","Col1a1","Col1a2","Col3a1","Col5a2","Col6a1","Cox4i2","Ctsk","Dcn","Dpt","Esam","Fstl1","Gja4","Gng11","Gsn","Higd1b","Lum","Mfap5","Mmp11","Myl9","Pdgfra","Postn","Rgs5","Tagln","Thbs2","Vcan","Stab2","Csf1r","Cd3g","Ebf1","Irf8","Sox9","Apoc3","Top2A","Jchain","Dcn","Ptrf","Sdpr","Murc","Prkcdbp","Ehd2","Cav1","Cav2","Cav3","Rarres2","Vwf","Clec4g","Wnt9b","Rspo3","Cdh13","Wnt2","Ednrb","Jag1","Lrg1","Efnb1","Ltbp4","Adgrg6","Fcgr2b","Gpr182","Kdr","Clec4g","Ehd3","Plpp3","Ednrb","Vwf","Colec11","Col14a1","Col3a1","Cxcl12","Aebp1","Prss23","Dcn","Rarres2","Stab2","Csf1r","Cd3g","Ebf1","Irf8","Sox9","Apoc3","Top2A","Jchain","Dcn","Wnt9b","Rspo3","Cdh13","Wnt2","Ednrb","Jag1","Lrg1","Efnb1","Ltbp4","Adgrg6","Fcgr2b","Gpr182","H2−Eb1","Wfdc2","Tmprss2","Krt8","Mlph","Fxyd3","Fam25c","Adig","Cd24a","Gipc2","Prtn3","Gata3","Krt19","Cidec","Chchd10","Cdo1","Areg","Epcam","Wfdc21","Hp","Steap4","Glul","Cd200","AW112010","Mrap","Tmem176a","Prlr","Clu","Cfd","Car3","Tmem176b","Retn","Wfdc18","Lcn2","Sncg","Krt18","Csn3","Krt7","Mgst1","Spint2","Rab25","Lpl","Apoc1","Cldn4","Ltc4s","Plin2","Rgcc","Trf","Stc2","Kdr","Clec4g","Ehd3","Plpp3","Ednrb","Vwf","Colec11","Col14a1","Col3a1","Cxcl12","Aebp1","Prss23","Dcn","Rarres2","Stab2","Csf1r","Cd3g","Ebf1","Irf8","Sox9","Apoc3","Top2A","Jchain","Dcn","Wnt9b","Rspo3","Cdh13","Wnt2","Ednrb","Jag1","Lrg1","Efnb1","Ltbp4","Adgrg6,Fcgr2b","Gpr182","Acta2","Aebp1","Ccdc80","Cfd","Col12a1","Col1a1","Col1a2","Col3a1","Col5a2","Col6a1","Cox4i2","Ctsk","Dcn","Dpt","Esam","Fstl1","Gja4","Gng11","Gsn","Higd1b","Lum","Mfap5","Mmp11","Myl9","Pdgfra","Postn","Rgs5","Tagln","Thbs2","Vcan")
endo_list<-intersect(rownames(pbmc),CAF_gene_list)
for (plotgene in endo_list)
  
{
  pdf(paste0("hsBCL9_Seurat_individualgenes_hsBCL9_pbmc_partgenes_umap_",plotgene,Sys.Date(),".pdf"))
  print(FeaturePlot(pbmc, features = plotgene,cols = c("grey","red") ))
  dev.off()
  
  pdf(paste0("hsBCL9_Seurat_individualgenes_hsBCL9_pbmc_partgenes_tsne_",plotgene,Sys.Date(),".pdf"))
  print(FeaturePlot(pbmc, reduction = "tsne",features = plotgene,cols = c("grey","red") ))
  dev.off()
  
}


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, num.cores = 64)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


################################################################
#########################Figure 1 D ############################
################################################################

pdf(paste0("Figure_1_D_Seurat_DoHeatmap_",Sys.Date(),".pdf"),10,40)
print(DoHeatmap(pbmc, features = top10$gene) + NoLegend())
dev.off()

################################################################
######################save endothelial cellss###################
################################################################

hsBCL9_pbmc_part_filtered_endo<- hsBCL9_pbmc_part_filtered_endo[,which(hsBCL9_pbmc_part_filtered_endo$seurat_clusters %in% c(7,8))]
saveRDS(hsBCL9_pbmc_part_filtered_endo, file = "hsBCL9_pbmc_part_filtered_endo.rds")
