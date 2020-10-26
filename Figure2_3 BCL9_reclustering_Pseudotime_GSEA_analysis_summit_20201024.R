

###############load functions##############

signal_gsea_with_geneset <- function(genelist)
  
{
  #genelist <- Angiogenesis
  genelist <- intersect(genelist, rownames(matrix_temp_scale))
  genelist_GSVA<-t(gsva(as.matrix(matrix_temp_scale),list(genelist)))
  return(data.frame(genelist_GSVA=genelist_GSVA))
}

convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F)
  humanx <- unique(genesV2[, 2])
  return(genesV2)
}


###############load packages##############

library("GSVA")
library(Seurat)
library("monocle")
library(dplyr)
library(reshape2)
library(tidyr)
library("ggpubr")
library(randomcoloR)
library(RColorBrewer)


###############set the working path##############

setwd("CT26 tumor shBCL9 scRNAseq results_out3")
hsBCL9_pbmc_part_filtered_endo <-readRDS("hsBCL9_pbmc_part_filtered_endo.rds")
hsBCL9_pbmc_part_filtered<- hsBCL9_pbmc_part_filtered_endo[,which(hsBCL9_pbmc_part_filtered_endo$seurat_clusters %in% c(1:4))]

################################################################
#########Figure 2A B D clustering of all cell population###########
################################################################

hsBCL9_pbmc_part_filtered <- ScaleData(hsBCL9_pbmc_part_filtered)
hsBCL9_pbmc_part_filtered <- FindVariableFeatures(hsBCL9_pbmc_part_filtered, selection.method = "vst", nfeatures = 2000)
hsBCL9_pbmc_part_filtered <- ScaleData(hsBCL9_pbmc_part_filtered, features = rownames(hsBCL9_pbmc_part_filtered),block.size = 100,min.cells.to.block = 1000)
hsBCL9_pbmc_part_filtered <- RunPCA(hsBCL9_pbmc_part_filtered)
hsBCL9_pbmc_part_filtered <- FindNeighbors(hsBCL9_pbmc_part_filtered, dims = 1:10)
hsBCL9_pbmc_part_filtered <- RunUMAP(hsBCL9_pbmc_part_filtered, dims = 1:10)
hsBCL9_pbmc_part_filtered <- RunTSNE(hsBCL9_pbmc_part_filtered, dims = 1:10)
hsBCL9_pbmc_part_filtered <- FindClusters(hsBCL9_pbmc_part_filtered, resolution = 0.1)
hsBCL9_pbmc_part_filtered@meta.data[["group2and3"]] <-hsBCL9_pbmc_part_filtered$seurat_clusters %in% c(2,3)


palette_seurat_clusters <- brewer.pal( length(unique(hsBCL9_pbmc_part_filtered$seurat_clusters)),name ="Set2")
palette_Sample <- brewer.pal( length(unique(hsBCL9_pbmc_part_filtered$Sample)),name ="Set3")
palette_Group <- brewer.pal( length(unique(hsBCL9_pbmc_part_filtered$Group)),name ="Dark2")
palette_group2and3 <- brewer.pal( length(unique(hsBCL9_pbmc_part_filtered$Group)),name ="Accent")

pdf(paste0("Figure2_3_hsBCL9_Seurat_umap_cluster_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part_filtered, reduction = "umap",label = TRUE ,cols = palette_seurat_clusters ))
dev.off()

pdf(paste0("Figure2_3_hsBCL9_Seurat_tsne_cluster_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part_filtered, reduction = "tsne",label = TRUE ,cols = palette_seurat_clusters))
dev.off()

pdf(paste0("Figure2_3_hsBCL9_Seurat_umap_Sample_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part_filtered, reduction = "umap", group.by = "Sample",label = TRUE,cols = palette_Sample ))
dev.off()

pdf(paste0("Figure2_3_hsBCL9_Seurat_tsne_Sample_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part_filtered, reduction = "tsne", group.by = "Sample",label = TRUE ,cols = palette_Sample))
dev.off()

pdf(paste0("Figure2_3_hsBCL9_Seurat_umap_Group_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part_filtered, reduction = "umap", group.by = "Group",cols = palette_Group ))
dev.off()

pdf(paste0("Figure2_3_hsBCL9_Seurat_tsne_Group_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part_filtered, reduction = "tsne", group.by = "Group" ,cols = palette_Group))
dev.off()

pdf(paste0("Figure2_3_hsBCL9_Seurat_umap_group2and3_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part_filtered, reduction = "umap", group.by = "group2and3",cols = palette_group2and3))
dev.off()

pdf(paste0("Figure2_3_hsBCL9_Seurat_tsne_group2and3_",Sys.Date(),".pdf"))
print(DimPlot(hsBCL9_pbmc_part_filtered, reduction = "tsne", group.by = "group2and3",cols =palette_group2and3 ))
dev.off()


endo_list<-c("Stab2","Csf1r","Cd3g","Ebf1","Irf8","Sox9","Apoc3","Top2A","Jchain","Dcn","Ptrf","Sdpr","Murc","Prkcdbp","Ehd2","Cav1","Cav2","Cav3","Rarres2","Vwf","Clec4g","Wnt9b","Rspo3","Cdh13","Wnt2","Ednrb","Jag1","Lrg1","Efnb1","Ltbp4","Adgrg6","Fcgr2b","Gpr182","Kdr","Clec4g","Ehd3","Plpp3","Ednrb","Vwf","Colec11","Col14a1","Col3a1","Cxcl12","Aebp1","Prss23","Dcn","Rarres2","Stab2","Csf1r","Cd3g","Ebf1","Irf8","Sox9","Apoc3","Top2A","Jchain","Dcn","Wnt9b","Rspo3","Cdh13","Wnt2","Ednrb","Jag1","Lrg1","Efnb1","Ltbp4","Adgrg6","Fcgr2b","Gpr182")
endo_list<-intersect(rownames(hsBCL9_pbmc_part_filtered),endo_list)

for (plotgene in endo_list)
  
{
  pdf(paste0("Figure2_3_hsBCL9_Seuratgenes_umap_",plotgene,Sys.Date(),".pdf"))
  print(FeaturePlot(hsBCL9_pbmc_part_filtered, features = plotgene,cols = c("grey","red") ))
  dev.off()
  
  pdf(paste0("Figure2_3_hsBCL9_Seuratgenes_tsne_",plotgene,Sys.Date(),".pdf"))
  print(FeaturePlot(hsBCL9_pbmc_part_filtered, reduction = "tsne",features = plotgene,cols = c("grey","red") ))
  dev.off()
  
}


hsBCL9_pbmc_part_filtered.markers <- FindAllMarkers(hsBCL9_pbmc_part_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, num.cores = 64)
hsBCL9_pbmc_part_filtered.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- hsBCL9_pbmc_part_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

pdf(paste0("Figure2_3_hsBCL9_Seurat_DoHeatmap_",Sys.Date(),".pdf"))
print(DoHeatmap(hsBCL9_pbmc_part_filtered, features = top10$gene) + NoLegend())
dev.off()



##############################################################
######################Figure 2C###############################
##############################################################


seurat_clusters_count<-data.frame(seurat_clusters=hsBCL9_pbmc_part_filtered$seurat_clusters,Group=hsBCL9_pbmc_part_filtered$Group)
seurat_clusters_count<-t(table(seurat_clusters_count))
barpt<-ggplot(melt(seurat_clusters_count), aes(fill=as.character(Group), y=value, x=seurat_clusters)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = palette_Group)+
  coord_flip()

pdf(paste0("Figure2C_seurat_clusters_barpt_",Sys.Date(),".pdf"))
barpt
dev.off()

##############################################################
######################Figure 2E plot GSEA#####################
##############################################################


FALSE_GSEA<-read.delim("CT26 tumor shBCL9 scRNAseq results_out3/ready_Ttouse/my_analysis.Gsea.1597976157032/gsea_report_for_FALSE_1597976157032.tsv",header = T)
TRUE_GSEA<-read.delim("CT26 tumor shBCL9 scRNAseq results_out3/ready_Ttouse/my_analysis.Gsea.1597976157032/gsea_report_for_TRUE_1597976157032.tsv",header = T)

FALSE_20<-subset(FALSE_GSEA, GS.DETAILS!="")
FALSE_20$Group<-c("FALSE")
TRUE_20<-subset(TRUE_GSEA, GS.DETAILS!="")
TRUE_20$Group<-c("TRUE")
All<-rbind(FALSE_20, TRUE_20)

GSEA_Plot<-ggplot(All, aes(reorder(NAME, NES), NES, fill=Group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c(palette_group2and3))+#
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA (Endothelial cells)")

pdf(file="Figure2E_hsBCL9_Pathways_NES_from_GSEA.pdf",width=10, height=10)
GSEA_Plot
dev.off()



################################################################
#########Figure 3 A-F Pseudotime analysis ######################
################################################################


patient_name<-"BCL9" 
Selected_population <- "Endo"

pbmc_part <- hsBCL9_pbmc_part_filtered
matrix_temp <- pbmc_part@assays[["RNA"]]@counts
matrix_temp_scale <- pbmc_part@assays[["RNA"]]@data
HSMM <- newCellDataSet(as.matrix(matrix_temp))


HSMM@phenoData@data[["Sample"]] <- pbmc_part@meta.data[["Sample"]]
HSMM@phenoData@data[["Group"]] <- pbmc_part@meta.data[["Group"]]
HSMM@phenoData@data[["seurat_clusters"]] <- pbmc_part@meta.data[["seurat_clusters"]]
HSMM@phenoData@data[["group2and3"]] <- pbmc_part@meta.data[["group2and3"]]


#rpc_matrix <- relative2abs(HSMM, method = "num_genes")
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))

HSMM <- detectGenes(HSMM, min_expr = 0.1)
L <- log(exprs(HSMM[expressed_genes,]))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.01 &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id


HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 20,
                        reduction_method = 'tSNE', verbose = T)

HSMM <- clusterCells(HSMM, num_clusters=2)
HSMM_myo <- estimateDispersions(HSMM)

disp_table <- dispersionTable(HSMM_myo)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.01 &
                           dispersion_empirical >=1 * dispersion_fit)$gene_id

HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2)
HSMM_myo <- orderCells(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo,root_state = 4)


pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_singlemonocleplot.pdf"))
print(plot_cell_trajectory(HSMM_myo))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_singlemonocleplotPseudotime.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by="Pseudotime"))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_singlemonocleplotSample.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "Sample")+ scale_color_manual(breaks = waiver(),values=palette_Sample))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_singlemonocleplotGroup.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "Group")+ scale_color_manual(breaks = waiver(),values=palette_Group))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_singlemonocleplotseurat_clusters.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "seurat_clusters")+ scale_color_manual(breaks = waiver(),values=palette_seurat_clusters))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_singlemonocleplotseurat_group2and3.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by= "group2and3")+ scale_color_manual(breaks = waiver(),values=palette_group2and3))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocleplot.pdf"))
print(plot_complex_cell_trajectory(HSMM_myo))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocleplotPseudotime.pdf"))
print(plot_complex_cell_trajectory(HSMM_myo, color_by="Pseudotime"))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocleplotSample.pdf"))
print(plot_complex_cell_trajectory(HSMM_myo, color_by= "Sample")+ scale_color_manual(breaks = waiver(),values=palette_Sample))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocleplotGroup.pdf"))
print(plot_complex_cell_trajectory(HSMM_myo, color_by= "Group")+ scale_color_manual(breaks = waiver(),values=palette_Group))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocleplotseurat_clusters.pdf"))
print(plot_complex_cell_trajectory(HSMM_myo, color_by= "seurat_clusters")+ scale_color_manual(breaks = waiver(),values=palette_seurat_clusters))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocleplotseurat_group2and3.pdf"))
print(plot_complex_cell_trajectory(HSMM_myo, color_by= "group2and3")+ scale_color_manual(breaks = waiver(),values=palette_group2and3))
dev.off()


endo_list<-c("Kdr","Clec4g","Ehd3","Plpp3","Ednrb","Vwf","Colec11","Col14a1","Col3a1","Cxcl12","Aebp1","Prss23","Dcn","Rarres2","Stab2","Csf1r","Cd3g","Ebf1","Irf8","Sox9","Apoc3","Top2A","Jchain","Dcn","Wnt9b","Rspo3","Cdh13","Wnt2","Ednrb","Jag1","Lrg1","Efnb1","Ltbp4","Adgrg6,Fcgr2b","Gpr182","Acta2","Aebp1","Ccdc80","Cfd","Col12a1","Col1a1","Col1a2","Col3a1","Col5a2","Col6a1","Cox4i2","Ctsk","Dcn","Dpt","Esam","Fstl1","Gja4","Gng11","Gsn","Higd1b","Lum","Mfap5","Mmp11","Myl9","Pdgfra","Postn","Rgs5","Tagln","Thbs2","Vcan")

endo_list<-intersect(rownames(hsBCL9_pbmc_part_filtered),endo_list)

for (plot_gene in endo_list )
  
{
  if (plot_gene %in% rownames(HSMM_myo) )
  {
    pdf(file=paste0(plot_gene,"_",patient_name,Selected_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocle_plotgene.pdf"))
    print(plot_cell_trajectory(HSMM_myo, color_by =(matrix_temp_scale[plot_gene,]))+scale_colour_gradientn(colours = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"),limits=c(min(matrix_temp_scale[plot_gene,]),max(matrix_temp_scale[plot_gene,])))+labs(title=paste("Cluster",Selected_population,plot_gene,patient_name)))
    dev.off()
  }
}


expressed_genes <- row.names(subset(fData(HSMM_myo), num_cells_expressed >= 100))
diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)",cores=8)
sig_gene_names <- row.names(subset(diff_test_res, (qval < 10**-100 & use_for_ordering=="TRUE")))

pdf(file=paste0(Selected_population,patient_name,"_finishImages.pdf"),10,20)
print(plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                              num_clusters = 6,
                              cores = 8,
                              show_rownames = T))
dev.off()


################################################################
#########Figure 3 A-F GSVA analysis ############################
################################################################


bcl9_up_50 <- c("IGFBP7","SPARC","RARRES2","BGN","LOXL1","COL5A2","FSTL1","COL6A2","DCN","MFAP5","SERPING1","AEBP1","GPX3","THY1","MMP2","BMP1","FBN1","ADAMTS2","COL1A1","COL6A3","RCN3","FBLN2","PLPP3","LOXL2","CD248","COL6A1","PDSSTN","RNASE4","COL3A1","COL1A2","COL5A3","C1QTNF6","MGST1","SERPINF1","SOD3","EBF1","EFEMP2","CYGB","SULF1","FXYD1","VCAN","NBL1","FN1","TGFBR2","SERPINA3","SELENOM","MMP14","RCN1","GPX7","BICC1")
bcl9_down_50 <-c("TCP1","RPL12","RPL3","RPL4","TMPO","RPL7","PRKG2","RRM2","LARS2","FCER1G","RAD21","EZR","MTAP","CD9","RPS6","TOP2A","HSP90AB1","HSPA9","MT-CYB","HBEGF","AMIGO2","ACTN4","ACTB","CAVIN2","PLA2G7","CENPF","ATP5F1B","HSPA8","EEF2","TUBA1C","RPS18","ANLN","RAN","WDR31","NOLC1","CPE","TM4SF1","HSPD1","SPP1","PHGDH","TUBA1B","S100A4","CD74","UBE2C","LGALS7","HMGB2","CAV2","ESM1","CCND1","HMGA1")
Angiogenesis <- c("ACVRL1","JAG1","ANGPT1","ANGPT2","CD34","CDC42","MAPK14","TYMP","EDN1","EFNB2","EGR3","EPHB4","PTK2B","FGF2","FGFR1","VEGFD","FOXC2","FLT1","FLT4","FN1","GPLD1","NR4A1","ID1","ITGA5","ITGAV","ITGB1","KDR","LOXL2","MMP14","NOTCH1","PDGFA","PDGFRB","PGF","PIK3CA","PTGS2","PTK2","ROBO1","SHC1","SRF","TAL1","TDGF1","TEK","VAV2","VEGFA","VEGFC","FGF18","NRP1","SEMA5A","RAMP2","CIB1","ESM1","JMJD6","HEY1","ADGRA2","GREM1","SRPX2","SOX18","PARVA","RNF213","E2F8","RSPO3","OTULIN","E2F7","CCBE1","BMPER","NRARP","TNFAIP6","VCAN","SPP1","CCND2","PIK3R1","STC1","JAG2")

plot_gene_list <- convertMouseGeneList(bcl9_up_50)
bcl9_up_50_lin_GSEA<-signal_gsea_with_geneset(plot_gene_list)
bcl9_up_50_GSVA <- bcl9_up_50_lin_GSEA$genelist_GSVA

plot_gene_list <- convertMouseGeneList(bcl9_down_50)
bcl9_down_50_lin_GSEA<-signal_gsea_with_geneset(plot_gene_list)
bcl9_down_50_GSVA <- bcl9_down_50_lin_GSEA$genelist_GSVA

plot_gene_list <- convertMouseGeneList(Angiogenesis)
Angiogenesis_GSEA<-signal_gsea_with_geneset(plot_gene_list)
Angiogenesis_GSVA <- Angiogenesis_GSEA$genelist_GSVA


GSEA_table <- data.frame(Angiogenesis_GSVA=Angiogenesis_GSVA,bcl9_up_50_GSVA=bcl9_up_50_GSVA,bcl9_down_50_GSVA=bcl9_down_50_GSVA)
colnames(GSEA_table) <- c("Angiogenesis_GSVA","bcl9_up_50_GSVA","bcl9_down_50_GSVA")

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocleplotseurat_Angiogenesis_GSVA.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by =(GSEA_table$Angiogenesis_GSVA))+scale_colour_gradientn(colours = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"),limits=c(min(GSEA_table$Angiogenesis_GSVA),max(GSEA_table$Angiogenesis_GSVA)))+labs(title=paste("Cluster",Selected_population,"Angiogenesis_GSVA",patient_name)))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocleplotseurat_Angiogenesis_GSVA.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by =(GSEA_table$bcl9_up_50_GSVA))+scale_colour_gradientn(colours = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"),limits=c(min(GSEA_table$bcl9_up_50_GSVA),max(GSEA_table$bcl9_up_50_GSVA)))+labs(title=paste("Cluster",Selected_population,"Angiogenesis_GSVA",patient_name)))
dev.off()

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocleplotseurat_Angiogenesis_GSVA.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by =(GSEA_table$bcl9_down_50_GSVA))+scale_colour_gradientn(colours = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"),limits=c(min(GSEA_table$bcl9_down_50_GSVA),max(GSEA_table$bcl9_down_50_GSVA)))+labs(title=paste("Cluster",Selected_population,"Angiogenesis_GSVA",patient_name)))
dev.off()


################################################################
#########Figure 3 K BCL9_endo_Score ############################
################################################################


bcl9_up_50_GSVA<-try(t(gsva(as.matrix(matrix_temp_scale),list(convertMouseGeneList(bcl9_up_50)))))
bcl9_down_50_GSVA<-try(t(gsva(as.matrix(matrix_temp_scale),list(convertMouseGeneList(bcl9_down_50)))))
BCL9_endo_Score <- bcl9_down_50_GSVA/(bcl9_up_50_GSVA+0.0000000001)

#Optimize visual adjustment
BCL9_endo_Score[BCL9_endo_Score>5] <-5
BCL9_endo_Score[BCL9_endo_Score<(-5)] <-(-5)
BCL9_endo_Score[BCL9_endo_Score>(-5) & BCL9_endo_Score<5] <-0

pdf(file=paste0(patient_name,Selected_population,Sys.Date(),"_plot_complex_cell_trajectorysinglemonocleplotseurat_BCL9_endo_Score.pdf"))
print(plot_cell_trajectory(HSMM_myo, color_by =(BCL9_endo_Score))+scale_colour_gradientn(colours = c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b"),limits=c(min(BCL9_endo_Score),max(BCL9_endo_Score)))+labs(title=paste("Cluster",Selected_population,"BCL9_endo_Score",patient_name)))
dev.off()



