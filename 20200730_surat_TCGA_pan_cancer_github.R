
library(Seurat)
library(dplyr)
library(scCATCH)

setwd("/TCGA/")
tcga_unused_allcancer_annot_RNA <- read.delim(gzfile("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"), header = T)

#####file downloaded from https://pancanatlas.xenahubs.net/download/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz

rownames(tcga_unused_allcancer_annot_RNA) <- make.unique(tcga_unused_allcancer_annot_RNA$sample) 
tcga_unused_allcancer_annot_RNA_clean <-  tcga_unused_allcancer_annot_RNA[,-1]


pbmc <- CreateSeuratObject(counts = tcga_unused_allcancer_annot_RNA_clean,project = "TCGA", min.cells = 3, min.features = 200)


pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#####
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc)
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(pbmc, dims = 1, reduction = "pca")
print(DimPlot(pbmc, reduction = "pca"))
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:20,n.neighbors=100)

pbmc <- RunTSNE(pbmc, dims = 1:10)

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 3)


tsne_result<-data.frame(pbmc@meta.data$seurat_clusters,pbmc@reductions[["pca"]]@cell.embeddings,pbmc@reductions[["umap"]]@cell.embeddings,pbmc@reductions[["tsne"]]@cell.embeddings)

save(tsne_result,file="TCGA_result")





pdf(paste0("TCGA_Seurat_testumap_",Sys.Date(),".pdf"))

print(DimPlot(pbmc, reduction = "umap",label = TRUE ))
print(DimPlot(pbmc, reduction = "tsne",label = TRUE ))
dev.off()


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, num.cores = 64)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)




save(pbmc, file = "final_one")


pdf(paste0("TCGA_Seurat_DoHeatmap_",Sys.Date(),".pdf"))

print(DoHeatmap(pbmc, features = top10$gene) + NoLegend())

dev.off()


#plot umap of TCGA pan-cancer

library(dbscan)
library("randomcoloR")


load("TCGA_result")

group_information <- read.delim("Survival_SupplementalTable_S1_20171025_xena_sp", header = T)

####Survival_SupplementalTable_S1_20171025_xena_sp file downloaded from https://pancanatlas.xenahubs.ne

tsne_result<-data.frame(sample=gsub("\\.", "-", rownames(tsne_result)),tsne_result)

tsne_result_group_information<-merge(group_information,tsne_result,by="sample")

dbscan_result <- dbscan(tsne_result_group_information[,c("UMAP_1","UMAP_2")], minPts = 0.5,eps =0.2)


n <- length(unique(dbscan_result$cluster ))
palette <- distinctColorPalette(n)
tsne_result_group_information <- data.frame(tsne_result_group_information,dbscan_cluster =as.character(dbscan_result$cluster) )


means_group <- tsne_result_group_information %>%
  group_by(cancer.type.abbreviation)%>%
  summarise(mean_t1 = mean(UMAP_1),
            mean_t2 = mean(UMAP_2))%>% print(n=1)

require("ggrepel")
pdf(paste0(Sys.Date(),"_result_group_information.pdf") )
ggplot(tsne_result_group_information,mapping = aes(x= UMAP_1, y= UMAP_2, color=  dbscan_cluster,label=cancer.type.abbreviation ))+
  geom_point(size = 1)+
  scale_color_manual(values =  palette)+
  geom_label_repel(data = means_group, aes(x = mean_t1,y = mean_t2),label =means_group$cancer.type.abbreviation, colour=1,size = 3,alpha = .5)+
  ggtitle("TCGA_UMAP_DBSCAN_ANALYSIS")

dev.off()




