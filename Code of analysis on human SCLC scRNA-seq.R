library(Seurat)
library(tidyverse)
library(DoubletFinder)
library(RColorBrewer)
library(infercnv)
library(ggplot2)
library(plotrix)
library(SCENIC)
library(SCopeLoomR)


####Data preprocessing and quality control####
sample_names <- c("LC51", "LC242", "LC250", "LC251", "LC254", "LC255", "LC270")
sample_paths <- paste0("D:/XCL/", sample_names, "/filtered_feature_bc_matrix")

seurat_list <- list()
for (i in 1:length(sample_names)) {
  data <- Read10X(sample_paths[i])
  seurat_object <- CreateSeuratObject(data, project = sample_names[i], min.cells = 3, min.features = 200)
  seurat_object$Sample <- sample_names[i]
  seurat_list[[sample_names[i]]] <- seurat_object
}

for (name in names(seurat_list)) {
  seurat_object <- seurat_list[[name]]
  
  seurat_object <- PercentageFeatureSet(seurat_object, "^MT-", col.name = "percent_mito")
  seurat_object <- PercentageFeatureSet(seurat_object, "^RP[SL]", col.name = "percent_ribo")
  seurat_object <- subset(seurat_object, subset = percent_mito < 15 & percent_ribo > 5)
  seurat_object <- seurat_object[!grepl("^MT-", rownames(seurat_object)), ]
  
  seurat_list[[name]] <- seurat_object
}

#remove doublets by DoubletFinder
for (i in 1:length(seurat_list)) {
  seurat_object <- NormalizeData(seurat_list[[i]])
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object, vars.to.regress = c("nFeature_RNA", "percent_mito"))
  seurat_object <- RunPCA(seurat_object, npcs = 30)
  seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:30)
  seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:30)              
  seurat_object <- FindClusters(seurat_object, resolution = 0.9)
  
  sweep.res <- paramSweep_v3(seurat_object)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  pk_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  DoubleRate <- ncol(seurat_object)*8*1e-6
  homotypic.prop <- modelHomotypic(seurat_object$seurat_clusters) 
  nExp_poi <- round(DoubleRate*ncol(seurat_object))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  seurat_object <- doubletFinder_v3(seurat_object, 
                                    PCs = 1:30, 
                                    pN = 0.25, 
                                    pK = pk_bcmvn, 
                                    nExp = nExp_poi.adj, 
                                    reuse.pANN = FALSE, 
                                    sct = FALSE)
  
  DF.name <- colnames(seurat_object@meta.data)[grepl("DF.classification", colnames(seurat_object@meta.data))]
  cowplot::plot_grid(ncol = 2, DimPlot(seurat_object, group.by = "seurat_clusters", label = T) + NoAxes(), 
                     DimPlot(seurat_object, group.by = DF.name) + NoAxes())
  VlnPlot(seurat_object, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
  
  seurat_object <- seurat_object[, seurat_object@meta.data[, DF.name] == "Singlet"]
}

for (name in names(seurat_list)) {
  seurat_object <- seurat_list[[name]]
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
  seurat_list[[name]] <- seurat_object
}


####Anchor-based sample integration####
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(seurat_list)
anchors <- FindIntegrationAnchors(seurat_list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.9)
DefaultAssay(combined) <- "RNA"

combined@meta.data <- combined@meta.data %>%
  mutate(new_cluster = factor(as.numeric(as.character(seurat_clusters)) + 1, levels = 1:26))
Idents(combined) <- "new_cluster"

combined@meta.data <- combined@meta.data %>%
  mutate(celltype = case_when(
    new_cluster %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 14, 19, 22, 25, 26) ~ "Tumor",
    new_cluster == 20 ~ "AT2",
    new_cluster %in% c(10, 16, 17, 24) ~ "T",
    new_cluster == 18 ~ "Fibroblast",
    new_cluster == 15 ~ "Myeloid",
    new_cluster == 23 ~ "Endothelial",
    new_cluster == 21 ~ "B",
    new_cluster == 12 ~ "Plasma",
    ))

#remove the low-quality immune and stromal cells with expressions of epithelial or tumor markers
Idents(combined) <- "celltype"
nonepi <- subset(combined, idents = c("T","Fibroblast","Myeloid","Endothelial","B","Plasma"))
nonepi2 <- subset(nonepi, EPCAM == 0 & INSM1 == 0 & ASCL1 == 0 & NEUROD1 == 0)

nonepi@meta.data <- nonepi@meta.data %>%
  mutate(celltype = if_else(rownames(nonepi@meta.data) %in% rownames(nonepi2@meta.data),
                            nonepi2@meta.data$celltype[match(rownames(nonepi@meta.data), rownames(nonepi2@meta.data))],
                            "Removed"))
combined@meta.data <- combined@meta.data %>%
  mutate(celltype = if_else(rownames(combined@meta.data) %in% rownames(nonepi@meta.data),
                            nonepi@meta.data$celltype[match(rownames(combined@meta.data), rownames(nonepi@meta.data))],
                            celltype))
Idents(combined) <- "celltype"
combined <- subset(combined, idents = "Removed", invert = TRUE)


####Dimension reduction, clustering and cell type annotation####
DefaultAssay(combined) <- "RNA"
combined <- FindVariableFeatures(combined)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.9)
DefaultAssay(combined) <- 'RNA'

combined@meta.data <- combined@meta.data %>%
  mutate(new_cluster = factor(as.numeric(as.character(seurat_clusters)) + 1, levels = 1:24))
Idents(combined) <- "new_cluster"

#cell type annotation
combined@meta.data <- combined@meta.data %>%
  mutate(celltype = case_when(
    new_cluster %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 17, 24) ~ "Tumor",
    new_cluster == 20 ~ "AT2",
    new_cluster %in% c(11, 15, 16, 23) ~ "T",
    new_cluster == 19 ~ "Fibroblast",
    new_cluster == 18 ~ "Myeloid",
    new_cluster == 22 ~ "Endothelial",
    new_cluster == 21 ~ "B",
    new_cluster == 14 ~ "Plasma",
  ))


####Copy number variation analysis####
combined_infercnv <- subset(combined, downsample = 300)
raw_counts_matrix <- as.data.frame(combined@assays[["RNA"]]@counts)
annotations_file <- combined@meta.data %>% rownames_to_column("cellname") %>% dplyr::select(c("cellname","new_cluster"))
write.table(annotations_file, "infercnv/annotations_file.txt", sep = '\t', row.names = F, col.names = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_counts_matrix,
                                    annotations_file="infercnv/annotations_file.txt",
                                    delim="\t",
                                    gene_order_file="gencode_v19_gene_pos.txt",
                                    ref_group_names=c("11","19","21"))

infercnv_obj = infercnv::run(infercnv_obj, cutoff=0.1, out_dir="infercnv", 
                             cluster_by_groups=TRUE, denoise=T, HMM=T)


####Visualization####
pdf("FigureS1A_dimplot_combined_sample.pdf", height = 8, width = 10)
DimPlot(combined, group.by = "Sample", pt.size = 0.1, cols = brewer.pal(7,"Set3"))
dev.off()


pdf("FigureS1B_dimplot_combined_cluster.pdf", height = 8, width = 10)
DimPlot(combined, group.by = "new_cluster",  pt.size = 0.1, label.box = T, label = T, 
        label.size = 2, repel = T, cols = c(brewer.pal(12,"Paired"),brewer.pal(12,"Set3")))
dev.off()


markers <- c("EPCAM","INSM1","NCAM1","UCHL1","ASCL1","NEUROD1","SFTPB","NAPSA",
           "MYC","VIM","PTPRC","CD3D","CD79A","MZB1","CD14","ACTA2","PECAM1")
Idents(combined) <- factor(combined$new_cluster, 
                            levels = c(20,1,2,3,4,5,6,7,8,9,10,12,17,24,13,11,15,16,23,21,14,18,19,22))
pdf("FigureS1C_dotplot_combined_marker.pdf",height = 8,width = 12)
DotPlot(combined, features = rev(markers), dot.scale = 6) +
  coord_flip() +
  scale_colour_gradient(low = "white", high = "#08519C")+
  theme(panel.border = element_rect(color = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()


Idents(combined) <- "new_cluster"
Cellratio <- prop.table(table(Idents(combined), 
                              factor(combined$Sample, 
                                     levels = rev(levels(combined$Sample)))), margin = 2) %>% as.data.frame()
pdf("FigureS1E_barplot_combined_cluster_split_by_sample.pdf", height = 5, width = 10)
ggplot(Cellratio) +
  geom_bar(aes(x=Var2, y=Freq, fill=Var1), stat = "identity", width = 0.8, size=0, colour = 'white') +
  theme_classic() +
  labs(title = "All cells", x='Samples', y='Proportions')+
  scale_fill_manual(values = c(brewer.pal(12,"Paired"),brewer.pal(12,"Set3")))+
  coord_flip() +
  theme(panel.border = element_rect(fill = NA, color="black", size = 0.5, linetype="solid"))
dev.off()


Idents(combined) <- factor(combined$celltype, 
                           levels = c("AT2","Tumor","T","B","Plasma","Myeloid","Fibroblast","Endothelial"))
names <- table(combined$celltype) %>% names()
ratio <- table(combined$celltype) %>% as.numeric()
pielabel <- paste0(names," (", round(ratio/sum(ratio)*100,2), "%)")
pdf("FigureS1F_Pieplot_combined_celltype.pdf", height = 5, width = 5)
pie3D(ratio,labels = pielabel,explode = 0.2, 
      col = brewer.pal(8,"Paired"), theta = pi/3,
      height = 0.1, labelcex = 0.9)
dev.off()


####Tumor cells subclustering####
tumor <- subset(combined, idents = "Tumor")

DefaultAssay(tumor) <- "RNA"
tumor <- FindVariableFeatures(tumor)
DefaultAssay(tumor) <- "integrated"
tumor <- ScaleData(tumor, verbose = FALSE)
tumor <- RunPCA(tumor, npcs = 30, verbose = FALSE)
tumor <- RunUMAP(tumor, reduction = "pca", dims = 1:30)
tumor <- FindNeighbors(tumor, reduction = "pca", dims = 1:30)
tumor <- FindClusters(tumor, resolution = 0.6)
DefaultAssay(tumor) <- 'RNA'

#reorder the tumor cluster
tumor@meta.data <- tumor@meta.data %>% 
  mutate(tumor_cluster = recode(seurat_clusters,
                            "0" = "1", "1" = "2", "2" = "3", "3" = "4", "4" = "5", 
                            "5" = "6", "7" = "7", "8" = "8", "9" = "9", "12" = "10", 
                            "11" = "11", "13" = "12", "15" = "13", "14" = "14", 
                            "6" = "15", "16" = "16", "10" = "17"))
tumor$tumor_cluster <- factor(tumor$tumor_cluster, levels = c(1:17))
Idents(tumor) <- "tumor_cluster"


####Visualization####
pdf("Figure1A_dimplot_tumor_cluster.pdf", height = 8, width = 10)
DimPlot(tumor, group.by = "tumor_cluster",  pt.size = 0.2, 
        cols = c(brewer.pal(10,"Paired"), brewer.pal(6, "Set3"), "#004D7A"),
        label.box = T, label = T, label.size = 4, repel = T)
dev.off()


markers <- c("EPCAM", "INSM1", "NCAM1", "UCHL1", "CALCA", 
           "ASCL1", "NEUROD1", "POU2F3", "YAP1", 
           "SFTPB", "SFTPC", "HLA-DRA", "SCGB1A1", "SCGB3A1", "NAPSA", 
           "MUC1", "KRT7", "KRT5", "TP63", 
           "CDH1", "VIM", "MYC", "MYCL", "MYCN")
pdf("Figure1B_dotplot_tumor_marker.pdf",height = 7,width = 8)
DotPlot(tumor, features = rev(markers), dot.scale = 6) +
  coord_flip() +
  scale_colour_gradient(low = "white", high = "#08519C")+
  theme(panel.border = element_rect(color = "black", size = 0.8))
dev.off()


Idents(tumor) <- "tumor_cluster"
Cellratio <- prop.table(table(Idents(tumor), 
                              factor(tumor$Sample, levels = rev(levels(tumor$Sample)))), margin = 2) %>%
  as.data.frame()
pdf("FigureS1G_barplot_tumor_cluster_split_by_sample.pdf", height = 5, width = 8)
ggplot(Cellratio) +
  geom_bar(aes(x=Var2, y=Freq, fill=Var1), stat = "identity", width = 0.8, size=0, colour = 'white') +
  theme_classic() +
  labs(title = "Tumor cells", x='Samples', y='Proportions')+
  scale_fill_manual(values = c(brewer.pal(10,"Paired"), brewer.pal(6, "Set3"), "#004D7A"))+
  coord_flip() +
  theme(panel.border = element_rect(fill = NA, color="black", size = 0.5, linetype="solid"))
dev.off()


####Transcription factor regulatory activity analysis####
#produce input file for pySCENIC analysis
Idents(tumor) <- "tumor_cluster"
set.seed(1234)
tumor_pyscenic <- subset(tumor, downsample = 300)
write.csv(t(as.matrix(tumor_pyscenic@assays$RNA@counts)), file = "pyscenic/pyscenic_input_tumor.CSV")

#extract regulon activity matrix from the output loom file
loom <- open_loom('pyscenic/pyscenic_output_tumor.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
close_loom(loom)

sub_regulonAUC <- regulonAUC[,match(colnames(tumor_pyscenic), colnames(regulonAUC))]
identical(colnames(sub_regulonAUC),colnames(tumor_pyscenic))
cellClusters <- data.frame(row.names = colnames(tumor_pyscenic), 
                           Cluster = as.character(tumor_pyscenic$tumor_cluster))

selectedResolution <- "Cluster"
cellsPerGroup <- split(rownames(cellClusters), cellClusters[,selectedResolution])
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),]

regulonactivity_bygroup <- sapply(cellsPerGroup, 
                                  function(cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
regulonactivity_bygroup_scaled <- t(scale(t(regulonactivity_bygroup), center = T, scale = T))
regulonactivity_bygroup_scaled <- regulonactivity_bygroup_scaled[]
regulonactivity_bygroup_scaled <- na.omit(regulonactivity_bygroup_scaled)

write.table(regulonactivity_bygroup_scaled, "pyscenic/tumor_scenic_scaled_regulon_activity.txt", sep = '\t')


#visualization of top TFs of tumor_cluster 17
rss <- regulonactivity_bygroup_scaled
df = do.call(rbind, 
             lapply(1:ncol(rss), function(i){
  dat = data.frame(
    path = rownames(rss),
    cluster = colnames(rss)[i],
    sd.1 = rss[,i],
    sd.2 = apply(rss[,-i], 1,mean)
    )
  }))
df$fc = df$sd.1 - df$sd.2
top <- df %>% group_by(cluster) %>% top_n(20, fc)
top_17 <- subset(top, cluster == "17")
top_17 <- top_17[order(top_17$fc, decreasing = T), ]
n = rss[top_17$path,]
sorted <- colnames(n)[order(as.numeric(colnames(n)))]
n_sorted <- n[, sorted]
pdf("Figure1C_heatmap_topTF_of_tumor_cluster17.pdf", height = 5, width = 5)
pheatmap::pheatmap(n_sorted, show_rownames = T, cluster_cols = F, cluster_rows = F, angle_col = "0") 
dev.off()


pdf("Figure1D_vlnplot_tumor_MYC_exp.pdf", height = 6, width = 8)
VlnPlot(tumor, features = "MYC", pt.size = 0, 
        cols = c(brewer.pal(10,"Paired"), brewer.pal(6, "Set3"), "#004d7a"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))
dev.off()