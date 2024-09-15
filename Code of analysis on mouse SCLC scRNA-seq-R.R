library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(readr)
library(stringr)
library(tidyverse)
library(DoubletFinder)
library(RColorBrewer)
library(homologene)
library(grid)
library(tibble)
library(decontX)
library(pheatmap)
### sample list
samples <- read_excel("metadata - KPM.xlsx", range = cell_cols("A:A")) %>% .$sample_id 

### import cellranger files from different data sets
for (i in seq_along(samples)){
  assign(paste0("scs_data", i), Read10X(data.dir = paste0(samples[i], "/filtered_feature_bc_matrix")))
}

### create seurat objects from cellranger files
for (i in seq_along(samples)){
  assign(paste0("seu_obj", i), CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), 
                                                  project = samples[i], min.cells = 3))
}

### merge data sets
seu_obj <- merge(seu_obj1, y = c(seu_obj2), add.cell.ids = samples, project = "lung")

### calculate mitochondrial, hemoglobin and ribosomal gene counts
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^mt-", col.name = "Mt")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^Hb[ba]-", col.name = "Hb")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^Rp[sl]", col.name = "Rp")

qcparams <- c("nFeature_RNA", "nCount_RNA", "Mt", "Hb", "Rp")
for (i in seq_along(qcparams)){
  print(VlnPlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident", pt.size = 0))
} 

for (i in seq_along(qcparams)){
  print(RidgePlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident"))
} 

RidgePlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "Mt","Hb", "Rp"), group.by = "orig.ident")
ggsave2("QC.pdf", path = "results", width = 30, height = 20, units = "cm")

### clear environment
remove(seu_obj1)
remove(seu_obj2)

remove(scs_data1)
remove(scs_data2)

## Before filtering
seu_obj_unfiltered <- seu_obj

nFeature_lower <- 200
nFeature_upper <- 6500
nCount_lower <- 500
nCount_upper <- 40000
Mt_lower <- 0
Mt_upper <- 25
Hb_lower <- 0
Hb_upper <- 5
## After filtering
seu_obj <- subset(seu_obj_unfiltered, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & 
                    nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & Mt < Mt_upper & Hb < Hb_upper)


save(seu_obj_unfiltered, file = "Before filtering.Rda")
save(seu_obj, file = "After filtering.Rda")

load("Before filtering.Rda")
####doublet####
table(seu_obj_unfiltered$orig.ident)
# KPM257 KPM261 
# 16385  16742  
table(seu_obj$orig.ident)
# KPM257 KPM261 
# 16221  16423

seu_split <- SplitObject(seu_obj, split.by = "orig.ident") 
for (i in 1:length(seu_split)) {
  seu_sample <- NormalizeData(seu_split[[i]])
  seu_sample <- FindVariableFeatures(seu_sample)
  seu_sample <- ScaleData(seu_sample)
  seu_sample <- RunPCA(seu_sample)
  

  stdv <- seu_sample[["pca"]]@stdev
  sum.stdv <- sum(seu_sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  # KPM257 [1] 13 KPM261 [1] 15
  
  seu_sample <- RunUMAP(seu_sample, dims = 1:min.pc)
  seu_sample <- FindNeighbors(object = seu_sample, dims = 1:min.pc)              
  seu_sample <- FindClusters(object = seu_sample, resolution = 0.1)
  
  # pK Identification (no ground-truth)
  sweep.list <- paramSweep(seu_sample, PCs = 1:min.pc)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- seu_sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(seu_sample@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  seu_sample <- doubletFinder(seu = seu_sample, 
                              PCs = 1:min.pc, 
                              pK = optimal.pk,
                              nExp = nExp.poi.adj)
  metadata <- seu_sample@meta.data
  colnames(metadata)[10] <- "doublet_finder"
  seu_sample@meta.data <- metadata
  
  
  samples <- read_excel("metadata - KPM.xlsx", range = cell_cols("A:A")) %>% .$sample_id 
  cowplot::plot_grid(ncol = 1, DimPlot(seu_sample, group.by = "doublet_finder"),
                     VlnPlot(seu_sample, features = "nFeature_RNA", group.by = "doublet_finder", pt.size = 0.1))
  ggsave2(paste0(samples[[i]],"_doublet_finder_result.pdf"), path = "results", width = 20, height = 32, units = "cm")
  
  # subset and save
  seu_singlets <- subset(seu_sample, doublet_finder == "Singlet")
  seu_split[[i]] <- seu_singlets
  remove(seu_singlets)
}

seu_singlets <- merge(x = seu_split[[1]],
                      y = c(seu_split[[2]]),
                      project = "lung")
save(seu_singlets, file = "seu_singlets.Rda")

table(seu_singlets$orig.ident)
# KPM257 KPM261 
# 14241  13580
remove(seu_split)
remove(seu_singlets)
#decontX

library(decontX)
seu_split <- SplitObject(seu_singlets, split.by = "orig.ident") 
for (i in 1:length(seu_split)) {
  counts <- seu_split[[i]]@assays$RNA@counts
  decontX_results <- decontX(counts) 
  seu_split[[i]]$Contamination =decontX_results$contamination
  seu_split[[i]] = seu_split[[i]][,seu_split[[i]]$Contamination < 0.2]
}


seu_singlets <- merge(x = seu_split[[1]],
                      y = c(seu_split[[2]]),
                      project = "lung")
seu_singlets_decontX <-seu_singlets
save(seu_singlets_decontX, file = "seu_singlets_decontX.Rda")
table(seu_singlets_decontX$orig.ident)
# KPM257 KPM261 
# 11247  10446 

seu_obj<-seu_singlets_decontX
seu_obj<- seu_obj %>%  
  NormalizeData(.) %>%
  FindVariableFeatures(.) %>%
  ScaleData(.) %>%
  RunPCA(.) 

stdv <- seu_obj[["pca"]]@stdev
sum.stdv <- sum(seu_obj[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                     percent.stdv[2:length(percent.stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min.pc <- min(co1, co2)
min.pc
# [1] 23
seu_obj<- seu_obj  %>%
  RunUMAP(., dims = 1:min.pc) %>%
  RunTSNE(., dims = 1:min.pc) %>%
  FindNeighbors(., dims = 1:min.pc) 
for (i in c(0.05,0.1, 0.2)) {
  seu_obj <- FindClusters(seu_obj, resolution = i)
  DimPlot(seu_obj, reduction = "umap", raster = F) + labs(title = paste0("resolution: ", i))
  ggsave2(paste0("Dimplot_findclusters_resolution_",i,".pdf"), path = "results", width = 15, height = 15, units = "cm")
}

table(seu_obj$orig.ident)
# KPM257 KPM261 
# 11247  10446

Idents(seu_obj)<-"RNA_snn_res.0.05"

mainmarkers <- c("Epcam", "Napsa", "Sftpc", "Ptprc","Pecam1", "Cd34", "Cdh5", "Acta2", "Pdgfra","Pdpn")
cowplot::plot_grid(ncol = 2, 
                   DimPlot(seu_obj, 
                           group.by = "RNA_snn_res.0.05",
                           cols = c(brewer.pal(9,"Paired")),
                           label = T,
                           label.box = T,
                           repel = F,
                   ),
                   DotPlot(seu_obj, 
                           features = rev(mainmarkers), 
                           group.by = "RNA_snn_res.0.05") +
                     coord_flip() +
                     scale_colour_gradient(low = "white", high = "#08519C"))


ggsave2("DotPlot_mainmarkers.pdf", path = "results", width = 16, height = 8)

abc <- seu_obj@meta.data
abc <- abc %>% mutate(mainmarkers = orig.ident)
abc$mainmarkers<-as.character(abc$mainmarkers)
class(abc$mainmarkers)
abc$mainmarkers[which(abc$RNA_snn_res.0.05==0)]="6"
abc$mainmarkers[which(abc$RNA_snn_res.0.05==1)]="2"
abc$mainmarkers[which(abc$RNA_snn_res.0.05==2)]="1"
abc$mainmarkers[which(abc$RNA_snn_res.0.05==3)]="3"
abc$mainmarkers[which(abc$RNA_snn_res.0.05==4)]="5"
abc$mainmarkers[which(abc$RNA_snn_res.0.05==5)]="4"


seu_obj@meta.data <- abc
Idents(seu_obj)<-"mainmarkers"

DimPlot(seu_obj, 
        group.by = "mainmarkers",
        cols = c("pink2","#bcbddc","#FE9929","#9FB6CD","coral3","lightcyan3"),
        label = T, label.size = 2,
        repel = F,label.box = T,pt.size =0.001)+theme(plot.title = element_text(size = 10, face = "bold"),
                                                      legend.text=element_text(size=6),
                                                      legend.key.height = unit(0.4, 'cm'),
                                                      legend.key.width = unit(0.2, 'cm'),
                                                      axis.text.x = element_text(size = 6),
                                                      axis.text.y = element_text(size = 6))
ggsave2("mainmarkers-umap.pdf", path = "results", width = 5, height = 5)
DotPlot(seu_obj, features = rev(mainmarkers), group.by = "mainmarkers",
        dot.scale = 4) +coord_flip() +
  scale_colour_gradient(low = "white", high = "#08519C")+theme(
    panel.border = element_rect(color = "black", size = 0.8))+theme(plot.title = element_text(size = 6, face = "bold"),
                                                                    legend.title=element_text(angle = 90,hjust = 0.5,vjust=0.5,size=5), 
                                                                    legend.text=element_text(size=5),
                                                                    legend.key.height = unit(0.4, 'cm'),
                                                                    legend.key.width = unit(0.2, 'cm'),
                                                                    axis.text.x = element_text(size = 6),
                                                                    axis.text.y = element_text(size = 6),
                                                                    axis.line.y=element_line(color="white",size=0),
                                                                    axis.line.x=element_line(color="white",size=0))
ggsave2("mainmarkers-dotplot.pdf", path = "results", width = 6, height = 9, units = "cm")


abc <- seu_obj@meta.data
abc <- abc %>% mutate(celltype = orig.ident)
abc$celltype<-as.character(abc$celltype)
class(abc$celltype)
abc$celltype[which(abc$mainmarkers==1)]="Epithelial"
abc$celltype[which(abc$mainmarkers==2)]="Epithelial"
abc$celltype[which(abc$mainmarkers==3)]="Immune"
abc$celltype[which(abc$mainmarkers==4)]="Endothelial"
abc$celltype[which(abc$mainmarkers==5)]="Fibroblast"
abc$celltype[which(abc$mainmarkers==6)]="Fibroblast"

seu_obj@meta.data <- abc
Idents(seu_obj)<-"celltype"

DimPlot(seu_obj, 
        group.by = "celltype",
        cols = c("salmon","lightcyan3","pink2","#bcbddc"),
        label = T, label.size = 2,
        repel = F,pt.size =0.001)+theme(plot.title = element_text(size = 10, face = "bold"),
                                        legend.text=element_text(size=6),
                                        legend.key.height = unit(0.4, 'cm'),
                                        legend.key.width = unit(0.2, 'cm'),
                                        axis.text.x = element_text(size = 6),
                                        axis.text.y = element_text(size = 6))
ggsave2("celltype-umap.pdf", path = "results-7-5", width = 5, height = 5)

save(seu_obj,file = "seu_obj_celltypemain-7-9.Rda")
load("seu_obj_celltypemain-7-9.Rda")
epi<-subset(seu_obj,idents = "Epithelial")

epi<- epi %>%  
  NormalizeData(.) %>%
  FindVariableFeatures(.) %>%
  ScaleData(.) %>%
  RunPCA(.) 

library(harmony)

epi<- epi %>% 
  RunHarmony(.,dims.use = 1:23,
             group.by.vars = "orig.ident",
             plot_convergence = TRUE) %>% 
  RunUMAP(.,  reduction = "harmony", dims = 1:23) %>% 
  RunTSNE(.,  reduction = "harmony", dims = 1:23) %>% 
  FindNeighbors(., reduction = "harmony", dims = 1:23)

for (i in c(0.1, 0.2)) {
  epi <- FindClusters(epi, resolution = i)
  DimPlot(epi, reduction = "umap", raster = F) + labs(title = paste0("resolution: ", i))
  ggsave2(paste0("Epi_Dimplot_findclusters_resolution_harmony",i,".pdf"), path = "results-7-31", width = 15, height = 15, units = "cm")
}


Idents(epi)<-"RNA_snn_res.0.1"

myc <- c('Foxj1', 'Ccdc78', 'Scgb1a1', 'Scgb3a1', 'Krt5', 'Trp63', 'Krt14', 'Creb5', 'Sox9', 'Sox2', 'Ager', 
         'Krt7', 'Epcam', 'Sftpa1', 'Sftpb', 'Sftpc', 'H2-Ab1', 'Napsa', 'Muc1', 'Nkx2-1', 'Uchl1', 'Insm1', 
         'Ncam1', 'Calca', 'Chga', 'Ascl1', 'Pou2f3', 'Yap1', "Hnf4a",'Cdh1', 'Vim', 'Myc', 'Kras', 'Trp53', 'Rb1')
DimPlot(epi, 
        group.by = "RNA_snn_res.0.1",
        cols = c("#FE9929","paleturquoise3","salmon","#9FB6CD","coral3","yellowgreen","pink2"),
        label = T, label.size = 2,
        repel = F)+theme(plot.title = element_text(size = 6, face = "bold"),
                         legend.text=element_text(size=6),
                         legend.key.height = unit(0.4, 'cm'),
                         legend.key.width = unit(0.2, 'cm'),
                         axis.text.x = element_text(size = 6),
                         axis.text.y = element_text(size = 6))


P1 <- DotPlot(epi, features = rev(myc), group.by = "RNA_snn_res.0.1",
              dot.scale = 4) +coord_flip() +
  scale_colour_gradient(low = "white", high = "#08519C")+theme(
    panel.border = element_rect(color = "black", size = 0.8))+theme(plot.title = element_text(size = 6, face = "bold"),
                                                                    legend.title=element_text(angle = 90,hjust = 0.5,vjust=0.5,size=5), 
                                                                    legend.text=element_text(size=5),
                                                                    legend.key.height = unit(0.4, 'cm'),
                                                                    legend.key.width = unit(0.2, 'cm'),
                                                                    axis.text.x = element_text(size = 6),
                                                                    axis.text.y = element_text(size = 6),
                                                                    axis.line.y=element_line(color="white",size=0),
                                                                    axis.line.x=element_line(color="white",size=0))
P1
ggsave2("dot-epi-umap.pdf", path = "results-7-31", width = 10, height = 15,units = "cm")

abc <- epi@meta.data
abc <- abc %>% mutate(celltype = orig.ident)
abc$celltype<-as.character(abc$celltype)
abc$celltype[which(abc$RNA_snn_res.0.1==0)]="LUAD"
abc$celltype[which(abc$RNA_snn_res.0.1==1)]="AT2"
abc$celltype[which(abc$RNA_snn_res.0.1==2)]="Intermediate"
abc$celltype[which(abc$RNA_snn_res.0.1==3)]="SCLC"
abc$celltype[which(abc$RNA_snn_res.0.1==4)]="Club"
abc$celltype[which(abc$RNA_snn_res.0.1==5)]="Ciliated"
abc$celltype[which(abc$RNA_snn_res.0.1==6)]="LUAD"

epi@meta.data <- abc
epi@meta.data$celltype <- factor(epi@meta.data$celltype, levels = c("SCLC", "Intermediate", "LUAD", "AT2", "Club", "Ciliated"))
Idents(epi)<-"celltype"
DimPlot(epi, 
        group.by = "celltype",
        cols = c("yellowgreen","#bcbddc","#9FB6CD","pink2","#FE9929","lightcyan3"),
        label = F, label.size = 5,
        repel = F)+theme(plot.title = element_text(size = 16, face = "bold"),
                         legend.text=element_text(size=12),
                         legend.key.height = unit(0.6, 'cm'),
                         legend.key.width = unit(0.3, 'cm'),
                         axis.text.x = element_text(size = 8),
                         axis.text.y = element_text(size = 8))
ggsave2("celltype_epi_7-28-not-ondata.pdf", path = "results-7-5", width = 18, height = 15,units = "cm")

P1 <- DotPlot(epi, features = rev(myc), group.by = "celltype",
              dot.scale = 4) +coord_flip() +
  scale_colour_gradient(low = "white", high = "#08519C")+theme(
    panel.border = element_rect(color = "black", size = 0.8))+theme(plot.title = element_text(size = 6, face = "bold"),
                                                                    legend.title=element_text(angle = 90,hjust = 0.5,vjust=0.5,size=5), 
                                                                    legend.text=element_text(size=5),
                                                                    legend.key.height = unit(0.4, 'cm'),
                                                                    legend.key.width = unit(0.2, 'cm'),
                                                                    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
                                                                    axis.text.y = element_text(size = 6),
                                                                    axis.line.y=element_line(color="white",size=0),
                                                                    axis.line.x=element_line(color="white",size=0))

P1
ggsave2("dotplot_epi_.pdf", path = "results-7-5", width = 6, height = 15,units = "cm")
save(epi,file = "seu_obj_epi-7-9.Rda")
load("seu_obj_epi-7-9.Rda")

####histogram####
library(scRNAtoolVis)
library(tidydr)
Idents(epi) <- "celltype"
table(Idents(epi), epi$orig.ident)
Cellratio <- prop.table(table(epi@meta.data$celltype,epi@meta.data$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype","sample","ratio") 
library(ggplot2)
library(ggalluvial)
Cellratio$celltype <- factor(Cellratio$celltype, levels = c("SCLC", "Intermediate", "LUAD", "AT2", "Club", "Ciliated"))
color_cluster=c("yellowgreen","#bcbddc","#9FB6CD","pink2","#FE9929","lightcyan3") 

names(color_cluster)=c()
library(scales)
show_col(color_cluster) 
library(Seurat)


Cellratio$sample <- factor(Cellratio$sample,levels = c("KPM257","KPM261"))  

p<-ggplot(Cellratio,aes(x=sample,y=ratio,fill=celltype,stratum=celltype,alluvium=celltype))+
  geom_col(width = 0.7,color=NA)+
  scale_fill_manual(values=color_cluster)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave(plot = p,"celltype_ratio4.pdf",width = 8,height = 8,units = "cm")
ggsave(plot = p,"celltype_ratio5.pdf",width = 15,height = 6,units = "cm")
p


####inferCNV####
library(infercnv)
library(rjags)
library(cowplot)
library(ggdendro)
library(dendextend)
seu_obj <- epi

epi_CNV <- subset(seu_obj, downsample = 500)
raw_counts_matrix <- as.data.frame(epi_CNV@assays[["RNA"]]@counts)
epi_CNV@meta.data$cell_id <- rownames(epi_CNV@meta.data)
annotations_file <- epi_CNV@meta.data %>% dplyr::select(c("cell_id","celltype"))
write.table(annotations_file, "inferCNV_epi_8-21/infercnv_annotations_file.txt", sep = '\t', row.names = F, col.names = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_counts_matrix,
                                    annotations_file="inferCNV_epi_8-21/infercnv_annotations_file.txt",
                                    delim="\t",
                                    gene_order_file="inferCNV_epi_8-21/mm_geneOrderingFile.txt",
                                    ref_group_names= c("Ciliated","Club","AT2"),
                                    chr_exclude = c("MT", "X", "Y")) 


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="inferCNV_epi_8-21", 
                             cluster_by_groups=TRUE, 
                             denoise=F,
                             HMM=T)


save(infercnv_obj, file = "inferCNV_epi_8-21/infercnv_results.Rda")

####GSVA####
library(GSVA)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(circlize)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(SCENIC)
library(SCopeLoomR)
library(GSEABase)

hallmark <- getGmt("mh.all.v2023.2.Mm.symbols.gmt")


LUAD_Inter_SCLC<-c("LUAD","Intermediate","SCLC")
LUAD_Inter_SCLC<-subset(epi,idents = LUAD_Inter_SCLC)



expr <- AverageExpression(LUAD_Inter_SCLC, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0, ] 
expr <- as.matrix(expr)


hallmark_result <- gsva(expr, hallmark, method="ssgsea")


rownames(hallmark_result) <- gsub(rownames(hallmark_result), pattern = "HALLMARK_", replacement = "", fixed = TRUE)

colnames(hallmark_result) <- gsub(colnames(hallmark_result), pattern = "-", replacement = "_", fixed = TRUE)

exprTable_t <- as.data.frame(t(hallmark_result))

col_dist = dist(exprTable_t)

hclust_1 <- hclust(col_dist)

pheatmap(hallmark_result, cluster_cols = hclust_1)

manual_order = c("LUAD", "Intermediate", "SCLC")

dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))))

dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))), agglo.FUN = max)

col_cluster <- as.hclust(dend)


pheatmap(hallmark_result, cluster_cols = col_cluster)

pheatmap(hallmark_result, cluster_cols = col_cluster, show_colnames = T, scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

####pyscenic####
library(SCENIC)
library(SCopeLoomR)
loom <- open_loom('LUAD_Inter_SCLC-7.scenic.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
close_loom(loom)

sub_regulonAUC <- regulonAUC[,match(colnames(LUAD_Inter_SCLC), colnames(regulonAUC))]
dim(sub_regulonAUC)
identical(colnames(sub_regulonAUC),colnames(LUAD_Inter_SCLC))

cellClusters <- data.frame(row.names = colnames(LUAD_Inter_SCLC), 
                           Cluster = as.character(LUAD_Inter_SCLC$celltype))

selectedResolution <- "Cluster"
cellsPerGroup <- split(rownames(cellClusters), 
                       cellClusters[,selectedResolution])

sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T))
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

regulonActivity_byGroup_Scaled <- regulonActivity_byGroup_Scaled %>% as.data.frame %>% dplyr::select(-2, everything()) %>%
  as.matrix()

BMC_tumor_regulon <- rownames(regulonActivity_byGroup_Scaled)


colnames(regulonActivity_byGroup_Scaled) <- paste0("Cluster_", colnames(regulonActivity_byGroup_Scaled))
write.table(regulonActivity_byGroup_Scaled, "pyscenic_on_LUAD_Inter_SCLC.txt", sep = '\t', row.names = T)
library(circlize)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)
library(reshape2)
colnames(regulonActivity_byGroup_Scaled)[colnames(regulonActivity_byGroup_Scaled) == 'Cluster_LUAD'] <- 'LUAD'
colnames(regulonActivity_byGroup_Scaled)[colnames(regulonActivity_byGroup_Scaled) == 'Cluster_Intermediate'] <- 'Intermediate'
colnames(regulonActivity_byGroup_Scaled)[colnames(regulonActivity_byGroup_Scaled) == 'Cluster_SCLC'] <- 'SCLC'

regulonActivity_byGroup_Scaled_1 = regulonActivity_byGroup_Scaled[,c(3,1:2)]
Heatmap(
  regulonActivity_byGroup_Scaled_1,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)


save(LUAD_Inter_SCLC,file = "LUAD_Inter_SCLC.Rda")

AT2_LUAD_Inter_SCLC<-c("AT2","LUAD","Intermediate","SCLC")
AT2_LUAD_Inter_SCLC<-subset(epi,idents = AT2_LUAD_Inter_SCLC)
system.time({fwrite(x = as.data.frame(AT2_LUAD_Inter_SCLC[["RNA"]]@counts), row.names=T,file = "AT2_LUAD_Inter_SCLC-7-9-.csv")})


