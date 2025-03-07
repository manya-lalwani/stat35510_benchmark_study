lapply(c("dplyr", "Seurat", "HGNChelper"), library, character.only = T)
library("openxlsx")
library("ggplot2")
library("presto")

# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')

pbmc <- LoadSeuratRds("C:/Users/mlalw/OneDrive/Desktop/Stat 35510/pbmc_tsned.rds")

markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cluster_gene_list <- markers %>%
  group_by(cluster) %>%
  summarise(genes = paste(gene, collapse = ", ")) %>%
  ungroup()

writeLines(paste0("cluster_", cluster_gene_list$cluster, ": ", cluster_gene_list$genes), 
           "pbmc_markers_by_cluster.txt")

res_names <- c('Central Memory CD4+ Alpha-Beta T Cell', 'Monocyte', 'Central Memory CD4+ Alpha-Beta T Cell', 'B cell', 'Effector Memory CD8+ Alpha-Beta T Cell', 'Unknown', 'Macrophage', 'Natural killer cell', 'Dendritic cell', 'Megakaryocyte')
res <- setNames(res_names, 0:(length(res_names)-1))
pbmc@meta.data$celltype <- as.factor(res[as.character(Idents(pbmc))])
p <- DimPlot(pbmc, group.by='celltype', pt.size = 1, label = TRUE, label.size = 9) +
  theme(
  plot.title = element_text(size = 40),
  axis.title = element_text(size = 30), 
  axis.text = element_text(size = 30),
  legend.text = element_text(size = 23)
  ) +
  ggtitle("ACT Cell Type Annotation of PBMC Data")
ggsave("ACT_pbmc_annotated.png", plot = p, width = 20, height = 12, dpi = 300)








# facs

facs <- LoadSeuratRds("C:/Users/mlalw/OneDrive/Desktop/Stat 35510/facs_tsned.rds")

markers <- FindAllMarkers(facs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cluster_gene_list <- markers %>%
  group_by(cluster) %>%
  summarise(genes = paste(gene, collapse = ", ")) %>%
  ungroup()                                                                                                  e 

writeLines(paste0("cluster_", cluster_gene_list$cluster, ": ", cluster_gene_list$genes), 
           "facs_markers_by_cluster.txt")




res_names <- c('Microglial', 'Alveolar type 2 fibroblast', 'Naive B Cell', 'Alveolar type 2 fibroblast', 'Endothelial', 'Unknown', 'Basal Cell', 'Hematopoietic Stem Cell', 'Hair Follicle', 'Kupffer Macrophage', 'Memory T Cell', 'Cholangiocyte Epithelial Cell', 'Vascular Smooth Muscle', 'Unknown', 'Endocardial', 'Alveolar type 2 fibroblast', 'Neutrophil', 'Unknown', 'Goblet epithelial', 'Proximal Tubule Brush Border Cell', 'Basal Cell', 'Langerhans Dendritic Cell', 'Parafollicular epithelial', 'Superficial epithelial', 'Microglial', 'Double-positive, alpha-beta thymocyte T cell', 'Pancreatic PP epithelial', 'Natural Killer', 'Immature B cell', 'Proximal Tubule Brush Border Cell', 'Stromal', 'Hepatocyte epithelial', 'Acinar epithelial', 'Enteroendocrine epithelial', 'Luminal cell', 'Fat cell', 'B cell', 'Mesangial', 'Chondrocyte stromal', 'Luminal cell', 'Neutrophil', 'Erythroblast', 'Microvascular endothelial', 'Cardiac muscle', 'Acinar epithelial', 'Double-positive, alpha-beta thymocyte T cell', 'Type II pneumocyte epithelial', 'Unknown', 'Acinar epithelial')
res <- setNames(res_names, 0:(length(res_names)-1))
facs@meta.data$celltype <- as.factor(res[as.character(Idents(facs))])
p <- DimPlot(facs, reduction = 'tsne', group.by='celltype', pt.size = 1, label = TRUE, label.size = 6) +
  theme(
  plot.title = element_text(size = 40),
  axis.title = element_text(size = 30), 
  axis.text = element_text(size = 30),
  legend.text = element_text(size = 23)
  ) +
  ggtitle("ACT Cell Type Annotation of FACS Data")
ggsave("ACT_facs_annotated.png", plot = p, width = 30, height = 16, dpi = 300)