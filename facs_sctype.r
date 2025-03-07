lapply(c("dplyr", "Seurat", "HGNChelper"), library, character.only = T)
library("openxlsx")
library("ggplot2")

facs_matrix_list <- list()

for (filename in c('Aorta', 'Bladder', 'Brain_Myeloid', 'Diaphragm', 'Fat', 'Heart', 'Kidney', 'Large_Intestine', 'Limb_Muscle', 'Liver', 'Lung', 'Mammary_Gland', 'Marrow', 'Pancreas', 'Skin', 'Spleen', 'Thymus', 'Tongue', 'Trachea')) {
  matrix <- read.csv(paste0("C:/Users/mlalw/OneDrive/Desktop/Stat 35510/FACS/", filename, "-counts.csv"), row.names = 1)
  facs_matrix_list[[filename]] <- matrix
}
facs_matrix <- do.call(cbind, facs_matrix_list)
facs <- CreateSeuratObject(counts = facs_matrix, project = "facs", min.cells = 3, min.features = 200)

facs[["percent.mt"]] <- PercentageFeatureSet(facs, pattern = "^MT-")
facs <- NormalizeData(facs, normalization.method = "LogNormalize", scale.factor = 10000)
facs <- FindVariableFeatures(facs, selection.method = "vst", nfeatures = 2000)

# scale and run PCA, visualize
facs <- ScaleData(facs, features = rownames(facs))
facs <- RunPCA(facs, features = VariableFeatures(object = facs)) # saved after here facs_pcad.rds
DimPlot(facs, reduction = "pca", pt.size = 1.5) +
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.position = "none"
  ) +
  ggtitle("PCA Reduction of FACS Data") +
  scale_x_continuous(limits = c(-20, 10)) +
  scale_y_continuous(limits = c(-20, 10))

# Check number of PC components (going with 19)
ElbowPlot(facs) +
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30)
  ) +
  ggtitle("PCA Elbow Plot of FACS Data")

# cluster and visualize
facs <- FindNeighbors(facs, dims = 1:19)
facs <- FindClusters(facs, resolution = 0.8) # louvain
facs <- RunTSNE(facs, dims = 1:19) # saving here facs_tsned.rds
p <- DimPlot(facs, reduction = "tsne", pt.size = 1.5) +
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 30)
  ) +
  ggtitle("TSNE Reduction of FACS Data")

ggsave("TSNE_reduction_facs.png", plot = p, width = 20, height = 20, dpi = 300)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";

gs_list <- list(
  gs_positive = list(),
  gs_negative = list()
)
# for (tissue in c('Heart', 'Brain', 'Muscle', 'Kidney', 'Intestine', 'Placenta', 'Immune system', 'Pancreas', 'Spleen', 'Thymus', 'Lung', 'Liver', 'Stomach')) {
#   gs <- gene_sets_prepare(db_, tissue)

#   if ('gs_positive' %in% names(gs)) {
#     gs_list$gs_positive <- c(gs_list$gs_positive, gs$gs_positive)
#   }

#   if ('gs_negative' %in% names(gs)) {
#     gs_list$gs_negative <- c(gs_list$gs_negative, gs$gs_negative)
#   }
# }

for (tissue in c('Heart', 'Brain', 'Muscle', 'Kidney', 'Intestine', 'Placenta', 'Immune system', 'Pancreas', 'Spleen', 'Thymus', 'Lung', 'Liver', 'Stomach')) {
  
  # Get the gene sets for the current tissue
  gs <- gene_sets_prepare(db_, tissue)
  
  # If gs_positive exists in the current tissue
  if ('gs_positive' %in% names(gs)) {
    # Append the positive gene sets to the gs_list list (without duplicating)
    for (cell_type in names(gs$gs_positive)) {
      if (!(cell_type %in% names(gs_list$gs_positive))) {
        gs_list$gs_positive[[cell_type]] <- gs$gs_positive[[cell_type]]  # Add the first time
      } else {
        gs_list$gs_positive[[cell_type]] <- unique(c(gs_list$gs_positive[[cell_type]], gs$gs_positive[[cell_type]]))  # Append without duplication
      }
    }
  }
  if ('gs_negative' %in% names(gs)) {
    # Append the negative gene sets to the gs_list list (without duplicating)
    for (cell_type in names(gs$gs_negative)) {
      if (!(cell_type %in% names(gs_list$gs_negative))) {
        gs_list$gs_negative[[cell_type]] <- gs$gs_negative[[cell_type]]  # Add the first time
      } else {
        gs_list$gs_negative[[cell_type]] <- unique(c(gs_list$gs_negative[[cell_type]], gs$gs_negative[[cell_type]]))  # Append without duplication
      }
    }
  }
}


seurat_package_v5 <- isFALSE('counts' %in% names(attributes(facs[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(facs[["RNA"]]$scale.data) else as.matrix(facs[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_results <- do.call("rbind", lapply(unique(facs@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(facs@meta.data[facs@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(facs@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3], n=49)


facs@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  facs@meta.data$sctype_classification[facs@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

p <- DimPlot(facs, reduction = "tsne", label = TRUE, repel = TRUE, group.by = 'sctype_classification', label.size = 9, pt.size = 0.4)+
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 20)
  ) +
  ggtitle("Cell Type Annotations - FACS Data")

ggsave("Annotated_facs.png", plot = p, width = 30, height = 20, dpi = 300)
