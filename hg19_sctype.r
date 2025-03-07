

# load libraries
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
library("openxlsx")
library("ggplot2")
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:/Users/mlalw/OneDrive/Desktop/Stat 35510/hg19")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# normalize data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# scale and run PCA, visualize
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca", pt.size = 1.5) +
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.position = "none"
  ) +
  ggtitle("PCA Reduction of PBMC Data") +
  scale_x_continuous(limits = c(-20, 10)) +
  scale_y_continuous(limits = c(-20, 10))

# Check number of PC components (going with 11)
ElbowPlot(pbmc) +
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30)
  ) +
  ggtitle("PCA Elbow Plot of PBMC Data")

# cluster and visualize
pbmc <- FindNeighbors(pbmc, dims = 1:11)
pbmc <- FindClusters(pbmc, resolution = 0.8) # louvain
pbmc <- RunUMAP(pbmc, dims = 1:11)
DimPlot(pbmc, reduction = "umap", pt.size = 1.5) +
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 30)
  ) +
  ggtitle("UMAP Reduction of PBMC Data")

# sctype cell assignment

# load gene set preparation function and cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
gs_list <- gene_sets_prepare(db_, tissue)

seurat_package_v5 <- isFALSE('counts' %in% names(attributes(pbmc[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(pbmc[["RNA"]]$scale.data) else as.matrix(pbmc[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_results <- do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])


pbmc@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  pbmc@meta.data$sctype_classification[pbmc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

p <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification', label.size = 8)+
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 20)
  ) +
  ggtitle("Cell Type Annotations - PBMC Data")

ggsave("UMAP_plot.png", plot = p, width = 15, height = 10, dpi = 300) # might not be right anymore

# make plot with annotation choices
cL_results <- cL_results[order(cL_results$cluster),]; edges = cL_results; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

nodes_lvl1 <- sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss <- c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_results$cluster))){
  dt_tmp = cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes <- rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db <- openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

gggr <- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
  
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)+ gggr +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.text = element_text(size = 20)
  ) +
  ggtitle("Annotation Choices - PBMC Data")
ggsave("gggr_plot.png", plot = p, width = 25, height = 10, dpi = 300)
