install.packages("openai")
remotes::install_github("Winnie09/GPTCelltype")
library("presto")

lapply(c("dplyr", "Seurat", "HGNChelper"), library, character.only = T)
library("openxlsx")
library("ggplot2")

pbmc <- LoadSeuratRds("C:/Users/mlalw/OneDrive/Desktop/Stat 35510/pbmc_tsned.rds")

markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

library(GPTCelltype)
library(openai)

prompt <- gptcelltype(markers, model = 'gpt-4')

res_names <- c('Ribosomal Cells', 'Myeloid Cells (Macrophages, Neutrophils)', 'T Cells', 'B Cells', 'Cytotoxic T Cells (CD8 + T Cells)', 'T Cells', 'Monocytes', 'Cytotoxic T Cells (CD8 + T Cells)', 'Dendritic Cells')
res <- setNames(res_names, 0:(length(res_names)-1)) # rerun this

pbmc@meta.data$celltype <- as.factor(res[as.character(Idents(pbmc))])
p <- DimPlot(pbmc, group.by='celltype', pt.size = 1) +
  theme(
  plot.title = element_text(size = 40),
  axis.title = element_text(size = 30),
  axis.text = element_text(size = 30),
  legend.text = element_text(size = 23)
  ) +
  ggtitle("GPT4 Cell Type Annotation of PBMC Data")
ggsave("GP4_pbmc_annotated.png", plot = p, width = 20, height = 12, dpi = 300)
