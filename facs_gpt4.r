install.packages("openai")
remotes::install_github("Winnie09/GPTCelltype")
library("presto")

lapply(c("dplyr", "Seurat", "HGNChelper"), library, character.only = T)
library("openxlsx")
library("ggplot2")

facs <- LoadSeuratRds("C:/Users/mlalw/OneDrive/Desktop/Stat 35510/facs_tsned.rds")

markers <- FindAllMarkers(facs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

library(GPTCelltype)
library(openai)

prompt <- gptcelltype(markers, model = 'gpt-4')

res_names <- c(
  "Microglia", "Fibroblasts", "B cells", "Fibroblasts", "Endothelial cells", 
  "NA", "Epidermal keratinocytes", "Neutrophils", "Epidermal keratinocytes", 
  "Myeloid cells", "T cells", "Epithelial cells", "Smooth muscle cells", "Secretory cells", 
  "Endothelial cells", "Muscle cells", "Macrophages", "Control genes", "Secretory cells", 
  "Secretory cells", "Epidermal keratinocytes", "Dendritic cells", "Fibroblasts", 
  "Epidermal keratinocytes", "Myeloid cells", "T cells", "Endocrine cells", "T cells", 
  "B cells", "Proximal tubule cells", "Ribosomal protein expressing cells", "Liver cells", 
  "Digestive enzymes", "Pancreatic beta cells", "Epithelial cells", "Endothelial cells", 
  "B cells", "Smooth muscle cells", "Fibroblasts", "Epithelial cells", "Neutrophils", 
  "Red blood cells", "Macrophages", "Cardiac muscle cells", "Pancreatic acinar cells", 
  "T cells", "Alveolar cells", "Pancreatic acinar cells"
)

res <- setNames(res_names, 0:(length(res_names)-1))

facs@meta.data$celltype <- as.factor(res[as.character(Idents(facs))])
p <- DimPlot(facs, reduction = 'tsne', group.by='celltype', pt.size = 1, label = TRUE, label.size = 6) +
  theme(
  plot.title = element_text(size = 40),
  axis.title = element_text(size = 30),
  axis.text = element_text(size = 30),
  legend.text = element_text(size = 23)
  ) +
  ggtitle("GPT4 Cell Type Annotation of FACS Data")
ggsave("GP4_facs_annotated.png", plot = p, width = 24, height = 12, dpi = 300)