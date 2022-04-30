rm(list = ls())

library(Seurat)
library(tidyverse)
library(SCPA)
library(msigdbr)
library(magrittr)
library(ComplexHeatmap)
library(SingleCellExperiment)
library(scGate)
library(tidySingleCellExperiment)
library(tidyseurat)
library(GGally)

library(AnnotationDbi)
library(annotables)
library(RColorBrewer)



# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# create working directory
workingDir <- "WorkingDirectory"
dirPath <- file.path(PrimaryDirectory, workingDir)
dir.create(dirPath)
setwd(dirPath)

saveDir <- file.path(dirPath, "results")
savePath <- saveDir
dir.create(saveDir)

filesPath <- file.path(dirPath, "files")
dir.create(filesPath)

plotDir <- file.path(saveDir, "plots")
dir.create(plotDir)

list.files(filesPath)

seurat_ext <- read_rds(file.path(filesPath, "GSE195937_Mouse_breast_ExtFig_scRNA_SeuratObj.rds"))
seurat_ext <- UpdateSeuratObject(seurat_ext)
seurat_ext

seurat <- read_rds(file.path(filesPath, "GSE195937_Mouse_breast_Fig1_scRNA_SeuratObj.rds"))
seurat <- UpdateSeuratObject(seurat)
seurat

###################### Plotting ##########################################
# Use colourblind-friendly colours
friendly_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC")

# Set theme
my_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
      )
  )

seurat %>%
  as_tibble() %>%
  dplyr::select(contains("PC"), everything()) %>%
  GGally::ggpairs(columns = 1:5, ggplot2::aes(color = Cluster)) +
  my_theme

# Identify top 10 markers per cluster
markers <-
  seurat %>%
  FindAllMarkers(only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)

# Plot heatmap
seurat %>%
  DoHeatmap(
    features = markers$gene,
    group.colors = friendly_cols
  )

seurat %>%
  ggplot(mapping = aes(UMAP_1, UMAP_2, color = Cluster)) +
  geom_point() +
  my_theme

  
levels(seurat)
C2_markers <- FindMarkers(seurat, ident.1 = "2", ident.2 = NULL)
C2_markers %>% 
  rownames_to_column(var = "feature") %>%
  as_tibble() %>%
  dplyr::inner_join(., grcm38, by = c("feature" = "ensgene")) %>%
  dplyr::select(symbol, feature:p_val_adj) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::select(symbol) %>%
  write_csv(file = file.path(dirPath, "abILTCK_signature.csv"))

TPEX_BM_sig <- read_csv(file = file.path(dirPath, "TPEX_BM_signature.csv")) %>%
  dplyr::inner_join(., grcm38, by = c("feature" = "symbol")) %>%
  dplyr::select(feature, ensgene) %>%
  dplyr::filter(ensgene %in% rownames(seurat))

grep("Fcer1g", TPEX_BM_sig$feature)
TPEX_BM_sig[grep("Fcer1g", TPEX_BM_sig$feature),]
TPEX_BM_sig[grep("Klrb1c", TPEX_BM_sig$feature),]
TPEX_BM_sig[grep("Cd7", TPEX_BM_sig$feature),]

# TPEX_BM_sig <- c("ENSMUSG00000058715", "ENSMUSG00000030325", "ENSMUSG00000025163")

TPEX_BM_sig <- TPEX_BM_sig %>%
  head(40) %>%
  pull(ensgene)

BM_sig <- read_csv(file = file.path(dirPath, "BM_sig.csv")) %>%
  dplyr::inner_join(., grcm38, by = c("feature" = "symbol")) %>%
  dplyr::select(feature, ensgene) %>%
  dplyr::filter(ensgene %in% rownames(seurat))
BM_symbol <- BM_sig$feature

grep("Klrb1c", BM_sig$feature)
grep("Cd7", BM_sig$feature)
grep("", BM_sig$feature)

BM_sig <- BM_sig %>% pull(ensgene)

rownames(seurat@assays$RNA) %>% head(10)
seurat <- seurat %>% AddModuleScore(., 
                                    features = list(BM_sig),
                                    name = "BM")
seurat
FeaturePlot(seurat, features = "BM1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(-0.2, 0.4), label = c("-0.2", "Maximum"))

abILTCK <- read_csv(file = file.path(dirPath, "abILTCK_signature.csv"))

BM_symbol %>% as_tibble() %>% rename(., symbol = value) %>% inner_join(., abILTCK, by = c("symbol" = "symbol"))
