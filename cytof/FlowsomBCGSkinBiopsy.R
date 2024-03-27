
library(readxl)
library(CATALYST)
library(flowCore)
library(cowplot)
library(diffcyt)
library(here)
library(FlowSOM)
library(tidyverse)
library(ggpubr)

outdir <- file.path(here::here(), "output")

# Load fcs files (n=10 at Day 3, n=8 at Day 15)
fcs <- dir("/home/emmabishop/workspace/human-bcg-challenge-cytof/cytof/data/BCG Skin Biopsy 1 Live CD45+ cells FCS", full.names=T, pattern="*.fcs$") 

# Load metadata and panel
md <- read_excel("/home/emmabishop/workspace/human-bcg-challenge-cytof/cytof/data/2022_human_BCG_challenge_metadata_for_flowsom.xlsx")
panel <- read_excel("/home/emmabishop/workspace/human-bcg-challenge-cytof/cytof/data/2022_human_BCG_challenge_CyTOF_panel.xlsx")

# Load manual annotations
merging_table1 <- read_excel("/home/emmabishop/workspace/human-bcg-challenge-cytof/cytof/data/BCG_challenge_cluster_annotations_24.xlsx")

# Create flowSet object
fs <- read.flowSet(fcs, transformation = FALSE, truncate_max_range = FALSE)

# Check that all panel columns are in the flowSet object
all(panel$fcs_colname %in% colnames(fs)) # TRUE

# Specify levels for ordering
md$timepoint <- factor(md$timepoint, levels = c("Day_3", "Day_15"))
md$sample_id <- factor(md$sample_id, levels = md$sample_id[order(md$timepoint)])

# Create SingleCellExperiment object
## By default, this arcsinh transformers marker expresssions with a cofactor of 5.
sce <- prepData(fs, panel, md, features = panel$fcs_colname,
                panel_cols = list(channel = "fcs_colname", antigen = "antigen"),
                md_cols = list(file = "fcs_name", id = "sample_id", 
                               factors = c("patient_id", "timepoint", "group", "age", "sex")))

## Diagnostic plots

p <- plotExprs(sce, color_by = "timepoint")
p$facet$params$ncol <- 10
p

ncells <- n_cells(sce)
ncells

# write.csv(ncells, "C:/Users/Steven Makatsa/OneDrive - UW/Shared Documents - SeshadriLab/Members/YuKrystle/Projects_and_Data/BCG_Skin_Biopsies/Mass Cytometry/ProcessedData/Rerun FlowSOM_SM/out/BCG Skin Biopsy total CD45 count.csv", row.names=FALSE)

# Bar plot of total CD45+ count per sample
cellcount_bar <- plotCounts(sce, group_by = "sample_id", color_by = "timepoint") +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() + 
  scale_fill_manual(values =  c("Day_3" = "white", "Day_15" = "gray50")) +
  theme(text = element_text(family="Arial"),
        axis.title.x = element_text(size = 9.5),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                                   color = "black", size =10),
        axis.text.y = element_text(color="black", size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = c(0.9, 0.9)) +
  ylab("Total CD45+ Cells") +
  xlab("Sample ID")
cellcount_bar

# Save
ggsave(file.path(outdir, "plots/SuppFig3_cellcount_bar.png"), plot = cellcount_bar,
       dpi = 300, width = 5, height = 3, device = "png")

cairo_pdf(file = file.path(outdir, "plots/SuppFig3_cellcount_bar.pdf"), 
          width = 5, height = 3, bg = "transparent", family = "Arial")
print(cellcount_bar)
dev.off()


# Non-redundancy scores (NRS) across samples for all markers
nrs_boxplot <- plotNRS(sce, color_by = "timepoint") +
  scale_color_manual(labels = c("Day_3", "Day_15"), values = c("#E66101", "gray50")) +
  theme_classic() +
  theme(text = element_text(family="Arial"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                                   color = "black", size = 10),
        axis.text.y = element_text(color="black", size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = c(0.9, 0.9),
        legend.spacing.x = unit(0.05, 'cm')) +
  ylab("Non-Redundancy Score") +
  xlab("Marker")
nrs_boxplot

# Save
ggsave(file.path(outdir, "plots/SuppFig3_nrs_boxplot.png"), plot = nrs_boxplot,
       dpi = 300, width = 7, height = 4.2, device = "png")

cairo_pdf(file = file.path(outdir, "plots/SuppFig3_nrs_boxplot.pdf"), 
          width = 7, height = 4.2, bg = "transparent", family = "Arial")
print(nrs_boxplot)
dev.off()

# MDS plot
pbMDS(sce, color_by = "timepoint", label_by = "sample_id")

## Perform unsupervised clustering

set.seed(1234)
sce <- cluster(sce, features = "type",
               xdim = 12, ydim = 12, maxK = 24, seed = 1234)

# Heatmap before annotation
plotExprHeatmap(sce, features = "type", 
                      by = "cluster_id", k = "meta24", 
                      bars = TRUE, perc = TRUE)

## Visual representation with UMAP

# run UMAP
set.seed(1234)
sce <- runDR(sce, "UMAP", cells = 10000, features = "type")

plotDR(sce, "UMAP", color_by = "CD3")

p1 <- plotDR(sce, "UMAP", color_by = "CD68") + 
  theme(legend.position = "none")
p2 <- plotDR(sce, "UMAP", color_by = "CD163") + 
  theme(legend.position = "none")
p3 <- plotDR(sce, "UMAP", color_by = "CD14") + 
  theme(legend.position = "none")
p4 <- plotDR(sce, "UMAP", color_by = "CD11c") + 
  theme(legend.position = "none")
p5 <- plotDR(sce, "UMAP", color_by = "CD11b") + 
  theme(legend.position = "none")
p6 <- plotDR(sce, "UMAP", color_by = "HLA_DR") + 
  theme(legend.position = "none")
p7 <- plotDR(sce, "UMAP", color_by = "CD20") + 
  theme(legend.position = "none")
p8 <- plotDR(sce, "UMAP", color_by = "CD3") + 
  theme(legend.position = "none")
p9 <- plotDR(sce, "UMAP", color_by = "MR1") + 
  theme(legend.position = "none")
p10 <- plotDR(sce, "UMAP", color_by = "CD66abce") + 
  theme(legend.position = "none")

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2)

# facet by sample
plotDR(sce, "UMAP", color_by = "meta24")

# facet by timepoint 
plotDR(sce, "UMAP", color_by = "meta24", facet_by = "timepoint")

plotDR(sce, "UMAP", color_by = "meta24", facet_by = "sample_id")

## Apply manual merging

merging_table1 <- here::here("cytof/data/BCG_challenge_cluster_annotations_24.xlsx")
merging_table1 <- read_excel(merging_table1)
head(data.frame(merging_table1))

# convert to factor with merged clusters
merging_table1$new_cluster <- factor(merging_table1$new_cluster)

# apply manual merging
sce <- mergeClusters(sce, k = "meta24", 
                     table = merging_table1, id = "Annotations", overwrite = TRUE)

plotDR(sce, "UMAP", color_by = "Annotations")

# Heatmap labelled by cluster and our annotations
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta24", m = "Annotations")

# Heatmap labelled by annotation with dendrogram
plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "Annotations", 
                bars = TRUE, perc = TRUE)

# Final heatmap
cytof_heat_plot <- plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "Annotations", 
                bars = TRUE, perc = TRUE,
                row_dend = FALSE,
                col_dend = FALSE)
cytof_heat_plot

cairo_pdf(file = file.path(outdir, "plots/Fig4_cytof_heatmap.pdf"), 
          width=9, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial")
print(cytof_heat_plot)
dev.off()


plotDR(sce, "UMAP", color_by = "Annotations", facet_by = "timepoint")

plotDR(sce, "UMAP", color_by = "Annotations", facet_by = "sample_id")

# remove a BCG07_02 sample
filterSCE(sce, sample_id != "BCG07_2")

FDR_cutoff <- 0.05

## Differential analysis

plotAbundances(sce, k = "Annotations", by = "cluster_id", 
                      shape_by = "group", group_by = "timepoint")

ei <- metadata(sce)$experiment_info
(da_formula1 <- createFormula(ei, 
                              cols_fixed = "timepoint", 
                              cols_random = "sample_id"))

(da_formula2 <- createFormula(ei, 
                              cols_fixed = "timepoint", 
                              cols_random = c("sample_id", "patient_id")))

contrast <- createContrast(c(1, 1))

da_res1 <- diffcyt(sce, 
                   formula = da_formula1, contrast = contrast,
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "Annotations", verbose = FALSE)
da_res2 <- diffcyt(sce, 
                   formula = da_formula2, contrast = contrast,
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "Annotations", verbose = FALSE)

names(da_res1)

rowData(da_res1$res) 

table(rowData(da_res1$res)$p_adj < FDR_cutoff)

table(rowData(da_res2$res)$p_adj < FDR_cutoff)

topTable(da_res2, show_props = TRUE, format_vals = TRUE, digits = 2)

table <- topTable(da_res2, show_props = TRUE, format_vals = TRUE, digits = 2)

table <- as.data.frame(table)

# Reformat for use by other scripts
table_out <- table %>%
  select(-c(p_val, p_adj)) %>%
  # Transpose
  t() %>%
  data.frame() %>%
  rownames_to_column("cluster_id") %>%
  janitor::row_to_names(1) %>%
  mutate(cluster_id = gsub("props_", "", cluster_id),
         patient_id = str_sub(cluster_id, 1, 5),
         timepoint = case_when(
           grepl("_1", cluster_id) ~ "Day_3",
           .default = "Day_15"),
         .after = "cluster_id")

write_csv(table_out, file.path(outdir, "processed_data/BCG_Skin_Biopsy_CD45_Subsets.csv"))

plotDiffHeatmap(sce, rowData(da_res2$res), all = TRUE, fdr = FDR_cutoff)

