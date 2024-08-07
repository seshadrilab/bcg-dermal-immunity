# Author: Emma Bishop

library(tidyverse)
library(readxl)
library(ggpubr)
library(rstatix)
library(ggbeeswarm)
library(scales)
library(harmony)
library(cowplot)
library(ComplexHeatmap)

set.seed(4)

###############
## Load data ##
###############

script_output_dir <- file.path(here::here(), "output")
outdir <- file.path(here::here(), "output/plots")

# Create folders for data if needed
if(!dir.exists(file.path(script_output_dir))) {
  cat(sprintf("Creating folder %s\n", file.path(script_output_dir, "processed_data")))
  dir.create(file.path(script_output_dir, "processed_data"), recursive = T)
  cat(sprintf("Creating folder %s\n", file.path(script_output_dir, "plots")))
  dir.create(file.path(script_output_dir, "plots"), recursive = T)
}

outdir <- file.path(here::here(), "output/plots")

# Heatmap data from Andrew Fiore-Gartland
heat_dat <- read_csv("/media/emmabishop/5TBSharedStorage/aws_bam_fastq/2022_BCGChallenge/Round1/deg_heatmap_values.csv")

# All cells
day3_clstr_filt <- readRDS(file.path(script_output_dir, "processed_data/9_final_annot_d3.rds"))
day15_clstr_filt <- readRDS(file.path(script_output_dir, "processed_data/9_final_annot_d15.rds"))

# Immune cells
d3_subset <- readRDS(file.path(script_output_dir, "processed_data/8_annot_sub_final_d3.rds"))
d15_subset <- readRDS(file.path(script_output_dir, "processed_data/8_annot_sub_final_d15.rds"))

# Micro data
micro_f <- "/media/emmabishop/5TBSharedStorage/2022_BCGChallenge/Round1/micro/Combined CFU MVT RS data KER.xlsx"

# CyTOF cell population frequencies
all_full_count <- read_csv(file.path(script_output_dir, "processed_data/BCG_Skin_Biopsy_CD45_Subsets.csv"))

# Colors
ptid_colors <- c("BCG01" = "salmon", "BCG05" = "orange", "BCG07" = "#7570B3",
                 "BCG08" = "#0CB702", "BCG09" = "#A6761D", "BCG10" = "#00BFC4",
                 "BCG11" = "#C77CFF", "BCG13" = "darkgreen",
                 "BCG16" = "darkblue")

ptid_colors2 <- c("1" = "salmon", "5" = "orange", "7" = "#7570B3",
                  "8" = "#0CB702", "9" = "#A6761D", "10" = "#00BFC4",
                  "11" = "#C77CFF", "12" = "#FF61CC", "13" = "darkgreen",
                  "16" = "darkblue")

############################
## Pre-process micro data ##
############################

micro <- read_excel(micro_f,
                    sheet = 2,
                    skip = 3,
                    n_max = 9,
                    na = c("N/A", "Detectable pre-rRNA in Day-5 and -8 culture."),
                    col_names = c("Treatment", "PTID", "CFU_D3", 
                                  NA, NA, NA, NA, 
                                  "MVT_D3", NA, "RS_D3", "CFU_D15", 
                                  NA, NA, NA, NA, 
                                  "MVT_D15", NA, "RS_D15")) %>%
  mutate(MVT_D3 = gsub("NDT", 0, MVT_D3),
         MVT_D15 = gsub("NDT", 0, MVT_D15)) %>%
  select(Treatment, PTID, CFU_D3, MVT_D3, RS_D3, CFU_D15, MVT_D15, RS_D15) %>%
  pivot_longer(cols = -c(Treatment, PTID), 
               names_to = c(".value", "Day"), 
               names_pattern = "(\\w+)_(\\w+)") %>%
  mutate(MVT = as.numeric(MVT)) %>%
  mutate(Day = str_replace_all(Day, c("D3"="Day 3", "D15"="Day 15"))) %>%
  mutate(Day = factor(Day, levels = c("Day 3", "Day 15")),
         Treatment = factor(Treatment, levels = c("Non-INH", "INH")))


# Prep data
cfu_log <- micro %>%
  # Replace 0 with 1 for log scale
  mutate(CFU = gsub(0, 1, CFU)) %>%
  mutate(CFU = as.numeric(CFU))


###############################################
# Correlations between micro and cytof data? ##
###############################################

# Prep data

micro2 <- micro %>%
  mutate(Day = gsub(" ", "_", Day)) %>%
  unite("sample", PTID:Day, remove = F) %>%
  # Exclude RS data for this since only one PTID has RS data at both timepoints
  select(-c(RS))
head(micro2)

cytof2 <- all_full_count %>%
  unite("sample", patient_id:timepoint, remove = F)
head(cytof2)

micro_cytof <- left_join(cytof2, micro2, by = "sample") %>%
  select(-c(cluster_id, Treatment, PTID, Day, sample)) %>%
  drop_na()
head(micro_cytof)

# From long to wide
micro_cytof_wide <- pivot_wider(micro_cytof, 
                                id_cols = patient_id, 
                                names_from = timepoint, 
                                values_from = `B cells`:MVT)
head(micro_cytof_wide)

wide_cols <- colnames(micro_cytof_wide)

d15_colnames <- wide_cols[grepl("Day_15", wide_cols)]
d3_colnames <- wide_cols[grepl("Day_3", wide_cols)]

# Do subtraction
changes <- micro_cytof_wide  %>%
  mutate(across(all_of(d15_colnames), .names = "Change_{.col}") - across(all_of(d3_colnames))) %>%
  dplyr::rename(PTID = patient_id) %>%
  select(PTID, starts_with("Change")) %>%
  drop_na()

# Clean up column names
colnames(changes) <- gsub("_Day_15", "", colnames(changes))
colnames(changes) <- gsub(" ", "_", colnames(changes))
colnames(changes) <- gsub("\\+", "", colnames(changes))

# Do correlation for every combination of changes (in B cells, CFU, etc)
to_corr <- changes %>%
  select(starts_with("Change")) %>%
  drop_na()

# Do correlations for everything, keeping those with some hint of positive or 
# negative correlation (abs(r) > 0.5)
correlations <- to_corr %>% 
  corrr::correlate() %>%
  corrr::stretch() %>%
  filter(abs(r) > 0.5) %>%
  filter(x %in% c("Change_CFU", "Change_RS", "Change_MVT"))

# Manually inspect correlations, sorting by greatest abs(R)
correlations %>%
  arrange(desc(abs(r)))

# Test which with an abs R >=0.6 are significant
cor.test(to_corr$Change_MVT, to_corr$Change_Tcm, method = "pearson")  # NOT sig
cor.test(to_corr$Change_MVT, to_corr$Change_CD11c_T_cells, method = "pearson")  # NOT sig

###############################
## CFU (on solid 7H10 media) ##
###############################

colors <- c("Non-INH" = "#c5c9c7", "INH" = "#1d5dec")

cfu_all <- micro %>%
  mutate(CFU = as.numeric(CFU)) %>%
  select(PTID, Day, Treatment, CFU) %>%
  mutate(cfu_new = CFU/1000)

cfu_no_outlier <- cfu_all %>%
  filter(CFU != 92800)

# Do stats (not significant with or without outlier)
stat_cfu <- cfu_all %>%
  group_by(Day) %>%
  wilcox_test(cfu_new ~ Treatment, paired = FALSE) %>%
  add_significance() %>%
  add_xy_position(x = "Treatment")
stat_cfu

stat_cfu_no_outlier <-cfu_no_outlier %>%
  group_by(Day) %>%
  wilcox_test(cfu_new ~ Treatment, paired = FALSE) %>%
  add_significance() %>%
  add_xy_position(x = "Treatment")
stat_cfu_no_outlier

stat_cfu_d3 <- stat_cfu_no_outlier %>%
  filter(Day == "Day 3")

# Make plots with linear scale
make_cfu_plot <- function(df, usestat) {
  out_cfu <- ggplot(df, aes(x = Treatment, y = cfu_new, group = Treatment)) +
    geom_beeswarm(aes(fill = Treatment), size = 4, shape = 21, cex = 3) +  # Use shape for outline appearance
    scale_fill_manual(values = colors) +
    theme_classic() +
    ylab(expression(paste("CFU x ", 10^{3}))) +
    theme(text = element_text(family="Arial"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14.5),
          axis.text.x = element_text(color="black", size=13.5),
          axis.text.y = element_text(color="black", size=14.5),
          legend.position = "none") +
    facet_grid(. ~ Day) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
    theme(strip.text.x = element_text(size = 14.5)) +
    stat_pvalue_manual(usestat, label.size = 3.5, bracket.size = 0.2)
  return(out_cfu)
}

out_cfu_all <- make_cfu_plot(cfu_all, stat_cfu)
out_cfu_all


# Make with log10 scale

cfu_log <- micro %>%
  mutate(CFU = as.numeric(CFU)) %>%
  mutate(cfu_logaxis = case_when(
    CFU == 0 ~ 1,
    .default = CFU)) %>%
  mutate(cfu_logaxis = log10(cfu_logaxis))

stat_log <- cfu_log %>%
  group_by(Day) %>%
  wilcox_test(cfu_logaxis ~ Treatment, paired = FALSE) %>%
  add_significance() %>%
  add_xy_position(x = "Treatment")
stat_log

out_cfu_log <- ggplot(cfu_log, aes(x = Treatment, y = cfu_logaxis, group = Treatment)) +
  geom_beeswarm(aes(fill = Treatment), size = 4, shape = 21, cex = 3) +  # Use shape for outline appearance
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylab(expression(paste(log[10], " CFU"))) +
  theme(text = element_text(family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14.5),
        axis.text.x = element_text(color="black", size=13.5),
        axis.text.y = element_text(color="black", size=14.5),
        legend.position = "none") +
  facet_grid(. ~ Day) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
  theme(strip.text.x = element_text(size = 14.5)) +
  stat_pvalue_manual(stat_log, label.size = 3.5, bracket.size = 0.2)
out_cfu_log

# Save
ggsave(file.path(outdir, "Fig1_dotplot_cfu_log10.png"), plot = out_cfu_log,
       dpi = 300, width = 3.75, height = 3.25, device = "png")

cairo_pdf(file = file.path(outdir, "Fig1_dotplot_cfu_log10.pdf"), 
          width=3.75, height=3.25, bg = "transparent", family = "Arial")
print(out_cfu_log)
dev.off()


##############################################
## pre-rRNA:rDNA ratios (aka MVT_max_ratio) ##
##############################################

# Do stats
stat_mvt <- micro %>%
  group_by(Day) %>%
  wilcox_test(MVT ~ Treatment, paired = FALSE) %>%
  add_significance() %>%
  add_xy_position(x = "Treatment")
stat_mvt

# Make dot plot
out_mvt <- ggplot(micro, aes(x = Treatment, y = MVT, group = Treatment)) +
  geom_beeswarm(aes(fill = Treatment), size = 4, shape = 21, cex = 3) +  # Use shape for outline appearance
  stat_summary(fun.y = "mean", geom = "crossbar", size = 0.5,
               mapping = aes(ymin = after_stat(y), ymax = after_stat(y)), width = 0.5,
               position = position_dodge(), show.legend = FALSE) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylim(-1, 150) +
  ylab('Maximum pre-rRNA:rDNA ratio') +
  theme(text = element_text(family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14.5),
        axis.text.x = element_text(color="black", size=13.5),
        axis.text.y = element_text(color="black", size=14.5),
        legend.position = "none") +
  facet_grid(. ~ Day) +
  theme(strip.text.x = element_text(size = 14.5)) +
  stat_pvalue_manual(stat_mvt, label.size = 3.5, bracket.size = 0.2)
out_mvt


## Save plots ##
ggsave(file.path(outdir, "Fig1_dotplot_mvt.png"), plot = out_mvt,
       dpi = 300, width = 3.75, height = 3.25, device = "png")

cairo_pdf(file = file.path(outdir, "Fig1_dotplot_mvt.pdf"), 
          width=3.75, height=3.25, bg = "transparent", family = "Arial")
print(out_mvt)
dev.off()

###############
## RS ratios ##
###############

# Do stats
stat_rs <- micro %>%
  group_by(Day) %>%
  wilcox_test(RS ~ Treatment, paired = FALSE) %>%
  add_significance() %>%
  add_xy_position(x = "Treatment")
stat_rs

# Make dot plot
out_rs <- ggplot(micro, aes(x = Treatment, y = RS, group = Treatment)) +
  geom_beeswarm(aes(fill = Treatment), size = 4, shape = 21, cex = 3) +  # Use shape for outline appearance
  stat_summary(fun.y = "mean", geom = "crossbar", size = 0.5,
               mapping = aes(ymin = after_stat(y), ymax = after_stat(y)), width = 0.5,
               position = position_dodge(), show.legend = FALSE) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylab('R:S ratio') +
  ylim(0, 310) +
  theme(text = element_text(family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14.5),
        axis.text.x = element_text(color="black", size=13.5),
        axis.text.y = element_text(color="black", size=14.5),
        legend.position = "none") +
  facet_grid(. ~ Day) +
  theme(strip.text.x = element_text(size = 14.5)) +
  stat_pvalue_manual(stat_rs, label.size = 3.5, bracket.size = 0.2)
out_rs

## Save plots ##
ggsave(file.path(outdir, "Fig1_dotplot_rsratio.png"), plot = out_rs,
       dpi = 300, width = 3.75, height = 3.25, device = "png")

cairo_pdf(file = file.path(outdir, "Fig1_dotplot_rsratio.pdf"), 
          width=3.75, height=3.25, bg = "transparent", family = "Arial")
print(out_rs)
dev.off()


##################################
## Correlation rRNA:rDNA vs CFU ##
##################################

# Kendall tau correlation
# Using Kendall tau per E. Chandler Church
cor.test(micro$MVT, micro$CFU, method=c("kendall"))

micro2 <- micro %>%
  mutate(cfu_log = case_when(
    CFU == 0 ~ 1,
    .default = CFU))

corr_mvt_cfu <- ggscatter(micro2, x = "MVT", y = "cfu_log", add = "reg.line") +
  geom_point(aes(fill = Treatment), size = 4, shape = 21, cex = 3) +
  scale_fill_manual(values = colors) +
  stat_cor(method = "kendall", cor.coef.name = "tau", size = 4) +
  xlim(0, 150) +
  xlab('Maximum pre-rRNA:rDNA ratio') +
  ylab('CFU') +
  scale_y_log10(limits = c(1,1500000),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(text = element_text(family="Arial"),
        axis.title.x = element_text(size = 14.5),
        axis.title.y = element_text(size = 14.5),
        axis.text.x = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", size=14),
        legend.position = "none")
corr_mvt_cfu

## Save plots ##
ggsave(file.path(outdir, "SuppFig1_corr_mvt_cfu.png"), plot = corr_mvt_cfu,
       dpi = 300, width = 4, height = 3, device = "png")

cairo_pdf(file = file.path(outdir, "SuppFig1_corr_mvt_cfu.pdf"), 
          width = 4, height = 3, bg = "transparent", family = "Arial")
print(corr_mvt_cfu)
dev.off()



##############################
## Correlation plots legend ##
##############################

corr_legend_plt <- ggscatter(micro2, x = "RS", y = "MVT", add = "reg.line") +
  geom_point(aes(fill = Treatment), size = 6, shape = 21, cex = 3) +
  scale_fill_manual(values = colors) +
  stat_cor(method = "kendall", cor.coef.name = "tau", size = 6) +
  xlim(0, 310) +
  xlab('R:S ratio') +
  ylim(0, 150) +
  ylab('Maximum pre-rRNA:rDNA ratio') +
  theme_classic() +
  guides(fill=guide_legend(title="Arm", byrow = TRUE)) +
  theme(text = element_text(family="Arial"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(color="black", size=16),
        axis.text.y = element_text(color="black", size=16),
        legend.spacing.y = unit(0.2, 'cm'),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16))
corr_legend_plt

# Save
ggsave(file.path(outdir, "SuppFig1_corr_legend.png"), plot = corr_legend_plt,
       dpi = 300, width = 4, height = 3, device = "png")

cairo_pdf(file = file.path(outdir, "SuppFig1_corr_legend.pdf"), 
          width = 4, height = 3, bg = "transparent", family = "Arial")
print(corr_legend_plt)
dev.off()

##########################
##  All cell type UMAP  ##
##########################

# Combine into one object
d3_d15_combined <- merge(day3_clstr_filt, y = day15_clstr_filt)
VariableFeatures(d3_d15_combined[["SCT"]]) <- rownames(d3_d15_combined[["SCT"]]@scale.data)
d3_d15_combined <- RunPCA(d3_d15_combined, verbose = FALSE)

# Integrate with harmony
d3_d15_integ <- RunHarmony(d3_d15_combined, "orig.ident")
d3_d15_integ <- RunUMAP(d3_d15_integ, reduction = "harmony", dims = 1:20)

d3_d15_integ@meta.data <- d3_d15_integ@meta.data %>%
  mutate(PTID = factor(PTID, levels = c("1", "5", "7", "8", "9", "10", "11", "12", "13", "16")),
         orig.ident = factor(orig.ident, levels = c("Day3", "Day15")))

# UMAP of high-level cell types
clstr_colors <- c("Lymphoid" = "#00A9FF", "Myeloid" = "#F8766D",
                  "Epithelial" = "#0CB702", "Connective tissue" = "#00BE67",
                  "Nervous tissue" = "#A52A2A", "Endothelial" = "#E68613",
                  "Muscle tissue" = "#C77CFF", "Keratinocytes" = "#FF61CC",
                  "Unknown" = "gray")

# Just have axis lines in the corner
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1.75, "cm")
)

# Do cells separate by time point in UMAP space?
umap_time <- DimPlot(d3_d15_integ, group.by = "orig.ident",
        pt.size = 0.005)+
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),
                                                     ends = "last")),
        axis.title = element_text(hjust = 0.01),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin = grid::unit(c(0,0,0,0), 'cm')) +
  scale_x_discrete("UMAP1") +
  scale_y_discrete("UMAP2") +
  ggtitle("Timepoint")
umap_time

# Do cells separate by PTID in UMAP space?
umap_ptid <- DimPlot(d3_d15_integ, group.by = "PTID",
                     pt.size = 0.005)+
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),
                                                     ends = "last")),
        axis.title = element_text(hjust = 0.01),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin = grid::unit(c(0,0,0,0), 'cm')) +
  scale_color_manual(values = ptid_colors2,
                     labels = c("BCG01", "BCG05", "BCG07", "BCG08", "BCG09",
                                "BCG10", "BCG11", "BCG12", "BCG13", "BCG16")) +
  scale_x_discrete("UMAP1") +
  scale_y_discrete("UMAP2")
umap_ptid

# Do cells separate by ARM in UMAP space?
umap_arm <- DimPlot(d3_d15_integ, group.by = "ARM",
                     pt.size = 0.005)+
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),
                                                     ends = "last")),
        axis.title = element_text(hjust = 0.01),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin = grid::unit(c(0,0,0,0), 'cm')) +
  scale_color_manual(values = c("Non-INH" = "red", "INH" = "blue")) +
  scale_x_discrete("UMAP1") +
  scale_y_discrete("UMAP2")
umap_arm

# Color UMAP by high-level annotations
all_umap <- DimPlot(d3_d15_integ,
                    pt.size = 0.005) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),
                                                     ends = "last")),
        axis.title = element_text(hjust = 0.01),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin = grid::unit(c(0,0,0,0), 'cm')) +
  scale_color_manual(values = clstr_colors) +
  scale_x_discrete("UMAP1") +
  scale_y_discrete("UMAP2") +
  NoLegend()

all_umap <- LabelClusters(all_umap, id = "ident", size = 4, repel = T, box.padding = 0.7)
all_umap

# Save PNGs
ggsave(file.path(outdir, "Fig5_umap_all.png"), plot = all_umap,
       dpi = 300, width = 4, height = 4, device = "png")
ggsave(file.path(outdir, "SuppFig4_umap_time.png"), plot = umap_time,
       dpi = 300, width = 4, height = 4, device = "png")
ggsave(file.path(outdir, "SuppFig4_umap_ptid.png"), plot = umap_ptid,
       dpi = 300, width = 4, height = 4, device = "png")
ggsave(file.path(outdir, "SuppFig4_umap_arm.png"), plot = umap_arm,
       dpi = 300, width = 4, height = 4, device = "png")


# Save PDFs
cairo_pdf(file = file.path(outdir, "Fig5_umap_all.pdf"), 
          width=5, height=5, bg = "transparent", family = "Arial")
print(all_umap)
dev.off()

cairo_pdf(file = file.path(outdir, "SuppFig4_umap_time.pdf"), 
          width=5, height=5, bg = "transparent", family = "Arial")
print(umap_time)
dev.off()

cairo_pdf(file = file.path(outdir, "SuppFig4_umap_ptid.pdf"), 
          width=5, height=5, bg = "transparent", family = "Arial")
print(umap_ptid)
dev.off()

cairo_pdf(file = file.path(outdir, "SuppFig4_umap_arm.pdf"), 
          width=5, height=5, bg = "transparent", family = "Arial")
print(umap_arm)
dev.off()


###############################
## All cell type FeaturePlot ##
###############################

make_marker_plot <- function(obj, marker) {
  out <- FeaturePlot(obj, features = marker, pt.size = 0.2) +
    theme(legend.position = "none",
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, 
                                    size = 20,
                                    face = "italic"))
  return(out)
}

all_1 <- make_marker_plot(d3_d15_integ, "IL7R")
all_2 <- make_marker_plot(d3_d15_integ, "LYZ")
all_3 <- make_marker_plot(d3_d15_integ, "S100A2")
all_4 <- make_marker_plot(d3_d15_integ, "DCN")
all_5 <- make_marker_plot(d3_d15_integ, "DCT")
all_6 <- make_marker_plot(d3_d15_integ, "FLT1")
all_7 <- make_marker_plot(d3_d15_integ, "COL3A1")
all_8 <- make_marker_plot(d3_d15_integ, "DMKN")

all_grid <- plot_grid(all_1, all_2, all_3, all_4, all_5, all_6, all_7, all_8, nrow = 2)
all_grid

# Save
ggsave(file.path(outdir, "SuppFig4_featureplot_all.png"), plot = all_grid,
       dpi = 300, width = 6, height = 4, device = "png")

# Save PDFs
cairo_pdf(file = file.path(outdir, "SuppFig4_featureplot_all.pdf"), 
          width=6, height=4, bg = "transparent", family = "Arial")
print(all_grid)
dev.off()

#############################
##  Immune cell type UMAP  ##
#############################

# Combine into one object
imm_combined <- merge(d3_subset, y = d15_subset)
VariableFeatures(imm_combined[["SCT"]]) <- rownames(imm_combined[["SCT"]]@scale.data)
imm_combined <- RunPCA(imm_combined, verbose = FALSE)

# Integrate with harmony
imm_integ <- RunHarmony(imm_combined, "orig.ident")
imm_integ <- RunUMAP(imm_integ, reduction = "harmony", dims = 1:20)

# Do cells seprate in UMAP space by timepoint, PTID, SEX, ARM, or AGE?
# DimPlot(imm_integ, group.by = "orig.ident")
# DimPlot(imm_integ, group.by = "PTID")
# DimPlot(imm_integ, group.by = "SEX")
# DimPlot(imm_integ, group.by = "ARM")
# FeaturePlot(imm_integ, "AGE")

imm_integ@meta.data <- imm_integ@meta.data %>%
  mutate(eb_idents = case_when(
    eb_idents %in% c("B intermediate", "B memory") ~ "B cells",
    eb_idents %in% c("NK Proliferating", "NK_CD56bright") ~ "NK", 
    .default = eb_idents
  ))

imm_integ@meta.data$eb_idents <- factor(imm_integ@meta.data$eb_idents,
                                        levels = c("CD4 TCM", "CD4 TEM", "CD4 Proliferating",
                                                   "CD8 TCM", "CD8 TEM", "CD8 Proliferating", "CD8 Naive",
                                                   "Treg", "dnT", "gdT", "MAIT", "CD14 Mono", "CD16 Mono",
                                                   "cDC1", "cDC2", "pDC","ILC", "NK", "Platelet", "B cells",
                                                   "HSPC"))
subtype_cols <- c(
  "B cells" = "darkorange",
  "cDC1" = "red",
  "cDC2" = "red",
  "pDC" = "red",
  "ILC" = "darkgreen",
  "CD14 Mono" = "lightsalmon",
  "CD16 Mono" = "lightsalmon",
  "NK" = "#0CB702",
  "Platelet" = "gold",
  "CD8 TCM" = "#00BFC4",
  "CD8 Naive" = "#00BFC4",
  "CD8 Proliferating" = "#00BFC4",
  "CD8 TEM" = "#00BFC4",
  "CD4 Proliferating" = "#8ab8fe",
  "CD4 TEM" = "#8ab8fe",
  "CD4 TCM" = "#8ab8fe",
  "Treg" = "cyan",
  "dnT" = "cyan",
  "gdT" = "cyan",
  "MAIT" = "cyan",
  "HSPC" = "navy")

Idents(imm_integ) <- "eb_idents"

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1.75, "cm")
)

umap_immune <- DimPlot(imm_integ,
                    pt.size = 0.005, label = FALSE) +
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),
                                                     ends = "last")),
        axis.title = element_text(hjust = 0.01),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin = grid::unit(c(0,0,0,0), 'cm')) +
  scale_color_manual(values = subtype_cols) +
  scale_x_discrete("UMAP1") +
  scale_y_discrete("UMAP2") +
  NoLegend()
umap_immune

# With legend
umap_immune_legend <- DimPlot(imm_integ,
                       pt.size = 0.005, label = FALSE) +
  guides(x = axis, y = axis,
         color = guide_legend(override.aes = list(size=3), ncol=1)) +
  theme(axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),
                                                     ends = "last")),
        axis.title = element_text(hjust = 0.01),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin = grid::unit(c(0,0,0,0), 'cm')) +
  scale_color_manual(values = subtype_cols) +
  scale_x_discrete("UMAP1") +
  scale_y_discrete("UMAP2")
umap_immune_legend

imm_legend <- ggpubr::get_legend(umap_immune_legend)
imm_legend2 <- as_ggplot(imm_legend)

# Save PNG
ggsave(file.path(outdir, "Fig5_umap_immune_revised.png"), plot = umap_immune,
       dpi = 300, width = 4, height = 4, device = "png")
# Save PDF
cairo_pdf(file = file.path(outdir, "Fig5_umap_immune_revised.pdf"), 
          width=5, height=5, bg = "transparent", family = "Arial")
print(umap_immune)
dev.off()

# Save legend
cairo_pdf(file = file.path(outdir, "Fig5_umap_imm_legend_revised.pdf"), 
          width=5, height=5, bg = "transparent", family = "Arial")
print(imm_legend2)
dev.off()

########################
## Immune FeaturePlot ##
########################

imm_1 <- make_marker_plot(d15_subset, "CD3D")
imm_2 <- make_marker_plot(d15_subset, "CTSS")
imm_3 <- make_marker_plot(d15_subset, "HLA-DPA1")
imm_4 <- make_marker_plot(d15_subset, "NKG7")

grid_imm <- plot_grid(imm_1, imm_2, imm_3, imm_4, nrow = 2)
grid_imm


# Save PNG
ggsave(file.path(outdir, "SuppFig4_featureplot_immune.png"), plot = grid_imm,
       dpi = 300, width = 4, height = 4, device = "png")
# Save PDF
cairo_pdf(file = file.path(outdir, "SuppFig4_featureplot_immune.pdf"), 
          width=4, height=4, bg = "transparent", family = "Arial")
print(grid_imm)
dev.off()

##################################
## Heatmap of bulk RNAseq genes ##
##################################

# Format relative expression data
heatmap_mat <- heat_dat %>%
  column_to_rownames("Gene") %>%
  data.matrix()

# Rows to highlight
enriched_neutr <- c('LIN7A','NCF4','EMR3','FCGR3A','S100P','FCGR3B','TNFRSF10C','C5AR1',
           'MGC31957','CHI3L1','NLRP12','CYP4F3','MXD1','SEPX1','MGAM','DGAT2',
           'RGL4','REPS2','VNN2','CXCR1','NFE2','KRT23','GPR109B','PYGL','FPR2',
           'G0S2','KCNJ15','LRRC4','FPR1','CREB5','FCGR2A','CMTM2','MANSC1',
           'CSF3R','FFAR2','RNF24','LRG1','ST6GALNAC2','ORM1','MME','CDA',
           'PROK2','PFKFB4','VNN3','SLC22A4','BASP1','TREM1','GLT1D1','GPR97',
           'PTAFR','STEAP4','ALPL','NPL','ARAP3','HSPA6','TYROBP','CFD','IMPA2',
           'DENND3','BTNL8','FRAT2','MBOAT7','TSPAN2','BEST1','FLJ10357','SLC40A1')

enriched_mono <-  c('TNFSF13','MYCL1','EPB41L3','DPYSL2','RTN1','SLC31A2','FES','LGALS3',
           'HCK','APLP2','LGALS1','PTGS2','EMR1','AMICA1','DOCK5','CD4','LY96',
           'ARHGEF10L','TNFSF13B','LTBR','PGD','TNFRSF1B','LRRK2','DPYD','MGAM',
           'PTX3','LYZ','IL1R2','DOK3','CARD9','EVI5','GPR109B','DUSP6','MYO1F',
           'FGD4','HHEX','HAL','ST3GAL6','DYSF','RNASE6','SLC24A4','VNN1','NAIP',
           'RHOU','CD68','CXCR2','NACC2','SMARCD3','PADI4','TMEM176B','LGALS3',
           'SAMHD1','CTSS','EMILIN2','ACPP','F5','STEAP4','C19orf59','ACSL1',
           'PAK1','C1orf162','MOSC1','TLR1','PID1','BCL6','HLA-DMB','MPP1','AGPAT9')

myRows <- unique(c(enriched_neutr, enriched_mono))

# Which enriched pathway genes show up in heatmap?
myRows[myRows %in% rownames(heatmap_mat)]

# Arbitrarily select 20 recognizable immune-related genes from myRows
myRows2 <- c("S100P", "FCGR3B", "TNFRSF10C", "C5AR1","CYP4F3", "FCGR2A", "LRG1", "TREM1",
             "HSPA6", "DENND3", "BTNL8", "DOCK5", "LY96", "TNFSF13B", "LRRK2",
             "IL1R2", "CXCR2", "TLR1", "BCL6", "MPP1")

# Create text annotation object for displaying select row names
row_idx <- which(rownames(heatmap_mat) %in% myRows2)

anno = anno_mark(at = row_idx, 
                 labels = rownames(heatmap_mat)[row_idx], 
                 which = "row",
                 padding = 1)

# Draw heatmap

bulk_rna_heatmap <- draw(
  Heatmap(heatmap_mat, 
        column_order = c("1", "3", "15", "56"),
        column_names_side = "bottom",
        column_title = "Study Day",
        column_title_side = "bottom",
        column_names_rot = 0,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        heatmap_legend_param = list(title = "LogFC", direction = "horizontal")
        ) +
  rowAnnotation(mark = anno),
  heatmap_legend_side = "top")
bulk_rna_heatmap

# Save
cairo_pdf(file = file.path(outdir, "Fig3_bulk_rna_heatmap.pdf"), 
          width = 4, height = 6, bg = "transparent", family = "Arial")
print(bulk_rna_heatmap)
dev.off()


#################
## SessionInfo ##
#################

sessionInfo()
