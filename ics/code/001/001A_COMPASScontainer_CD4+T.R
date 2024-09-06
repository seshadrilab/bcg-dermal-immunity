###--- R script
# CD4+ T
library(flowWorkspace)
library(tidyverse)
library(doMC)
library(here)
library(data.table)
library(COMPASS)

# ml fhR/4.1.0-foss-2020b


############################################
##############--- BCG TICE ---##############
############################################

####################################################
####################################################
###--- Open flowWorkspace objects
dt.FCS <- read_csv(file = here("data", "BCG-TICE_dt_FCS-info.csv"))
folders <- list.files(path = "/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/data-raw/tmpdata", full.names = TRUE, pattern = "22")
gs_list <- foreach(i = folders) %do%
{
  cat(i, "\n")
  gs <- load_gs(path = i)
  # Processing (TB WCl -> Mtb213)
  if(str_detect(string = i, pattern = "2240")){
    # .FCS files
    fcs <- c("Specimen_025_A7_A07_116.fcs", "Specimen_026_B7_B07_118.fcs", "Specimen_027_C7_C07_120.fcs", "Specimen_028_D7_D07_122.fcs",
             "Specimen_029_E7_E07_124.fcs", "Specimen_030_F7_F07_126.fcs", "Specimen_031_G7_G07_128.fcs", "Specimen_032_H7_H07_130.fcs")
    fcs <- paste0("/fh/fast/_VIDD/HVTN/CHIL_AWS/ICS/RawFCSData/BCG TICE/2240-H-BCG TICE/", fcs)
    dt.FCS %>%
      filter(FILENAME %in% fcs) %>%
      pull(FIL) -> fcs.files
    # pData
    pd <- pData(gs)
    pd <- pd %>%
      mutate(STIM = case_when(name %in% fcs.files ~ "Mtb 213", .default = STIM))
    pData(gs) <- pd
    cols <- c("STIM", "Batch", "PROTOCOL", "CTRSAMPNAME", "STDY_DESC", "SAMP_ORD", "TESTDT", "ASSAY", "SPECROLE", "Run Num", "ASSAYID", "name", "COUNT_METHOD", "PTIDTYPE", "Replicate", "DRAWDT", "PTID", "NETWORK", "VISITNO", "sample_name", "ASSAY_METHOD", "GUAVA_MUSE_ID", "LABID", "Comments", "NUMVIALS", "EXPERIMENT NAME", "roworder", "GUSPEC")
    pData(gs) <- pData(gs)[, cols]
  }
  # Output
  gs %>% return()
}
dt <- enframe(gs_list, value = "flowWorkspace")
dt$name <- folders
dt <- dt %>%
  mutate(phenoData = map(.x = flowWorkspace, .f = function(x) return(x %>% pData)))


####################################################
####################################################
###--- Filtering - Remove duplicate samples
#- Duplicates
# Find duplicate samples
batchList <- lapply(1:nrow(dt), function(i)
{
  unique(pData(dt$flowWorkspace[[i]])[, c("Batch", "PTID", "VISITNO", "STIM", "Run Num", "Replicate")])
})
batchDF <- do.call(rbind, batchList)
batchDF$`Run Num` <- as.numeric(as.character(batchDF$`Run Num`))
# Create ID tag for PTID:VISITNO:STIM
batchDF <- batchDF %>%
  mutate(PTID.VISITNO.STIM = paste(batchDF$PTID, batchDF$VISITNO, batchDF$STIM, sep = ":")) %>%
  arrange(PTID.VISITNO.STIM)
batchDF <- batchDF %>%
  group_by(PTID.VISITNO.STIM) %>%
  mutate(duplicated = ifelse(`Run Num` == max(`Run Num`), FALSE, TRUE))
# Remove duplicates
batchDF <- subset(batchDF, duplicated == 1 & !PTID %in% "CTAC0022")
duplicated.samples <- apply(batchDF[, c("Batch", "PTID", "VISITNO", "STIM", "Run Num", "Replicate")], 1, function(x) paste(x, collapse = ":")) # 0 sample

#- Unreliable PTIDs
unreliable.s <- readxl::read_xlsx(path = here("misc-docs", "BCG Tice Unreliable data_Feb24.xlsx"))
unreliable.s <- unreliable.s[10:13, ]
unreliable.s %>%
  filter(`Unreliable Stim` == "All conditions for all stims") %>%
  mutate(SAMPLE = paste(`Batch #`, PTID, `Visit no.`, RUNNUM, sep = "_")) %>%
  pull(SAMPLE) -> unreliable.s_1
unreliable.s %>%
  filter(!(`Unreliable Stim` == "All conditions for all stims")) %>%
  mutate(SAMPLE = paste(`Batch #`, PTID, `Visit no.`, `Unreliable Stim`, RUNNUM, sep = "_")) %>%
  pull(SAMPLE) -> unreliable.s_2

#- Filtering
# The default filtering is to filter out unreliable as well as out of immunogenicity flag samples before processing.
bind_rows(dt$phenoData) %>%
  filter(!(PTID %in% "CTAC0022")) %>% # Remove control samples
  filter(!(paste(Batch, PTID, VISITNO, STIM, `Run Num`, Replicate, sep = ":") %in% duplicated.samples)) %>% # Remove duplicate samples (RUNNUM)
  filter(!(paste(Batch, PTID, VISITNO, `Run Num`, sep = "_") %in% unreliable.s_1)) %>% # Remove unreliable samples (1)
  filter(!(paste(Batch, PTID, VISITNO, STIM, `Run Num`, sep = "_") %in% unreliable.s_2)) %>% # Remove unreliable samples (2)
  pull(sample_name) -> selected.samples
dt <- dt %>%
  mutate(filtering = map_chr(.x = phenoData, .f = function(x) ifelse(length(intersect(x$sample_name, selected.samples)) != 0, "Y", "N"))) %>%
  filter(filtering == "Y") %>%
  mutate(GatingSet.filtered = map2(.x = flowWorkspace, .y = phenoData, .f = function(x, y) x[which(y$sample_name %in% selected.samples), ] %>% return())) %>%
  mutate(phenoData.filtered = map(.x = GatingSet.filtered, .f = function(x) return(x %>% pData)))


####################################################
####################################################
###--- IFN-g and/or IL-2
#- Extraction
#output_nodes <- c("/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+",
#                  "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/IFNg_OR_IL2")
#results.list <- foreach(i = dt$name) %do%
#{
#  #- Processing
#  cat("Processing", i, "\n")
#  gs <- dt %>%
#    filter(name == i) %>%
#    pull(GatingSet.filtered)
#  gs <- gs[[1]]
#
#  #- Processing & Output
#  exprs.list <- lapply(gs, function(x)
#  {
#    cat("\t", pData(x) %>% rownames(), "\n")
#
#    #- Get boolean positivity call for markers
#    options(warn = 0)
#    marker_response <- try(lapply(output_nodes, function(mrkr){gh_pop_get_indices(x, mrkr)}))
#    while(class(marker_response) == "try-error"){
#      marker_response <- try(lapply(output_nodes, function(mrkr){gh_pop_get_indices(x, mrkr)}))
#    }
#    names(marker_response) <- output_nodes
#    marker_response <- bind_rows(marker_response)
#    colnames(marker_response) <- c("boolean_CD4+", "boolean_IFNg_and_or_IL2+")
#
#    #- Expression, boolean positivity call & add meta.data
#    dt.exprs <- marker_response %>%
#      data.table() %>%
#      mutate(FCS = rownames(pData(x))) %>%
#      mutate(BATCH = pData(x)$Batch) %>%
#      mutate(PTID = pData(x)$PTID) %>%
#      mutate(STIM = pData(x)$STIM) %>%
#      mutate(VISITNO = pData(x)$VISITNO) %>%
#      mutate(RUNNUM = pData(x)$`Run Num`) %>%
#      mutate(REPLICATE = pData(x)$Replicate)
#    dt.exprs$nb_cytokines <- apply(dt.exprs[, c("boolean_IFNg_and_or_IL2+")], 1, function(x) sum(x == TRUE))
#    dt.exprs <- dt.exprs %>%
#      filter(`boolean_CD4+` == TRUE) %>%
#      mutate(nsub = n()) %>%
#      mutate(boolean_cytokines = factor(x = ifelse(nb_cytokines >= 1, TRUE, FALSE), levels = c(TRUE, FALSE))) %>%
#      group_by(BATCH, PTID, STIM, VISITNO, RUNNUM, REPLICATE, nsub, boolean_cytokines, .drop = FALSE) %>%
#      summarize(cytnum = n())
#    dt.exprs %>% return()
#  })
#  bind_rows(exprs.list) %>% return()
#}
#dt.results <- do.call("rbind", results.list)
#saveRDS(object = dt.results, file = here("data", "001_dt_CD4+T_IFNg-and-or-IL2.rds"))

#- Fig.
#dt.results$cytnum %>% summary() # min: 8
#dt.results <- dt.results %>%
#  group_by(BATCH, PTID, VISITNO, STIM, RUNNUM, boolean_cytokines) %>%
#  summarise(nsub = sum(nsub), cytnum = sum(cytnum)) %>%
#  mutate(SAMPLE = paste(PTID, VISITNO, RUNNUM)) %>%
#  ungroup()
# Background subtraction
#dt.tmp_1 <- dt.results %>%
#  filter(boolean_cytokines == TRUE) %>%
#  filter(STIM == "negctrl") %>%
#  rename(nsub_neg = "nsub", cytnum_neg = "cytnum") %>%
#  select(SAMPLE, nsub_neg, cytnum_neg)
#dt.tmp_2 <- dt.results %>%
#  filter(boolean_cytokines == TRUE) %>%
#  filter(STIM != "negctrl")
#dt.ggplot <- merge(x = dt.tmp_2, y = dt.tmp_1, by = "SAMPLE", all.x = TRUE) %>%
#  mutate(pctpos = (cytnum / nsub) * 100) %>%
#  mutate(pctneg = (cytnum_neg / nsub_neg) * 100) %>%
#  mutate(pctpos_adj = pctpos - pctneg) %>%
#  select(BATCH, PTID, STIM, VISITNO, nsub, cytnum, pctpos, nsub_neg, cytnum_neg, pctneg, pctpos_adj) %>%
#  arrange(BATCH, PTID, STIM, VISITNO) %>%
#  mutate(VISITNO = factor(x = VISITNO, levels = c("2", "9", "11", "12", "13", "17")))
# Figure
#dt.ggplot %>%
#  filter(STIM != "sebctrl") %>%
#  mutate(pctpos_adj = case_when(pctpos_adj < 0.001 ~ 0.001, .default = pctpos_adj)) %>%
#  arrange(PTID, VISITNO) %>%
#  ggplot() +
#    geom_boxplot(aes(x = VISITNO, y = pctpos_adj, fill = STIM),
#                alpha = .25, width = .75, outlier.shape = NA) +
#    geom_line(aes(x = VISITNO, y = pctpos_adj, group = PTID),
#              position = position_jitter(width = .1, seed = 123),
#              color = "grey80", linewidth = .25, alpha = .5) +
#    geom_point(aes(x = VISITNO, y = pctpos_adj, fill = STIM),
#               position = position_jitter(width = .1, seed = 123),
#               size = 2, alpha = 1, pch = 21, color = "white") +
#    facet_wrap(~ STIM, ncol = 1) +
#    scale_y_continuous(trans = "log10",
#                       breaks = c(0.001, 0.01, 0.1, 1, 10),
#                       labels = c("<0.001", "0.01", "0.1", "1", "10"),
#                       limits = c(0.0009, 1)) +
#    scale_fill_manual(values = c("#1F78B4", "#33A02C", "#E31A1C")) +
#    labs(x = NULL,
#         y = "% of CD4+ T cells expressing IL-2 and/or IFN-\u03B3\nBackground-subtracted",
#         fill = NULL) +
#    theme_bw() +
#    theme(axis.text = element_text(size = 15),
#          axis.title = element_text(size = 16),
#          title = element_text(size = 16),
#          strip.text.x = element_text(size = 16),
#          strip.text.y = element_text(size = 9),
#          panel.grid = element_blank(),
#          legend.text = element_text(size = 15),
#          legend.position = "none") -> plot
#cairo_pdf(file = here("output", "001_output", "CD4+T_boxplots_IFN-g_and-or_IL-2.pdf"), width = 4, height = 8)
#plot
#dev.off()


####################################################
####################################################
##########--- COMPASS
###--- Create COMPASS container
#- Function
# gs_pop_get_count_fast(x)$Population %>% unique()

# "Ki67 BB660" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/Ki67+"
# "Perforin Alx700" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/Perf+"
# "NKG2C BV650" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/NKG2C+"
# "CCR6 BV711" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/CCR6+"
# "CCR7 BV785" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/R7+"
# "CD45RA BUV496" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/4+RA+"
# "GranLys PE" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/Granulysin+"
# "GM CSF PE-Dazzle594" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/GM-CSF"
# "CXCR3 PE-Cy5" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/CXCR3+"
# "HLA DR PE-Cy5 5" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/DR+"
# "CD153 APC" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/153+"
# "CD154 BUV737" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/154+"
# "IFNy V450" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/IFNg+"
# "IL17a PE-Cy7" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/IL17A_OR_IL17F"
# "IL2 BB700" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/IL2+"
# "IL4/IL13 BB630" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/IL4_OR_IL13"
# "TNFa BUV395" : "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+/TNFa+"

create_compass_container <- function(gs,
                                     node,
                                     markers = c("CD153 APC",
                                                 "CD154 BUV737",
                                                 "IFNy V450",
                                                 "IL17a PE-Cy7",
                                                 "IL2 BB700",
                                                 "IL4/IL13 BB630",
                                                 "GM CSF PE-Dazzle594",
                                                 "TNFa BUV395"),
                                     mp = list("153+" = "CD153 APC",
                                               "154+" = "CD154 BUV737",
                                               "IFNg+" = "IFNy V450",
                                               "IL17A_OR_IL17F" = "IL17a PE-Cy7",
                                               "IL2+" = "IL2 BB700",
                                               "IL4_OR_IL13" = "IL4/IL13 BB630",
                                               "GM-CSF" = "GM CSF PE-Dazzle594",
                                               "TNFa+" = "TNFa BUV395"))
{
  names(mp) <- paste0(node, "/", names(mp))

  # meta.data
  pData(gs)$ptid.visitno <- paste(pData(gs)$PTID, pData(gs)$VISITNO, sep = "_")

  # COMPASS Container
  cc <- COMPASSContainerFromGatingSet(gs,
                                      node = node,
                                      markers = markers,
                                      mp = mp,
                                      individual_id = "ptid.visitno",
                                      countFilterThreshold = 10000)
  cc
}

#- GatingSet list
gs_list <- dt %>%
  pull(GatingSet.filtered)

#- COMPASS container for all GatingSet objects
cc4 <- lapply(gs_list,
              create_compass_container,
              node = "/Time/S/K1/K2/K3/K4/Lv/14-SSlo/S2/L/NotDRhi/3+/3+excl16br/MR1-Va24-/gd-/4+")
# Warning message:
#   In if (!is.na(markers)) { :
#       the condition has length > 1 and only the first element will be used
cc4 <- Reduce(merge, cc4)
cc4

#- Output
saveRDS(cc4, here("data", "001_COMPASS-container_CD4+T.rds"))

