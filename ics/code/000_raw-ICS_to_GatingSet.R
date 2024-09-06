library(here)
library(tidyverse)
library(data.table)
library(doMC)

# ml fhR/4.1.0-foss-2020b

# ------------------------------------------------------------------------------
# Valentin Voillet
# 12th FEB 2024
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
### STEP 1: Process Master Thaw List / Batch File ------------------------------
#
# This is done outside of R. Ask the lab for the Master Thaw List for this project
# and copy it this project's directory. The file will typically have three header rows.
# Process this file as follows:
# 1. Delete all header rows EXCEPT the one with Batch #, PTID, VISITNO, etc.
# 2. Rename Batch # to Batch
# 3. Save as a csv to your project's local directory
#
# Once you have a properly-formatted master thaw list in your local directory,
# you can continue to parsing the XML file for keyword information.
batchData <- readxl::read_xlsx(path = here("misc-docs", "BCG Tice MTL_07Aug23_modified-by-VV_12FEB24.xlsx"))

##-- Useful diagnostics
table(batchData$Batch, batchData$ASSAYID, useNA = "always") # Ensure Batch matches ASSAYID
tapply(batchData$PTID, batchData$VISITNO, function(x) length(unique(x))) # Number of PTIDs by VISITNO
reshape2::dcast(PTID ~ VISITNO, data = batchData) # Each PTID has expected number of samples
# VISITNO: 2, 9, 10.2 (CTAC0022), 11, 12, 13 & 17



# ------------------------------------------------------------------------------
### STEP 2: Match XML workspace files with FCS batch directories ---------------
#
# For each batch, a Gating Set is created using flowWorkpace::parseWorkspace()
# These objects are a combination of raw FCS data combined with the assay-specific
# metadata contained in the workspace XML file. When you run parseWorkspace,
# it gates (read as defines subsets) each FCS file and combines them into a single object.

# Make a vector of XML workspace file paths. You might have to do some work
# if the files are spread across several directories.

##-- Copy-paste XML files in data-raw
# .XML files are manullay copied-pasted
# /fh/fast/_VIDD/HVTN/CHIL_AWS/ICS/Clinical Trial Workspaces/BCG TICE


##-- .xml files & batch
xml_files <- list.files(path = here("data-raw", "xml_files"))
xml.path_files <- list.files(path = here("data-raw", "xml_files"), full.names = TRUE)
dt <- data.table(assayid = sapply(X = xml_files, FUN = function(x) str_split(string = x, pattern = "-")[[1]][1]),
                 xml = xml.path_files)


##-- .fcs directories, .xml files & batch
ROOT <- "/fh/fast/_VIDD/HVTN/CHIL_AWS/ICS/RawFCSData/BCG TICE"
folders <- list.files(path = ROOT, pattern = "BCG TICE", full.names = TRUE)
folders <- folders[!str_detect(string = folders, pattern = "QC")]
folders <- folders[-2] # Batch 2241 is excluded
dt <- dt %>%
  mutate(fcs = folders)


##-- Output -  Save .csv - workspace
# Write to .csv. After performing step 3, it might be worthwhile to edit this file by
# hand to add XML keyword and group names so we have it on record.
write_csv(x = dt, file = here("data-raw", "tmpdata", "batch_to_workspace_map.csv"))



# ------------------------------------------------------------------------------
### STEP 3: Parse XML files for keywords ---------------------------------------
#
# The values to arguments of parseWorkspace may change from study to study, so
# this script helps automate parsing of the XML workspace files to identify
# key arguments. The following information needs to extracted from the XML file:
#
# workspace: This is the XML file parsed by openWorkspace()
# name: Sample Group to process. Typically it's "Test samples"
# path: Path to directory of ICS files corresponding to the XML batch
# keywords: Information to extract from the XML files to use as metadata.
#   This is typically c("$FIL", "Stim", "Sample Order", "EXPERIMENT NAME", "Replicate")
#   Yes, case matters!
#
# Create dataset for each ASSAYID, with all information to generate compass gating set ---
# Need ASSAYID, XML workspace file, ICS/FCS dir, sample name, and keywords ---

##-- Functions
get_workspace_keywords <- function(xml, keyword_vars) {

  library(xml2)
  ws <- read_xml(xml)

  # Find Group Names, the name argument to parseWorkspace.
  gn <- xml_find_all(ws, ".//GroupNode")
  all_groups <- sapply(gn, function(x){
    groups <- xml_attrs(x)[c("nodeName", "owningGroup")]
    unique(groups)
  })

  # Check that our keywords exist in the XML file.
  kw <- xml_find_all(ws, ".//KeywordInfo")
  kwc <- xml_children(kw)
  all_keywords <- sapply(xml_attrs(kwc), function(x) x[["name"]])
  keywords_keep <- all_keywords[tolower(all_keywords) %in% tolower(keyword_vars)]

  # Return a character vector of keywords and names to quickly check.
  c(xml, keywords_keep, all_groups)
}


##-- Gating set info
dt.workspace <- read_csv(file = here("data-raw", "tmpdata", "batch_to_workspace_map.csv"))
keyword_vars <- c("$FIL", "Stim", "Sample Order", "EXPERIMENT NAME", "Replicate")
gatingset_info <- foreach(i = dt.workspace$xml) %do% {
  cat(i)
  get_workspace_keywords(i, keyword_vars)
}
gatingset_info

##-- Output
saveRDS(gatingset_info, here("data-raw", "tmpdata", "gating_set_info.rds"))



# ------------------------------------------------------------------------------
### STEP 4: Create gating set objects ------------------------------------------
library(here)
library(tidyverse)
library(data.table)
library(doMC)
library(XML)
library(flowWorkspace)
library(CytoML)
library(parallel)


##-- Batch & Data
batchData <- readxl::read_xlsx(path = here("misc-docs", "BCG Tice MTL_07Aug23_modified-by-VV_12FEB24.xlsx"))
dt.workspace <- read_csv(file = here("data-raw", "tmpdata", "batch_to_workspace_map.csv"))


##-- GatingSet function
create_gating_set <- function(assayid,
                              xml_path,
                              fcs_path,
                              sample_name = "Samples",
                              xml_keywords = c("$FIL", "$TOT", "EXPERIMENT NAME", "Sample Order", "Stim", "Replicate")) {
  cat(assayid, "\n")

  #- GatingSet
  workspace <- CytoML::open_flowjo_xml(xml_path)
  G <- CytoML::flowjo_to_gatingset(workspace,
                                   name = sample_name,
                                   keywords = xml_keywords,
                                   path = fcs_path,
                                   additional.sampleID = TRUE)

  #- Add annotation data to gating set
  batchdf <- batchData[which(batchData$Batch == assayid), ]
  pd <- flowWorkspace::pData(G)
  pd$Batch <- unique(batchdf$Batch)
  # Find Sample Order column and convert to numeric. SAMP_ORD is key for merge with batch data
  pd$SAMP_ORD <- as.numeric(as.character(pd[, which(tolower(names(pd)) %in% "sample order")]))
  pd$STIM <- pd[, which(tolower(names(pd)) %in% "stim")]
  pd$roworder <- 1:nrow(pd)
  pd$sample_name <- rownames(pd)

  #- Add batchData information to pData, then apply to GatingSet.
  pd <- merge(x = pd[,c("name", "SAMP_ORD", "Batch", "STIM", "Replicate", "EXPERIMENT NAME", "roworder", "sample_name")],
              y = batchdf,
              by = c("SAMP_ORD", "Batch"), all.x = TRUE)
  rownames(pd) <- pd[["sample_name"]]
  pd <- pd[order(pd$roworder), ]
  flowWorkspace::pData(G) <- pd

  #- Output
  flowWorkspace::save_gs(G, path = here::here("data-raw", "tmpdata", assayid), overwrite = TRUE)
}


##-- GatingSet (investigation / test)
#- BATCH: 2240
assayid <- dt.workspace$assayid[1]
xml_path <- dt.workspace$xml[1]
fcs_path <- dt.workspace$fcs[1]
sample_name <- "Samples"
xml_keywords <- c("$FIL", "$TOT", "EXPERIMENT NAME", "Sample Order", "Stim", "Replicate")
#- GatingSet
workspace <- CytoML::open_flowjo_xml(xml_path)
G <- CytoML::flowjo_to_gatingset(workspace,
                                 name = sample_name,
                                 keywords = xml_keywords,
                                 path = fcs_path,
                                 includeGates = TRUE,
                                 additional.sampleID = TRUE)
#- Add annotation data to gating set
batchdf <- batchData[which(batchData$Batch == assayid), ]
pd <- flowWorkspace::pData(G)
pd$Batch <- unique(batchdf$Batch)
# Find Sample Order column and convert to numeric. SAMP_ORD is key for merge with batch data
pd$SAMP_ORD <- as.numeric(as.character(pd[, which(tolower(names(pd)) %in% "sample order")]))
pd$STIM <- pd[, which(tolower(names(pd)) %in% "stim")]
pd$roworder <- 1:nrow(pd)
pd$sample_name <- rownames(pd)
#- Add batchData information to pData, then apply to GatingSet.
pd <- merge(x = pd[,c("name", "SAMP_ORD", "Batch", "STIM", "Replicate", "EXPERIMENT NAME", "roworder", "sample_name")],
            y = batchdf,
            by = c("SAMP_ORD", "Batch"), all.x = TRUE)
rownames(pd) <- pd[["sample_name"]]
pd <- pd[order(pd$roworder), ]
flowWorkspace::pData(G) <- pd
#- Output
flowWorkspace::save_gs(G, path = here::here("data-raw", "GatingSet_2240_INVESTIGATION"), overwrite = TRUE)


##-- Running
cl <- makeCluster(3)
clusterExport(cl, c("batchData", "dt.workspace", "create_gating_set"))
parApply(cl, dt.workspace, 1, function(x) {
  create_gating_set(
    assayid = x[["assayid"]],
    xml_path = x[["xml"]],
    fcs_path = x[["fcs"]]
  )
})
stopCluster(cl)
# Error - BATCH 2240
# Error in checkForRemoteErrors(val) :
# 2 nodes produced errors; first error: Error in .cpp_saveGatingSet(gs@pointer, path = path, backend_opt = backend_opt,  :
#                                                                     Not a valid GatingSet archiving folder! /fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/data-raw/tmpdata/2240
#                                                                   gs file not matched to GatingSet uid: /fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/data-raw/tmpdata/2240/43d7746e-ccd3-4eec-a17c-2cc61e6b2e66.gs
# Error - BATCH 2242
# Error in checkForRemoteErrors(val) :
#  one node produced an error: Error in .cpp_saveGatingSet(gs@pointer, path = path, backend_opt = backend_opt,  :
#                                                            Not a valid GatingSet archiving folder! /fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/data-raw/tmpdata/2242
#                                                          gs file not matched to GatingSet uid: /fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/data-raw/tmpdata/2242/74e156b3-0852-46a2-89f1-aa7aaa525ce3.gs



# ------------------------------------------------------------------------------
### STEP 5: .FCS processing ----------------------------------------------------
#- Data
dt.workspace <- read_csv(file = here("data-raw", "tmpdata", "batch_to_workspace_map.csv"))

#- Processing
# For each batch, get .fcs files related to Samples
registerDoMC(3)
dt.FCS <- foreach(i = 1:nrow(dt.workspace)) %dopar% {
  # Folders, Files & ASSAYID
  folder.fcs <- dt.workspace$fcs[i]
  assayid <- dt.workspace$assayid[i]
  files <- list.files(path = folder.fcs, pattern = ".fcs", full.names = TRUE)
  # data.frame() with .FCS & description
  cat(assayid, "\n")
  dt.tmp <- foreach(file = files) %do%
  {
    fcs.tmp <- flowCore::read.FCS(filename = file)
    data.table(ASSAYID = assayid,
               FILENAME = fcs.tmp@description$FILENAME,
               FIL = fcs.tmp@description$`$FIL`) %>%
      return()
  } %>% bind_rows()
  # Output
  dt.tmp %>% return()
} %>% bind_rows()

#- Output
dt.FCS %>%
  write_csv(file = here("data", "BCG-TICE_dt_FCS-info.csv"))

