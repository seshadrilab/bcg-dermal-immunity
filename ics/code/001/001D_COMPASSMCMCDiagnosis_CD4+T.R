###--- R script
# CD4+ T
library(flowWorkspace)
library(tidyverse)
library(doMC)
library(here)
library(data.table)
library(COMPASS)

# ml fhR/4.1.0-foss-2020b

####################################################
####################################################
###--- COMPASSMCMCDiagnosis
#- BCG TICE
# COMPASSMCMCDiagnosis
model_files <- list.files("/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/output/001_output/tuning/BCG TICE", full.names = TRUE, pattern = "CD4[+]T_negctrl_BCG TICE")
models <- lapply(model_files, readRDS)
diagnosis <- COMPASSMCMCDiagnosis(models)
# Output
saveRDS(object = diagnosis, file = "/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/output/001_output/tuning/BCG TICE/CD4+T_BCG TICE_COMPASS-dx.rds")

#- Mtb 213
# COMPASSMCMCDiagnosis
model_files <- list.files("/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/output/001_output/tuning/Mtb 213", full.names = TRUE, pattern = "CD4[+]T_negctrl_Mtb 213")
models <- lapply(model_files, readRDS)
diagnosis <- COMPASSMCMCDiagnosis(models)
# Output
saveRDS(object = diagnosis, file = "/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/output/001_output/tuning/Mtb 213/CD4+T_Mtb 213_COMPASS-dx.rds")

#- TB WCL
# COMPASSMCMCDiagnosis
model_files <- list.files("/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/output/001_output/tuning/TB WCL", full.names = TRUE, pattern = "CD4[+]T_negctrl_TB WCL")
models <- lapply(model_files, readRDS)
diagnosis <- COMPASSMCMCDiagnosis(models)
# Output
saveRDS(object = diagnosis, file = "/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/output/001_output/tuning/TB WCL/CD4+T_TB WCL_COMPASS-dx.rds")

