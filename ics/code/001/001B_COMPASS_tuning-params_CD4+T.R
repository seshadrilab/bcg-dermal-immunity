###--- R script
# CD4+ T
library(flowWorkspace)
library(tidyverse)
library(doMC)
library(here)
library(data.table)
library(COMPASS)
source(file = here("misc-docs", "R_COMPASS", "tuning_fns.R"))

# ml fhR/4.1.0-foss-2020b


####################################################
####################################################
##########--- COMPASS
#- Data
cc4 <- readRDS(file = here("data", "001_COMPASS-container_CD4+T.rds"))

#- Set up control & stimulations
ctrl_stim <- "negctrl"
stims <- unique(cc4$meta$STIM)
stims <- stims[stims != ctrl_stim & !is.na(stims)]
stims <- setdiff(stims, "sebctrl")
# Number of iterations per round
n_iter <- 1000
# Number of replicates per round
n_rep <- 2
# Number of rounds to run
n_rounds <- 20
# Reset limits for tuning so we don't get stuck
forget_every <- 5

#- Tuning
for (stim in stims)
{
  stimdir <- glue("/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/output/001_output/tuning/{stim}")
  dir.create(stimdir, recursive = TRUE, showWarnings = FALSE)
  trt_stim <- stim

  # Tune the parameters - output are saved to <stimdir>
  tuning_results <- tune_parameters(cc4,
                                    trt_stim,
                                    ctrl_stim,
                                    min_count = 9,
                                    min_cells = 4,
                                    remove_markers = NULL,
                                    n_iter,
                                    n_rep,
                                    n_rounds,
                                    forget_every,
                                    save_results = TRUE,
                                    out_folder = stimdir,
                                    file_prefix = glue("tune_cdr_{stim}"),
                                    verbose = TRUE)

  # Save the final set of parameters for our runs
  write_rds(x = tuning_results$params[[length(tuning_results$params)]], file = glue("/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/output/001_output/tuning/tuned_params_CD4+T_{stim}.rds"))
}

