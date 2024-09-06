###--- R script
# CD8+ T
library(flowWorkspace)
library(tidyverse)
library(doMC)
library(here)
library(data.table)
library(COMPASS)
library(glue)
library(readr)

# ml fhR/4.1.0-foss-2020b


####################################################
####################################################
##########--- COMPASS
#- COMPASS
#- Data & Params
cc8 <- readRDS(file = here("data", "002_COMPASS-container_CD8+T.rds"))

#- Foreach stim
ctrl_stim <- "negctrl"
stims <- unique(cc8$meta$STIM)
stims <- stims[stims != ctrl_stim & !is.na(stims)]
stims <- setdiff(stims, "sebctrl")
for(trt_stim in stims){
  # Parameter
  params <- read_rds(glue("/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/output/002_output/tuning/tuned_params_CD8+T_{trt_stim}.rds"))

  # For 5-10 runs
  for(i in 1:10){
    # Variable
    n_iter <- 60000
    n_rep <- 10
    s <- sample(1:100000, 1)

    # Messages
    message(glue("CD8 {trt_stim} COMPASS run"))
    message(glue("control: {ctrl_stim}"))
    message(glue("{n_iter} iterations, {n_rep} replications, seed {s}"))

    # Output files
    outfolder <- glue("/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/BCG_TICE/ICS/output/002_output/tuning/{trt_stim}")
    outfile <- glue("CD8+T_{ctrl_stim}_{trt_stim}_tuned_run_{n_iter}iter_{n_rep}rep_seed{s}.rds")
    outfilepattern <- stringr::str_replace_all(outfile, "\\+", "\\\\+")
    message(glue("Checking for file {outfile} in {outfolder}"))

    # COMPASS
    set.seed(s)
    run_results <- COMPASS(cc8,
                           treatment = STIM == trt_stim,
                           control = STIM == ctrl_stim,
                           iterations = n_iter,
                           replications = n_rep,
                           varp_s = params$svar,
                           varp_u = params$uvar,
                           pp = params$gamma_mix,
                           pb1 = params$gamma_p1,
                           pb2 = params$gamma_p2,
                           category_filter = function(x) colSums(x > 9) > 4,
                           verbose = TRUE)

    # Output
    write_rds(run_results, glue("{outfolder}/{Sys.Date()}_{outfile}"))
  }
}

