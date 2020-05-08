# !diagnostics off

library(tidyverse)
library(randomizr)
library(extraDistr)
library(parallel)
library(tictoc)
library(rjags)
library(coda)
library(digest)
library(rstanmodels)

source("data.R")
source("util.R")
source("models.R")
source("randomise.R")
source("trials.R")
source("decision.R")


idxscen <- read.csv("scenarios.cfg", row.names = NULL, header = F)[, 1]

#' Simulate over a range of scenarios
#'
#' @return
#' @export
#'
#' @examples
#' i=2; x=1
#' cfg_interface = cfg_fixed4(trial_interface = trial_Fixed_BR, outdir = "outfix4", nsim = 3)
#' cfg_interface = cfg_fixed(trial_interface = trial_Fixed_BR, outdir = "outfix", nsim = 25)
#' cfg_interface = cfg_br_3(trial_interface = trial_GS_BR, outdir = "outbr", nsim = 1)
#' cfg_interface = cfg_rar_3(trial_interface = trial_GS_RAR, outdir = "outrar", nsim = 10)
#' cfg_interface = cfg_rar_4(trial_interface = trial_GS_RAR, outdir = "outrar4", nsim = 10)
sim <- function(cfg_interface = NULL){

  # convenience var
  lpar <- cfg_interface

  dir.create(file.path(lpar$outdir), showWarnings = FALSE)
  

  if(.Platform$OS.type == "unix") {
    ncores <- parallel::detectCores()-1
  } else {
    ncores <- 1
  }
  
  jags_setup()
  

  message(Sys.time(), " Scenarios for current machine ", 
          paste0(idxscen, collapse = " "))

  
  # Only run the scenarios that are defined in config.
  # for(i in 1:nrow(dpars)){
  i=1
  for(i in idxscen){
    
    if(!(i %in% 1:length(lpar$scenarios))){
      message(Sys.time(), " idxscen ", i, " not in scenario range ", 
              paste0(1:length(lpar$scenarios), collapse = " "))
      next
    }
    
    message(Sys.time(), " ", lpar$outdir, ", starting scenario ", lpar$scenarios[[i]]$scen)

    lpar$scen_idx <- i
    lpar$Nmax <- lpar$scenarios[[i]]$Nc
    lpar$mu_n_household <- lpar$scenarios[[i]]$mu_n_household
    lpar$sig_u0 <- lpar$scenarios[[i]]$sig_u0
    lpar$prob_symp <- lpar$scenarios[[i]]$prob_symp
    lpar$interim_at <- lpar$scenarios[[i]]$interim_at
    
    # analyses occur based on number of clusters enrolled
    if(lpar$interim_at < lpar$Nmax ){
      nclustanalys <- seq(lpar$interim_start, lpar$Nmax, by = lpar$interim_at)
      if(max(nclustanalys) == lpar$Nmax){
        nclustanalys <- nclustanalys[-length(nclustanalys)]
      }
      lpar$nclustanalys <- nclustanalys
    } else {
      lpar$nclustanalys <- numeric(0)
    }
    
    if(!is.null(lpar$scenarios[[i]]$arms_activ_at_start)){
      lpar$arms_activ_at_start <- lpar$scenarios[[i]]$arms_activ_at_start
      lpar$arms_enabled_for_anly <- lpar$scenarios[[i]]$arms_enabled_for_anly
      lpar$arms_start_idx <- lpar$scenarios[[i]]$arms_start_idx
    }
    
    if(!is.null(lpar$scenarios[[i]]$thresh_sup)){
      message("Updating superiority threshold to ", lpar$scenarios[[i]]$thresh_sup)
      lpar$thresh_sup <- lpar$scenarios[[i]]$thresh_sup
    }
    
 
    # Some no-brainer checks.
    stopifnot(length(lpar$prob_symp) == length(lpar$arms_activ_at_start))
    stopifnot(lpar$arms_start_idx >= 8)
    

    ld <- mclapply(X=1:lpar$nsim, mc.cores=ncores, FUN=function(x) {
      
      # Response variable for all arms are generated regardless of whether
      # arm is active or not.
      d <- getdat1(lpar, x) 
      
      n_arms <- length(lpar$prob_symp)
      # Always ensure that p_rand, p_best, p_beat_soc, active have length equal to n_arms
      arm_status <- list(K = n_arms,
                         active = lpar$arms_activ_at_start,
                         enabled_for_anly = lpar$arms_enabled_for_anly,
                         is_best = rep(NA, n_arms),
                         is_sup = rep(NA, n_arms),
                         is_inf = rep(NA, n_arms),
                         is_equ = rep(NA, n_arms),
                         sup_at = rep(NA, n_arms),
                         inf_at = rep(NA, n_arms),
                         equ_at = rep(NA, n_arms),
                         arms_in_post = rep(NA, n_arms),
                         p_rand = rep(0, n_arms),
                         p_best = rep(0, n_arms),
                         p_beat_soc = rep(0, n_arms),
                         p_equiv_soc = rep(0, n_arms),
                         var_k = rep(0, n_arms),
                         nk = rep(0, n_arms),
                         nki = rep(0, n_arms),
                         p_emp = rep(0, n_arms))
      
      list(d=d, arm_status = arm_status)
    })
    
    ok <- lapply(1:length(ld), function(x){
      
      # Works out the number of clusters available at the first analysis
      if(length(lpar$nclustanalys)>0){
        recs <- sum(ld[[x]]$d$clustid <= lpar$nclustanalys[1])
      } else {
        recs <- length(ld[[x]]$d$clustid)
      }
      
      if(ld[[x]]$d$entry_time[recs] < lpar$test_at_day){
        message(lpar$outdir, " ERROR pre-data check failure", Sys.time()); 
        message("Entry time of last obs at first interim is lower than test day.")
        message("This implies that the first analysis would have no data.")
        message("Saving relevant dataset, number ", x)
        tim <- format(Sys.time(), "%Y%m%d_%H_%M_%S")
        fil <- paste0("Dataset_at_crash_", tim, ".RDS")
        message("Saving ", fil)
        saveRDS(ld[[x]], fil)
        stop()
      }
    })

    # ldat = ld[[x]]; 

    message("ncores ", ncores)
    message("length(ld) ", length(ld))
    
    
    lres = mclapply(X=1:length(ld), mc.cores=ncores, FUN=function(x) {
      
      lmet <- tryCatch( lpar$trial_interface(scen = lpar$scenarios[[i]]$scen,
                                             d = ld[[x]]$d, 
                                             arm_status = ld[[x]]$arm_status, 
                                             lpar, 
                                             x) ,
                        error=function(e) { 
                          message(get_hash(), " ", lpar$outdir, " HELLO ERROR " , e); 
                          tim <- format(Sys.time(), "%Y%m%d_%H_%M_%S")
                          fil <- paste0("Dataset_at_crash_", tim, ".RDS")
                          message("Saving ", fil)
                          saveRDS(ld[[x]], fil)
                          stop() })
      
      lmet
    })
    
    message(get_hash(), " Simulation ", lpar$outdir, " complete, saving to file.")
    
    # lres contains dataset
    lout <- list(lpar = lpar, lres = lres)

    # Need to create output directory if it does not exist.
    if(i < 10){
      saveRDS(lout, paste0(lpar$outdir, "/sim_dparid_0", i, ".RDS"))
    } else {
      saveRDS(lout, paste0(lpar$outdir, "/sim_dparid_", i, ".RDS"))
    }
  }
}

# Unfortunately, you have to be really careful here. idxscen (see main loop)
# is an index supplied by a config file that tells each VM which set of 
# scenarios to run.

# sim(cfg_fixed4(trial_interface = trial_Fixed_BR, outdir = "outfix4", nsim = 2000))

sim(cfg_rar_3(trial_interface = trial_GS_RAR, outdir = "outrar", nsim = 5000))
#sim(cfg_fixed(trial_interface = trial_Fixed_BR, outdir = "outfix", nsim = 1000))
sim(cfg_rar_4(trial_interface = trial_GS_RAR, outdir = "outrar4", nsim = 5000))






























