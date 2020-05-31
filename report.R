library(tidyverse)

source("decision.R")
source("util.R")



#' Loads and evaluates the results using the decision thresholds provided in
#' the simulation configuration description \code{lpar}.
#' 
#' Todo provide override for decision thresholds.
#'
#' @return
#' @export
#'
#' @examples
get_results <- function(dirc = NULL, debug = F){
  
  
  if(is.null(dirc)){
    message(Sys.time(), " need to specify source dir")
    stop()
  } else {
    message(Sys.time(), " reading files from ", dirc)
  }
  
  if(.Platform$OS.type == "unix") {
    
    dirc <- paste0("//data//", dirc)
    message(Sys.time(), " updated dirc to ", dirc)
  } 
  debug = F
  
  thefiles <- list.files(dirc, "*.RDS")
  
  
  # Each i corresponds to a simulation study.
  dtrial <- tibble()
  darm <- tibble()
  for(i in 1:length(thefiles)){
    
    f <- readRDS(paste0(dirc,"//", thefiles[i]))  
    
    fpre <- substr(thefiles[i], 6, nchar(thefiles[i])-4)
    
    checkdat <- unlist(lapply(1:length(f$lres), 
                              function(i){ 
                                is.null(f$lres[[i]]$d)
                              }))
    
    if(sum(checkdat) > 0){
      message("File: ", thefiles[i])
      message("One or more null datasets. Located at indexes:")
      message(paste0(which(checkdat == 1), collapse=", "))
      next
    }
    
    # Each f$lres contains the results from a series of simulated trails from 
    # a single scenario.
    
    # Just take f$lres[[1]] to get the scenario name.
    
    dtmp1 <- data.frame(scen = f$lres[[1]]$trial_status$scen, 
                        fname = as.character(fpre),
                        nmax = f$lpar$Nmax,
                        nint = length(f$lpar$nclustanalys),
                        intat = f$lpar$interim_at,
                        sig = f$lpar$sig_u0,
                        thresh_sup = f$lpar$thresh_sup,
                        thresh_fut = f$lpar$thresh_fut,
                        thresh_equ = f$lpar$thresh_equ,
                        equ_delta = f$lpar$eq_delta)
    dtmp2 <- extract_trials(f$lres, f$lpar)
    dtmp3 <- extract_arms(f$lres, f$lpar)
    dtmp4<- extract_posts(f$lres, f$lpar)
    
    
    dtmp5 <- merge(dtmp1, t(data.frame(colMeans(dtmp2[,-1])))) %>%
      dplyr::mutate(muc = nki/nk)
    dtrial <- rbind(dtrial, dtmp5)

    dtmp6 <- merge(dtmp4, dtmp3, by = c("sim", "arm")) %>%
      dplyr::arrange(sim, arm) %>%
      dplyr::mutate(scen = f$lres[[1]]$trial_status$scen,
                    fname = as.character(fpre)) 
    darm <- rbind(darm, dtmp6)
    
  }
  
  dtrial <- as_tibble(dtrial)
  darm <- as_tibble(darm)
  
  return(list(dtrial = dtrial, darm = darm))   

}



#' Retrieves data from the actual trials to get a sense of accrual and 
#' stopping under the different scenarios.
#'
#' @param dirc 
#' @param debug 
#'
#' @return
#' @export
#'
#' @examples
get_data <- function(dirc = NULL, debug = F){
  
  if(is.null(dirc)){
    message(Sys.time(), " need to specify source dir")
    stop()
  } else {
    message(Sys.time(), " reading files from ", dirc)
  }
  
  if(.Platform$OS.type == "unix") {
    
    dirc <- paste0("//data//", dirc)
    message(Sys.time(), " updated dirc to ", dirc)
  } 
  debug = F
  
  thefiles <- list.files(dirc, "*.RDS")
  
  dacc <- tibble()
  # Now do the whole thing again to look at accrual
  for(i in 1:length(thefiles)){
    
    f <- readRDS(paste0(dirc,"//", thefiles[i]))  
    
    if(length(f) < 1000){
      message("Note, less than 1000 trials in current file ", thefiles[i])
    }
    
    fpre <- substr(thefiles[i], 6, nchar(thefiles[i])-4)
    
    # Each f$lres contains the results from a series of simulated trails
    dacc1 <- do.call(rbind, lapply(1:length(f$lres), function(j){

        if(f$lres[[j]]$trial_status$sup == 1){
          stopreas <- "Found eff trt(s)"
        } else if(f$lres[[j]]$trial_status$fut == 1){
          stopreas <- "All trt(s) fut"
        } else if(f$lres[[j]]$trial_status$fa == 1){
          stopreas <- "Did not stop until fa"
        } else {
          stop("Unknown")
        }
      
        if(f$lres[[j]]$trial_status$stop == T){
  
          ntmp1 <- f$lres[[j]]$d %>%
            dplyr::mutate(wk = floor(entry_time / 7))  %>%
            dplyr::filter(clustid <= f$lres[[j]]$trial_status$nk) %>%
            dplyr::select(wk, arm, clustid, entry_time, id) %>%
            dplyr::group_by(wk, arm) %>%
            dplyr::summarise(n_k = length(unique(clustid))) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(arm) %>%
            dplyr::arrange(arm, wk) %>%
            dplyr::mutate(n_k = cumsum(n_k)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(scen = f$lpar$scenarios[[i]]$scen,
                          sim = j,
                          stopreas = stopreas,
                          win = f$lres[[j]]$trial_status$win,
                          fa = f$lres[[j]]$trial_status$fa,
                          p = f$lpar$prob_symp[arm]) %>%
            dplyr::select(scen, sim, arm, wk, everything())
   
          
        } else {
          
          ntmp1 <- f$lres[[j]]$d %>%
            dplyr::mutate(wk = floor(entry_time / 7))  %>%
            dplyr::select(wk, arm, clustid, entry_time, id) %>%
            dplyr::group_by(wk, arm) %>%
            dplyr::summarise(n_k = length(unique(clustid))) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(arm) %>%
            dplyr::arrange(arm, wk) %>%
            dplyr::mutate(n_k = cumsum(n_k)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(scen = f$lpar$scenarios[[i]]$scen,
                          sim = j,
                          stopreas = stopreas,
                          win = f$lres[[j]]$trial_status$win,
                          fa = f$lres[[j]]$trial_status$fa,
                          p = f$lpar$prob_symp[arm]) %>%
            dplyr::select(scen, sim, arm, wk, everything())
  
        }
        return(ntmp1)
  
    }))
    
    dacc <- bind_rows(dacc,
                      dacc1 )
    
  }
  dacc
}


#' Extracts the results from a series of simulations performed under a 
#' trial configuration specified by \code{lpar}.
#'
#' @param lres 
#' @param lpar 
#'
#' @return
#' @export
#'
#' @examples
#' lres <- f$lres
#' lpar <- f$lpar
extract_trials <- function(lres = NA, 
                          lpar = NA){
  
  trials <- do.call(rbind, lapply(1:length(lres), function(i){
    
    if(length(unlist(lres[[i]]$trial_status)) != 11){
      stop("len not 11, ", i, " len = ", length(unlist(lres[[i]]$trial_status)))
    }
    
    c(simidx = i, unlist(lres[[i]]$trial_status[-1]), nsim = length(lres))
    
  }))
  trials[is.na(trials)] <- 0
  
  trials
  
}


#' Extracts the arm_status data
#'
#' @param lres 
#' @param lpar 
#'
#' @return
#' @export
#'
#' @examples
#' lres <- f$lres
#' lpar <- f$lpar
extract_arms <- function(lres = NA, 
                           lpar = NA){
  
  mynam <- c("sim", "arm", "p_tru", names(lres[[1]]$arm_status)[-1])
  
  arms <- do.call(rbind, lapply(1:length(lres), function(i){
    
    armt <- do.call(cbind, lapply(lres[[i]]$arm_status, function(x){
      if(length(x)>1){
        matrix(x, ncol = 1)
      }
      }))
    
    armt <- cbind(scen = i, 
                  arm = 1:length(lres[[i]]$arm_status$active), 
                  p_tru = lpar$prob_symp, 
                  armt)

    armt
    }))
  
  colnames(arms) <- mynam
  
  # arms[, "p_best"][arms[, "p_best"]==DUMVAL] = 0
  # arms[, "p_beat_soc"][arms[, "p_beat_soc"]==DUMVAL] = 0
  # arms[, "p_equiv_soc"][arms[, "p_equiv_soc"]==DUMVAL] = 0
  # arms[, "var_k"][arms[, "var_k"]==DUMVAL] = 0
  
  arms[,"is_sup"][is.na(arms[,"is_sup"])] <- 0
  arms[,"is_inf"][is.na(arms[,"is_inf"])] <- 0
  arms[,"is_equ"][is.na(arms[,"is_equ"])] <- 0
  
  arms
  
}

#' Extracts data from posterior associated with each simulation.
#'
#' @param lres 
#' @param lpar 
#'
#' @return
#' @export
#'
#' @examples
#' lres <- f$lres
#' lpar <- f$lpar
extract_posts <- function(lres = NA, 
                         lpar = NA){
  
  mynam <- c("sim", names(lres[[1]]$arm_status)[-1])
  
  arms <- do.call(rbind, lapply(1:length(lres), function(i){
    
    p <- apply(lres[[i]]$post[, -ncol(lres[[i]]$post)], 2, expit)
    p_mu <- matrix(colMeans(p), ncol = 1)
    p_var <- apply(p, 2, var)
    armt <- cbind(scen = i, arm = 1:ncol(p), p_mu = p_mu, p_var = p_var)
    colnames(armt) <- c("sim", "arm", "p_est", "p_var")
    armt
    
  }))
  
  arms
  
}









#' Retrieves all the data filenames.
#'
#' @param dirc 
#'
#' @return
#' @export
#'
#' @examples
get_fnames <- function(dirc = NULL){
  
  if(is.null(dirc)){
    message(Sys.time(), " need to specify source dir")
    stop()
  } else {
    message(Sys.time(), " reading files from ", dirc)
  }
  
  if(.Platform$OS.type == "unix") {
    
    dirc <- paste0("//data//", dirc)
    message(Sys.time(), " updated dirc to ", dirc)
  } 
  debug = F
  
  thefiles <- list.files(dirc, "*.RDS")
  thefiles
}




