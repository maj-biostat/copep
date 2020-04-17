#' Fixed size K arm trial
#'
#' @param d 
#' @param lpar 
#'
#' @return
#' @export
#'
#' @examples
trial_Fixed_BR <- function(d, arm_status, lpar, x){
  
  set_hash()
  
  message(get_hash(), " trial_Fixed_BR, starting with dataset number ", x)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility/equivalence - all arms futile/equivalent
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       fa = F,  # ran all the way to the final analysis
                       na = 0,  # current analysis number
                       nk = 0,  # current total number of clusters
                       nki = 0, # current total number of participants 
                       durn = 0)# elapsed duration of trial  

  # Final analysis
  message(get_hash(), " trial_Fixed_BR, starting final")
  
  # final analysis
  # randomise per the last interim
  
  d$arm <- randomizr::cluster_ra(d$clustid, 
                                 prob_each = rep(1/arm_status$K, arm_status$K), 
                                 num_arms = arm_status$K)
  d$arm <- as.numeric(substr(as.character(d$arm), 2, 2))
  
  d$y <- unlist(lapply(1:nrow(d), function(j){
    y <- d[j, lpar$arms_start_idx + d$arm[j]]
    y
  }))
  
  trial_status$fa <- T
  trial_status$nki <- nrow(d)
  trial_status$nk <- d$clustid[nrow(d)]
  trial_status$durn <- d$entry_time[nrow(d)]
  
  message(get_hash(), " trial_Fixed_BR, modelling final ")
  
  post <- glmm_jags(d)
  
  
  # tpost <- post[, 1:3]
  # 
  # colMeans(tpost)
  # hist(tpost[, 1])
  # 
  # quantile(tpost[, 1], probs = c(0.025, 0.975))
  # 
  # diff1_2 <- tpost[, 1]-tpost[, 2]
  # quantile(diff1_2, probs = c(0.025, 0.5, 0.975))
  # diff1_3 <- tpost[, 1]-tpost[, 3]
  # hist(diff1_2)
  # hist(abs(diff1_2))
  # 
  # p <- apply(tpost, 2, expit)
  # quantile(p[, 1], probs = c(0.025, 0.5, 0.975))
  # hist(p[, 1])
  # apply(p, 2, sd)
  # 
  # diff1_2 <- p[, 1]-p[, 2]
  # round(quantile(diff1_2, probs = c(0.025, 0.5, 0.975)), 3)
  # diff1_3 <- p[, 1]-p[, 3]
  # 
  # hist(abs(diff1_2))
  
  decis <- decision_fa(post, 
                       a_s = arm_status, 
                       t_s = trial_status, 
                       lpar)
  
  trial_status <- decis$trial_status
  arm_status <- decis$arm_status
  
  lret <- list(post = post,
               d = d,
               trial_status = trial_status,
               arm_status = arm_status)
  
  lret
}





















#' Group sequential trial with K arms and balanced randomisation (but no blocking)
#' and with early stopping over interims configured via \code{lpar}. Invokes 
#' analyses at each interim and final analysis and returns the posterior samples 
#' from each analysis.
#'
#' @param d 
#' @param lpar 
#'
#' @return
#' @export
#'
#' @examples
trial_GS_BR <- function(ldat, lpar, x){
  
  if(x==1 | x %% 50 == 0){
    message(get_hash(), " trial_GS_BR, starting with dataset number ", x)
  }
  
  fa <- list()
  ia <- list()
  
  # Interim analyses
  if(length(lpar$nclustanalys) > 0){
    
    for(i in 1:length(lpar$nclustanalys)){
      
      # allocate to arms, run interim and break if we triggered a stopping rule
      # Note - balanced allocation for first interim
      message(get_hash(), " trial_GS_BR, starting interim ", i)
      
      idxstart <- ifelse(i==1, 1, nrow(ldat$d) - sum(ldat$d$clustid > lpar$nclustanalys[i-1])+1  )
      idxend <- sum(ldat$d$clustid <= lpar$nclustanalys[i])
      
      message(get_hash(), " trial_GS_BR, alloc records ", idxstart, " to ", idxend)
      
      ldat$d$arm[idxstart:idxend] <- BR_alloc(ldat, idxstart, idxend, interim_idx = i)
      
      ldat$d$y[idxstart:idxend] <- unlist(lapply(idxstart:idxend, function(j){
        y <- ldat$d[j, lpar$arms_start_idx + ldat$d$arm[j]]
        y
      }))
      
      tmp <- intrm_GS(dint = ldat$d[1:idxend,], 
                      cur_arm_status = ldat$arm_status, 
                      lpar)
      ia[[i]] <- list(lpost_ia = tmp$lpost_ia,
                      ldecis_ia = tmp$ldecis_ia$vret)
      ldat$arm_status <- tmp$ldecis_ia$arm_status
      
      # if any stop then break
      if(any(ia[[i]]$ldecis_ia$stop == 1)){
        message("Triggered interim stopping rule.")
        ia[[i]]$ldecis_ia$ss <- lpar$nclustanalys[i]
        ia[[i]]$ldecis_ia$durn <- ldat$d$entry_time[lpar$nclustanalys[i]]
        break
      }
      
    }
    
  } 
  
  # Final analysis
  # This only occurs if we did not stop at an interim
  
  # This gives you the assessment of p_best at the final analysis if we 
  # did not stop at an interim.
  # If we did stop at an interim then the arm_status reflects the state of
  # the evidence when we stopped (ie prior to the final analysis).
  if(all(ia[[i]]$ldecis_ia$stop == 0)){
    
    message(get_hash(), " trial_GS_RAR, starting final ")
    
    # final analysis
    # randomise per the last interim

    ldat$d$arm[(idxend+1):nrow(ldat$d)] <- BR_alloc(ldat, idxend+1, nrow(ldat$d), interim = 99)
    
    ldat$d$y[(idxend+1):nrow(ldat$d)] <- unlist(lapply((idxend+1):nrow(ldat$d), function(j){
      y <- ldat$d[j, lpar$arms_start_idx + ldat$d$arm[j]]
      y
    }))
    
    fa$post_fa <- glmm_jags(ldat$d)
    
    tmp <- decision_fa(post = fa$post_fa, 
                       cur_arm_status = ldat$arm_status, lpar)
    
    fa$ldecis_fa <- tmp$vret
    fa$ldecis_fa$durn <- max(ldat$d$entry_time)

    # This gives you the assessment of p_best at the final analysis if we 
    # did not stop at an interim.
    # If we did stop at an interim then the arm_status reflects the state of
    # the evidence when we stopped (ie prior to the final analysis).
    ldat$arm_status <- tmp$arm_status
  }
  
  
  

  # Now I get rid of all the redundant interim analyses and just keep the 
  # last one.
  if(lpar$keep_intrms == F & length(lpar$nclustanalys) > 0){
    
    ia <- list(lpost_ia = ia[[i]]$lpost_ia,
               ldecis_ia = ia[[i]]$ldecis_ia)
    
  }
  
  # It is quite valid for fa to be of zero length - i.e. if we stopped at 
  # an interim analysis we would not undertake a final analysis.
  lret <- list(fa = fa,
               ia = ia,
               ldat = ldat)
  
  lret
}




#' Runs interim analyses
#'
#' @param d 
#' @param lpar 
#'
#' @return
#' @export
#'
#' @examples
intrm_GS <- function(dint, a_s, t_s, lpar){
  
  # Adjust the data available based on time elapsed from entry.
  # Only those that have hit the endpoint are included in the analysis.
  tbl <- table(dint$arm)
  
  message(get_hash(), " interim data on arms ", 
          paste0(names(tbl), collapse = " "), " with ",
          paste0(as.numeric(tbl), collapse = " "), " with ",
          nrow(dint), " recs in dataset and max entry time ",
          max(dint$entry_time, na.rm = T), " days")
  
  dtmp <- dint %>% dplyr::group_by(arm) %>% dplyr::summarise(nk = length(unique(clustid)))
  a_s$nk[dtmp$arm] <- dtmp$nk
  
  dobs <- dint %>% 
    dplyr::mutate(maxentry = max(entry_time)) %>%
    dplyr::filter(entry_time <= maxentry - (lpar$test_at_day)) %>%
    dplyr::select(-maxentry)
  
  tbl <- table(dobs$arm)
  
  message(get_hash(), " observed data on arms ", 
          paste0(names(tbl), collapse = " "), " with ",
          paste0(as.numeric(tbl), collapse = " "), " with ",
          nrow(dobs), " recs in dataset and max entry time ",
          max(dobs$entry_time, na.rm = T), " days (fu is " , lpar$test_at_day, " days).")
  

  post <- glmm_jags(dobs)
  
  decis <- decision_intrm(post, a_s, t_s, lpar)

  return(list(post = post, 
              decis = decis))
  
}



#' Group sequential trial with K arms and response adaptive randomisation 
#' and with early stopping over interims configured via \code{lpar}. Invokes 
#' analyses at each interim and final analysis and returns the posterior samples 
#' from each analysis.
#'
#' @param d 
#' @param lpar 
#'
#' @return
#' @export
#'
#' @examples
trial_GS_RAR <- function(d, arm_status, lpar, x){
  
  set_hash()
  
  message(get_hash(), " trial_GS_RAR, starting trial ", x)

  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility/equivalence - all arms futile/equivalent
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       fa = F,  # ran all the way to the final analysis
                       na = 0,  # current analysis number
                       nk = 0,  # current total number of clusters
                       nki = 0, # current total number of participants 
                       durn = 0)# elapsed duration of trial  
  
  
  # Interim analyses
  if(length(lpar$nclustanalys) > 0){
    
    i=1
    for(i in 1:length(lpar$nclustanalys)){
      
      # allocate to arms, run interim and break if we triggered a stopping rule
      # Note - balanced allocation for first interim
      message(get_hash(), " trial_GS_RAR, starting interim ", i, 
              " p_rand ", 
              paste0(round(arm_status$p_rand, 3), collapse = " "))
      
      trial_status$na <- i
      
      
      # checking whether arms should be activated
      for(m in 2:arm_status$K){
        if(arm_status$enabled_for_anly[m] == i){
          message(get_hash(), " trial_GS_BR, activating arm ", m)
          arm_status$active[m] <- T
        }
      }
      
      
      # Number of participants between end of the previous to this interim.
      idxstart <- ifelse(i==1, 1, nrow(d) - sum(d$clustid > lpar$nclustanalys[i-1])+1  )
      idxend <- sum(d$clustid <= lpar$nclustanalys[i])
      
      message(get_hash(), " trial_GS_RAR, alloc records ", idxstart, " to ", idxend)
      trial_status$nki <- idxend
      trial_status$nk <- d$clustid[idxend]
      trial_status$durn <- d$entry_time[idxend]
      
      # This is the allocation that happened up to this interim. I know that is a bit confusing.
      # Basically, I am say the analysis up to this time point was based on the following 
      # randomisation scheme. At the first analysis, we didn't have any previous data so all
      # we can do is balanced allocation.
      tmp <- RAR_alloc(clustid = d$clustid, 
                       a_s = arm_status,
                       idxstart, 
                       idxend, 
                       interim_idx = i)

      
      d$arm[idxstart:idxend] <- tmp$rand_arm
      arm_status <- tmp$arm_status
      d$y[idxstart:idxend] <- unlist(lapply(idxstart:idxend, function(j){
        y <- d[j, lpar$arms_start_idx + d$arm[j]]
        y
      }))
      
      
      tmp <- intrm_GS(dint = d[1:idxend,], 
                      a_s = arm_status, 
                      t_s = trial_status,
                      lpar)
      
      post <- tmp$post
      arm_status <- tmp$decis$arm_status
      trial_status <- tmp$decis$trial_status
      
      
      # if any stop then break
      if(trial_status$stop){
        
        message(get_hash(), " Triggered interim stopping rule!")
        break
        
      }
      
    }
    
  } 
  
  
  
  # Final analysis
  # This only occurs if we did not stop at an interim
  if(!trial_status$stop){

    message(get_hash(), " trial_GS_RAR, starting final interims ended at ", 
            lpar$nclustanalys[i], " max ss ", lpar$Nmax)
 
    tmp <- RAR_alloc(clustid = d$clustid, 
                     a_s = arm_status,
                     idxend+1, 
                     nrow(d), 
                     interim = 99) # bogus interim id just needs to be bigger than max interim
    
    d$arm[(idxend+1):nrow(d)] <- tmp$rand_arm
    arm_status <- tmp$arm_status
  
    d$y[(idxend+1):nrow(d)] <- unlist(lapply((idxend+1):nrow(d), function(j){
      y <- d[j, lpar$arms_start_idx + d$arm[j]]
      y
    }))
    
    trial_status$fa <- T
    trial_status$nki <- nrow(d)
    trial_status$nk <- d$clustid[nrow(d)]
    trial_status$durn <- d$entry_time[nrow(d)]
    
    message(get_hash(), " trial_GS_RAR, modelling final ")
    
    post <- glmm_jags(d)
    
    decis <- decision_fa(post, 
                         a_s = arm_status, 
                         t_s = trial_status, 
                         lpar)

    trial_status <- decis$trial_status
    arm_status <- decis$arm_status

  }
  
 
  # It is quite valid for fa to be of zero length - i.e. if we stopped at 
  # an interim analysis we would not undertake a final analysis.
  lret <- list(post = post,
               d = d,
               trial_status = trial_status,
               arm_status = arm_status)
  
  lret
  

}








