

#' Balanced cluster randomisation across all arms
#'
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
BR_alloc <- function(ldat, idxstart, idxend, interim_idx){


  # Number arms
  K <- ldat$arm_status$K

  # # Number of clusters
  J <- length(unique(clustid[idxstart:idxend]))
  
  # Say if we have four arms 1 2 3 4,
  # and arms 1 and 3 are active, ie ldat$arm_status$active = T F T F,
  # then the active_arms variable will contain 1 3
  active_arms <- (1:K)[ldat$arm_status$active]
  

  # Add the randomised treatment allocation.
  # Enforce allocation of first K clusts to set of distinct arms.
  if(interim_idx == 1){
    
    message(get_hash(), " BR_alloc first interim enforcing all active arms to have some allocation ", 
            paste0(active_arms, collapse = " "))

    rand_arm <- c(sample(active_arms, size = length(active_arms), replace = F),
                  sample(active_arms, size = J-length(active_arms), replace = T))
    
  } else {
    
    message(get_hash(), " BR_alloc using active arms ")
    
    rand_arm <- c(sample(active_arms, size = J, replace = T))
  }
  
  rep(rand_arm, times = table(clustid[idxstart:idxend]))

}



#' Response adaptive randomisation
#'
#' @param d 
#' @param lpar 
#'
#' @return
#' @export
#'   
#' @examples
#' 
#' d, pbest, var_k, Jprev, lpar
RAR_alloc <- function(clustid, a_s, idxstart, idxend, interim_idx){
  
  
  # # Number of clusters
  J <- length(unique(clustid[idxstart:idxend]))
  
  # Say if we have four arms 1 2 3 4,
  # and arms 1 and 3 are active, ie ldat$arm_status$active = T F T F,
  # then the active_arms variable will contain 1 3
  
  arms_for_p_best <- which((a_s$arms_in_post & a_s$active) == T)
  
  # An arm might be active but this might be the first time it is randomised
  # to so it will not be in the posterior yet.
  active_arms <- (1:a_s$K)[a_s$active]
  
  a_s$p_rand <- rep(0, a_s$K)
  
  
  if(interim_idx == 1){
    
    message(get_hash(), " first interim balanced alloc total arms ", 
            a_s$K, " active_arms ",
            paste0(active_arms, collapse = " "))
    
    rand_arm <- c(sample(active_arms, size = length(active_arms), replace = F),
                  sample(active_arms, size = J-length(active_arms), replace = T))
    rand_arm <- rep(rand_arm, times = table(clustid[idxstart:idxend]))
    a_s$p_rand[active_arms] <- 1/length(active_arms)
    
    return(list(rand_arm = rand_arm, arm_status = a_s))
    
  } else {
    
    Jprev <- clustid[idxstart-1]
    pbest <- numeric(a_s$K)
    vark <- numeric(a_s$K)
    pbest[arms_for_p_best] <- a_s$p_best[arms_for_p_best]
    vark[arms_for_p_best] <- a_s$var_k[arms_for_p_best]
    # r <- sqrt(pbest * vark / Jprev) 
    
    r <- sqrt(pbest * vark / (a_s$nk+1) )
    
    message(get_hash(), " arms_enabled_for_anly ", 
            paste0(a_s$enabled_for_anly, collapse = " "), " a_s$active ",
            paste0(a_s$active, collapse = " "))

    # This deals with the situation where we have an activated arm but it was not included
    # in the previous analysis and so does not appear in posterior. This could happen because
    # the outcome is not available immediately. Thus even though we gaurantee that all arms
    # are randomised during activation, it is possible that the new arm contains no observed data.
    apply_balanced_alloc <- rep(F, a_s$K)
    apply_balanced_alloc[interim_idx >= a_s$enabled_for_anly] <- T
    apply_balanced_alloc[!a_s$active] <- F
    apply_balanced_alloc[arms_for_p_best] <- F
    

    mass_remaining <- 1
    if(sum(apply_balanced_alloc)>0){

      message(get_hash(), " apply_balanced_alloc to arms ",
              paste0(apply_balanced_alloc, collapse = " "))

      a_s$p_rand[apply_balanced_alloc] <- 1/sum(a_s$active)

      mass_remaining <- 1 - sum(a_s$p_rand[apply_balanced_alloc])
    }
    


    a_s$p_rand[arms_for_p_best] <- (r/sum(r))[arms_for_p_best] * mass_remaining
    
    message(get_hash(), " rand probs p_rand ", 
            paste0(round(a_s$p_rand, 3), collapse = " "))
    message(get_hash(), " pbest ", 
            paste0(round(pbest, 3), collapse = " "))
    message(get_hash(), " apply_balanced_alloc ", 
            paste0(apply_balanced_alloc, collapse = " "))
    
    
    if(sum(apply_balanced_alloc)>0){
      
      message(get_hash(), " apply_balanced_alloc ", 
              paste0(apply_balanced_alloc, collapse = " ") , 
              " as only arms ",
              paste0(a_s$arms_in_post, collapse = " "),
              " have post draws")
      
      # Note that in the first case p_rand is irrelevant and in the 
      # second case p_rand = 0 takes care of inactive arms.
      rand_arm <- c(sample(active_arms, 
                           size = length(active_arms), 
                           replace = F),
                    sample(1:a_s$K, 
                           size = J-length(active_arms), 
                           replace = T, 
                           prob = a_s$p_rand))
    }
    else {
      
      rand_arm <- c(sample(1:a_s$K, 
                           size = J, 
                           replace = T, 
                           prob = a_s$p_rand))
    }
 
    

    tbl <- table(rand_arm)
    message(get_hash(), " allocated ", 
            paste0(as.numeric(tbl), collapse = " "), " to active arms ", 
            paste0(active_arms, collapse = " "))

    

    return(list(rand_arm = rep(rand_arm, 
                               times = table(clustid[idxstart:idxend])), 
                arm_status = a_s))
  }


}
