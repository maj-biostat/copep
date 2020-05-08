
source("util.R")




#' Parses the column names from the posterior matrix to deduce which
#' arms are present.
#'
#' @param post 
#' @param a_s 
#'
#' @return
#' @export
#'
#' @examples
arm_in_post <- function(post, a_s){
  
  message(get_hash(), " cols in posterior matrix ", 
          paste0(colnames(post), collapse=" "))
  
  arms <- substr(colnames(post), 7, 7)
  arms <- as.numeric(arms)
  
  a_in_post <- logical(a_s$K)
  a_in_post[arms] <- T
  
  message(get_hash(), " found arms in posterior draws ", 
          paste0(a_in_post, collapse=" "))
  a_in_post
}


#' Computes the probabilities associated with the hypotheses of interest,
#' pbest, psup, pinf, pequiv
#'
#' @param post 
#' @param a_s 
#' @param arms_in_post 
#' @param arms_to_check 
#'
#' @return
#' @export
#'
#' @examples
compute_decision_probs <- function(
  post, 
  a_s,
  eq_delta){
  
  arms_in_post <- which(a_s$arms_in_post)
  arms_to_check <- which((a_s$arms_in_post & a_s$active) == T)
  arms_to_check <- arms_to_check[arms_to_check != 1]
  
  message(get_hash(), " compute_decision_probs arms_in_post ",
          paste0(arms_in_post, collapse = " "), " arms_to_check ",
          paste0(arms_to_check, collapse = " "))
  
  a_s$p_best <- rep(DUMVAL, a_s$K)
  a_s$p_beat_soc <- rep(DUMVAL, a_s$K)
  a_s$p_equiv_soc <- rep(DUMVAL, a_s$K)
  a_s$var_k <- rep(DUMVAL, a_s$K)
  a_s$is_best <- rep(F, a_s$K)
  
  # Compute all the necessary probabilities.
  if(length(arms_to_check) == 0){
    
    message(get_hash(), " arms_to_check is of length 0, i.e only soc active.")

    # p_best always based on all arms in posterior - only used for tracking.
    a_s$p_best[arms_in_post] <- prob_min(-post[, -ncol(post)])
    bestidx <- which.max(a_s$p_best)
    a_s$is_best[bestidx] <- T
    
    # reqd for randomisation and decisions
    a_s$var_k[c(1, arms_to_check)] <- var(post[, 1])
    a_s$p_beat_soc[c(1, arms_to_check)] <- 0
    a_s$p_equiv_soc[c(1, arms_to_check)] <- 1
    
  } else {
    # p_best always based on all arms in posterior - only used for tracking.
    a_s$p_best[arms_in_post] <- prob_min(-post[, -ncol(post)])
    bestidx <- which.max(a_s$p_best)
    a_s$is_best[bestidx] <- T
    
    # reqd for randomisation and decisions
    a_s$var_k[c(1, arms_to_check)] <- apply(post[, c(1, arms_to_check)], 2, var)
    a_s$p_beat_soc[c(1, arms_to_check)] <- prob_lt(post[, 1], 
                                             post[, arms_to_check])
    a_s$p_equiv_soc[c(1, arms_to_check)] <- prob_equiv(post[, 1], 
                                                 post[, arms_to_check], 
                                                 eq_delta)
  }
  
  return(a_s)
}


#' Checks for arms being superior to soc
#'
#' @return
#' @export
#'
#' @examples
effective_v_soc <- function(a_s, t_s, lpar){
  
  res_win <- which(a_s$p_beat_soc > lpar$thresh_sup)
  if(length(res_win)!=0){
    
    message(get_hash(), 
            "   EFFECTIVE arm(s), deactivating soc, res_win is ",
            paste0(res_win, collapse = " "), " threshold ", 
            lpar$thresh_sup, " p_best ", 
            paste0(round(a_s$p_best,3), collapse = " "), " p_beat_soc ", 
            paste0(round(a_s$p_beat_soc,3), collapse = " "))
    
    # soc arm receives no allocation from here on in - deactivate soc
    a_s$active[1] <- F
    # a_s$is_inf[1] <- T
    # If we haven't previously found a sup arm then set the sample size, 
    # otherwise leave as it.
    if(is.na(a_s$inf_at[1])){
      a_s$inf_at[1] <- a_s$nk[1]
    }

    a_s$is_sup[res_win] <- T
    a_s$sup_at[res_win] <- a_s$nk[res_win]
    
    # Can consider trial as a success (counts towards power)
    if(is.na(t_s$sup)){
      t_s$sup <- T
      t_s$inf_or_equ <- F
    }
    
    t_s$pb <- picked_best_arm(a_s, lpar)
    
    message(get_hash(), 
            " active arms ", 
            paste0(a_s$active, collapse = " "))
    
    if(sum(a_s$active) <= 1){
      message(get_hash(), 
              " stop trial, found effective, no further arms to compare")
      t_s$stop <- T
    }

    
    message(get_hash(), 
            " effective arms ", 
            paste0(a_s$is_sup, collapse = " "))
    
    if(a_s$active[1] == F & 
       sum(a_s$is_sup[2:length(a_s$is_sup)], na.rm = T) == length(a_s$is_sup)-1){
      message(get_hash(), 
              " stop trial, all arms effective, no further arms to compare")
      t_s$stop <- T
    }
    
  }
  tmp <- list(a_s = a_s, t_s = t_s)
  return(tmp)
}


#' Inferiority check relative to soc
#'
#' @param a_s 
#' @param t_s 
#' @param lpar 
#'
#' @return
#' @export
#'
#' @examples
ineffective_v_soc <- function(a_s, t_s, lpar){
  
  res_off <- which(a_s$p_beat_soc[-1] < lpar$thresh_fut & 
                     a_s$p_beat_soc[-1] != DUMVAL)
  
  if(length(res_off)!=0){
    
    # res_off index will always be off by 1 since we excluded the control
    res_off <- res_off + 1
    
    for(j in res_off){
      if(is.na(a_s$inf_at[j])){
        a_s$inf_at[j] <- a_s$nk[j]
      }
    }

    a_s$is_inf[res_off] <- T
    a_s$active[res_off] <- F

    message(get_hash(), 
            "   INEFFECTIVE arm(s) ",
            paste0(res_off, collapse = " "), " with threshold ", 
            lpar$thresh_fut, " p_best ", 
            paste0(round(a_s$p_best,3), collapse = " "), " p_beat_soc ", 
            paste0(round(a_s$p_beat_soc,3), collapse = " "))
    
    if(is.na(t_s$sup)){
      message(get_hash(), 
              " set current trial status to ineffective")
      t_s$inf_or_equ <- T
    }
    
    message(get_hash(), 
            " active arms ", 
            paste0(a_s$active, collapse = " "))
    
    if(sum(a_s$active) <= 1){
      
      message(get_hash(), 
              " stopping trial, ineffective, no further arms to compare")
      
      # If no arms have been identified as being superior then consider all
      # treatments inferior or equivalent and trial is considered failure.
      
      # Note that this prevents reversing an initial decision that an arm
      # was superior to soc.
      
      
      t_s$pb <- picked_best_arm(a_s, lpar)
      t_s$stop <- T
 
    } 
  }
  tmp <- list(a_s = a_s, t_s = t_s)
  return(tmp)
  
}


#' Equivalent check relative to soc
#'
#' @param a_s 
#' @param t_s 
#' @param lpar 
#'
#' @return
#' @export
#'
#' @examples
equivalent_v_soc <- function(a_s, t_s, lpar){
  
  res_equiv <- which(a_s$p_equiv_soc[-1] > lpar$thresh_equ)
  if(length(res_equiv)!=0){
    
    res_equiv <- res_equiv + 1
    
    message(get_hash(), 
            "   EQUIVALENT arm(s) ", 
            paste0(res_equiv, collapse = " "), " p_equiv_soc ", 
            paste0(round(a_s$p_equiv_soc,3), collapse = " "), " res_equiv ", 
            paste0(res_equiv, collapse = " "), " thresh ", 
            lpar$thresh_equ)
    
    
    for(j in res_equiv){
      if(is.na(a_s$equ_at[j])){
        a_s$equ_at[j] <- a_s$nk[j]
      }
    }
    
    a_s$is_equ[res_equiv] <- T
    a_s$active[res_equiv] <- F
    
    if(is.na(t_s$sup)){
      
      message(get_hash(), 
              " set current trial status to equivalent")
      t_s$inf_or_equ <- T
    }
    
    message(get_hash(), 
            " active arms ", 
            paste0(a_s$active, collapse = " "))
    
    if(sum(a_s$active) <= 1){
 
      # If no arms have been identified as being superior then consider all
      # treatments inferior or equivalent and trial is considered failure.
      
      # Note that this prevents reversing an initial decision that an arm
      # was superior to soc.
      message(get_hash(), 
                " stop trial, equivalent, no further arms to compare")
   
      t_s$pb <- picked_best_arm(a_s, lpar)
      t_s$stop <- T

    } 
  }
  tmp <- list(a_s = a_s, t_s = t_s)
  return(tmp)
  
}


#' Check on whether active arms inferior to arm deemed sup to soc
#'
#' @param a_s 
#' @param t_s 
#' @param lpar 
#'
#' @return
#' @export
#'
#' @examples
ineffective_v_effective <- function(post, a_s, t_s, lpar){
  
  # Inferior check relative to arms that are deemed sup to soc
  if(!is.na(t_s$sup) & t_s$sup == T){
    
    arms_sup <- which(a_s$is_sup)
    arms_to_check <- which((a_s$arms_in_post & a_s$active) == T)
    
    # Usually res_win is only length 1 but it is possible to find two winning arms.
    for(isup in arms_sup){
      
      # For each arm not in the set {K_1, K_sup} where K_1 is the soc arm and K_sup
      for(j in arms_to_check){
        
        if(a_s$active[j]==F){
          next
        }
        if(j %in% arms_sup){
          next
        }
        if(a_s$arms_in_post[j]==F){
          next
        }
        
        # If we have arm(s) deemed superior to soc (coz they have high prob of lower
        # response) then assess whether the sup arms is also sup to the other remaining
        # active arms. If so then we deactive the relevant arm.
        if(mean(post[, isup] < post[, j]) > lpar$thresh_sup){
          
          message(get_hash(), 
                  "   INFERIOR arm ", j, 
                  " inferior to effective arm ", isup)
          
          a_s$inf_at[j] <- a_s$nk[j]
          a_s$is_inf[j] <- T
          a_s$active[j] <- F
          
          message(get_hash(), 
                  " active arms ", 
                  paste0(a_s$active, collapse = " "))
          
          # And if we have only got 1 arm left then stop trial
          if(sum(a_s$active) <= 1){
            
            message(get_hash(), 
                    " stop trial, no further arms to compare with res_win.")
            
            # Other stuff (sup, inf etc) already set.
            
            t_s$pb <- picked_best_arm(a_s, lpar)
            t_s$stop <- T

          } 
        }
      }
    }
  }
  
  tmp <- list(a_s = a_s, t_s = t_s)
  return(tmp)
  
}


#' Decision logic for interim analysis. Determines whether a stopping rule
#' has been triggered.
#'
#' @param post 
#' @param lpar 
#'
#' @return
#' @export
#'
#' @examples
decision <- function(post, a_s, t_s, lpar){
  
  message(get_hash(), " decision called ")
  if(is.null(post)){
    message(get_hash(), " decision post is null ")
    stop()
  }

  # Probability that arm K has the lowest (best) proportion of cases
  a_s$arms_in_post <- arm_in_post(post[, -ncol(post)], a_s)
  
  a_s <- compute_decision_probs(post, a_s, lpar$eq_delta)

  message(get_hash(), " active arms ", 
          paste0(a_s$active, collapse = " "), " p_beat_soc ", 
          paste0(round(a_s$p_beat_soc, 3), collapse = " "), " p_best ",
          paste0(round(a_s$p_best, 3), collapse = " "))

  # Effectiveness test
  tmp <- effective_v_soc(a_s, t_s, lpar)
  a_s <- tmp$a_s
  t_s <- tmp$t_s
  
  # Ineffectiveness test (relative to soc)
  tmp <- ineffective_v_soc(a_s, t_s, lpar)
  a_s <- tmp$a_s
  t_s <- tmp$t_s

  # Equivalence test (relative to soc)
  tmp <- equivalent_v_soc(a_s, t_s, lpar)
  a_s <- tmp$a_s
  t_s <- tmp$t_s
  
  # Inferiority check (relative to arms deemed effective v soc)
  tmp <- ineffective_v_effective(post, a_s, t_s, lpar)
  a_s <- tmp$a_s
  t_s <- tmp$t_s
  
  return(list(trial_status = t_s,
              arm_status = a_s))
}






picked_best_arm <- function(a_s, lpar){
  
  all_arms_equal <- T
  first_arm <- lpar$prob_symp[1]
  for(i in 2:length(lpar$prob_symp)){
    if(first_arm != lpar$prob_symp[i]){
      all_arms_equal <- F
      break
    }
  }

  if(all_arms_equal){
    return(F)
  }
  
  # return first value which is min
  min_arm <- which.min(lpar$prob_symp)
  
  if(!is.na(a_s$is_sup[min_arm]) & a_s$is_sup[min_arm] == T){
    return(T)
  }
  
  return(F)
}
