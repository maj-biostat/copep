
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
  
  message(get_hash(), " cols in post matrix ", 
          paste0(colnames(post), collapse=" "))
  
  arms <- substr(colnames(post), 7, 7)
  arms <- as.numeric(arms)
  
  a_in_post <- logical(a_s$K)
  a_in_post[arms] <- T
  
  message(get_hash(), " arms derived from posterior draws ", 
          paste0(a_in_post, collapse=" "))
  a_in_post
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
decision_intrm <- function(post, a_s, t_s, lpar){
  
  if(is.null(post)){
    message(get_hash(), " decision_intrm post is null ")
    stop()
  }

  # Probability that arm K has the lowest (best) proportion of cases
  
  a_s$arms_in_post <- arm_in_post(post[, -ncol(post)], a_s)
  
  arms_to_check <- which((a_s$arms_in_post & a_s$active) == T)
  
  message(get_hash(), " computing p_best on arms ",
          paste0(arms_to_check, collapse = " "))
  
  a_s$p_best <- rep(DUMVAL, a_s$K)
  a_s$var_k <- rep(DUMVAL, a_s$K)
  a_s$p_beat_soc <- rep(DUMVAL, a_s$K)
  a_s$p_equiv_soc <- rep(DUMVAL, a_s$K)
  
  if(length(arms_to_check) == 1 & arms_to_check[1] == 1){
    
    message(get_hash(), " arms_to_check is of length 1 casting vector to matrix.")
    m <- matrix(post[, arms_to_check], ncol = 1)
    
    a_s$p_best[arms_to_check] <- prob_min(-m)
    a_s$var_k[arms_to_check] <- apply(m, 2, var)
    a_s$p_beat_soc[arms_to_check] <- prob_lt(m)
    a_s$p_equiv_soc[arms_to_check] <- prob_equiv(m, lpar$eq_delta)
    
  } else {
    a_s$p_best[arms_to_check] <- prob_min(-post[, arms_to_check])
    a_s$var_k[arms_to_check] <- apply(post[, arms_to_check], 2, var)
    a_s$p_beat_soc[arms_to_check] <- prob_lt(post[, arms_to_check])
    a_s$p_equiv_soc[arms_to_check] <- prob_equiv(post[, arms_to_check], lpar$eq_delta)
    
    
    # diff1_2 <- post[, 1]-post[, 2]
    # qtmp <- quantile(diff1_2, probs = c(0.025, 0.5, 0.975))
    # message(get_hash(), " quantiles of difference between arm 1 and arm 2 ", 
    #         paste0(round(qtmp, 3), collapse = " "))
    
  }
  

  
  message(get_hash(), " decision_intrm active ", 
          paste0(a_s$active, collapse = " "), " p_best ",
          paste0(round(a_s$p_best, 3), collapse = " "))
  
  message(get_hash(), " decision_intrm active ", 
          paste0(a_s$active, collapse = " "), " p_beat_soc ", 
          paste0(round(a_s$p_beat_soc, 3), collapse = " "))
  
  
  
  
  
  
  res_win <- which(a_s$p_beat_soc > lpar$thresh_sup)
  if(length(res_win)!=0){

    a_s$is_sup[res_win] <- T
    a_s$eff_at[res_win] <- a_s$nk[res_win]
    
    t_s$win <- T
    t_s$eff <- T
    t_s$arm <- T
    t_s$pb <- picked_best_arm(a_s, lpar)
    t_s$stop <- T
    
    message(get_hash(), 
            " decision_intrm, effective trt arm(s) ", 
            paste0(res_win, collapse = " "), " threshold ", 
            lpar$thresh_sup, " p_best ", 
            paste0(a_s$p_best, collapse = " "), " p_beat_soc ", 
            paste0(a_s$p_beat_soc, collapse = " "))
    
    return(list(trial_status = t_s,
                arm_status = a_s))
    
  }

  # Turn off poor arms (cannot currently turn off control)
  res_off <- which(a_s$p_beat_soc[-1] < lpar$thresh_fut & 
                     a_s$p_beat_soc[-1] != DUMVAL)
  
  if(length(res_off)!=0){
    
    # res_off index will always be off by 1 since we excluded the control
    res_off <- res_off + 1
    
    a_s$fut_at[res_off] <- a_s$nk[res_off]
    a_s$is_fut[res_off] <- T
    a_s$active[res_off] <- F
    a_s$p_rand[res_off] <- 0
    
    t_s$arm <- T
    
    message(get_hash(), 
            " decision_intrm, stopping poorly performing arm(s) ",
            paste0(res_off, collapse = " "), " with threshold ", 
            lpar$thresh_fut, " p_best ", 
            paste0(a_s$p_best, collapse = " "), " p_beat_soc ", 
            paste0(a_s$p_beat_soc, collapse = " "))

    if(sum(a_s$active) == 1){
      
      message(get_hash(), 
              " decision_intrm, stop trial, all trt arms now inactive, p_best ", 
              paste0(a_s$p_best, collapse = " "), " p_beat_soc ", 
              paste0(a_s$p_beat_soc, collapse = " "))
 
      t_s$win <- F
      t_s$fut_or_equ <- T
      t_s$pb <- picked_best_arm(a_s, lpar)
      t_s$stop <- T

    } 
    
    return(list(trial_status = t_s,
                arm_status = a_s))
 
  }
  
  
  
  
  # Identify clear winning ARM (PLURAL!!) (based on the definition of the test, 
  # we CANNOT stop if we think control is best)
  
  res_equiv <- which(a_s$p_equiv_soc[-1] > lpar$thresh_equ)

  if(length(res_equiv)!=0){
    
    res_equiv <- res_equiv + 1
    
    message(get_hash(), 
            " decision_intrm, equivalent trt arm(s) ", 
            paste0(res_equiv, collapse = " "), " p_equiv_soc ", 
            paste0(a_s$p_equiv_soc, collapse = " "))
    
    for(l in 1:length(res_equiv)){
      
      qtmp <- quantile(expit(post[, 1]) - expit(post[, res_equiv[l]]), 
                       probs = c(0.025, 0.5, 0.975))
      
      message(get_hash(), 
              " decision_intrm, equivalent trt arm(s) diff in p_est SOC - Arm_", l+1, " 95CI ", 
              paste0(round(qtmp, 3), collapse = " "))
    }

   
    a_s$equ_at[res_equiv] <- a_s$nk[res_equiv]
    a_s$is_equ[res_equiv] <- T
    a_s$active[res_equiv] <- F
    a_s$p_rand[res_equiv] <- 0
    
    t_s$arm <- T
    
    message(get_hash(), " decision_intrm active ", 
            paste0(a_s$active, collapse = " "), " p_equiv_soc ", 
            paste0(round(a_s$p_equiv_soc, 3), collapse = " "), " res_equiv ", 
            paste0(res_equiv, collapse = " "), " thresh ", 
            lpar$thresh_equ)
    
    if(sum(a_s$active) == 1){
      
      message(get_hash(), 
              " decision_intrm, stop trial, all trt arms now inactive, p_best ", 
              paste0(a_s$p_best, collapse = " "), " p_beat_soc ", 
              paste0(a_s$p_beat_soc, collapse = " "))
      
      t_s$win <- F
      t_s$fut_or_equ <- T
      t_s$pb <- picked_best_arm(a_s, lpar)
      t_s$stop <- T
      
    } 
    
  }
  
  
  

  
  return(list(trial_status = t_s,
              arm_status = a_s))
}


#' Assesses posterior from final analysis. In a fixed-size trial, this is the
#' only assessment that will be done.
#'
#' @param post 
#' @param lpar 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
decision_fa <- function(post, a_s, t_s, lpar){
  
  if(is.null(post)){
    message(get_hash(), "decision_fa post is null ")
  }
  
  
  # Probability that arm K has the lowest (best) proportion of cases
  
  a_s$arms_in_post <- arm_in_post(post[, -ncol(post)], a_s)
  
  arms_to_check <- which((a_s$arms_in_post & a_s$active) == T)
  
  message(get_hash(), " decision_fa computing p_best on arms ",
          paste0(arms_to_check, collapse = " "))
  
  a_s$p_best <- rep(DUMVAL, a_s$K)
  a_s$var_k <- rep(DUMVAL, a_s$K)
  a_s$p_beat_soc <- rep(DUMVAL, a_s$K)
  
  a_s$p_best[arms_to_check] <- prob_min(-post[, arms_to_check])
  a_s$var_k[arms_to_check] <- apply(post[, arms_to_check], 2, var)
  a_s$p_beat_soc[arms_to_check] <- prob_lt(post[, arms_to_check])
  
  
  message(get_hash(), " decision_fa active ", 
          paste0(a_s$active, collapse = " "), " p_best ",
          paste0(round(a_s$p_best, 3), collapse = " "))
  
  message(get_hash(), " decision_fa active ", 
          paste0(a_s$active, collapse = " "), " p_beat_soc ", 
          paste0(round(a_s$p_beat_soc, 3), collapse = " "))
  
  message(get_hash(), " decision_fa active ", 
          paste0(a_s$active, collapse = " "), " p_equiv_soc ", 
          paste0(round(a_s$p_equiv_soc, 3), collapse = " "))
  
  
  res_win <- which(a_s$p_beat_soc > lpar$thresh_sup)
  if(length(res_win)!=0){
  
    a_s$is_sup[res_win] <- T
    a_s$eff_at[res_win] <- a_s$nk[res_win] 
    
    t_s$win <- T
 
  } else {
    
    t_s$win <- F
    
  }
  
  res_off <- which(a_s$p_beat_soc[-1] < lpar$thresh_fut)
  if(length(res_off)!=0){
    
    res_off <- res_off + 1
    
    a_s$is_fut[res_off] <- T
    a_s$fut_at[res_off] <- a_s$nk[res_off] 
    
  } 
  

  res_equiv <- which(a_s$p_equiv_soc[-1] > lpar$thresh_equ)
  
  if(length(res_equiv)!=0){
    
    res_equiv <- res_equiv + 1
    
    a_s$is_equ[res_equiv] <- T
    a_s$equ_at[res_equiv] <- a_s$nk[res_equiv] 
 
  }
 
  t_s$pb <- picked_best_arm(a_s, lpar)
  
  message(get_hash(), " decision_fa, win ", t_s$win)
  
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
  
  min_arm <- which.min(lpar$prob_symp)
  
  if(a_s$is_sup[min_arm] == T){
    return(T)
  }
  
  return(F)
}
