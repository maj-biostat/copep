# !diagnostics off


library(tidyverse)
library(randomizr)
library(extraDistr)
library(tictoc)
library(rstan)
# Try jags
library(rjags)
library(coda)
# library(qs)

source("data.R")
source("util.R")
source("models.R")
source("randomise.R")
source("trials.R")
source("decision.R")



# test_RAR_alloc()
# test_decision_intrm()
# test_picked_best()
# test_prob_min()
# test_prob_lt()


post_samples <- function(){
  
  lpost <- list()
  
  set.seed(1)
  # Note that max.col used by prob_min uses random as the ties method
  # and therefore we need quite a lot of samples before the prob_min are
  # approx all equal.
  a1 <- rnorm(10000, qlogis(0.1), 1)
  a2 <- rnorm(10000, qlogis(0.1), 1)
  a3 <- rnorm(10000, qlogis(0.1), 1)
  a4 <- rnorm(10000, qlogis(0.1), 1)
  junk <- rnorm(10000, qlogis(0.1), 1)
  
  lpost$post3_all_equal <- cbind(a1, a2, a3, junk)
  lpost$post4_all_equal <- cbind(a1, a2, a3, a4, junk)
  
  set.seed(1)
  a1 <- rnorm(1000, qlogis(0.1), 1)
  a2 <- rnorm(1000, qlogis(0.1), 1)
  a3 <- rnorm(1000, qlogis(0.1), 1)
  a4 <- rnorm(1000, qlogis(0.1), 1)
  junk <- rnorm(1000, qlogis(0.1), 1)
  
  lpost$post3_1_best <- cbind(a1-1.5, a2, a3, junk)
  lpost$post3_3_best <- cbind(a1, a2, a3-1.5, junk)
  
  lpost$post3_2_3_best <- cbind(a1, a2-1.5, a3-1.5, junk)
  
  lpost$post3_1_worst <- cbind(a1+1.3, a2, a3, junk)
  lpost$post3_2_worst <- cbind(a1, a2+1.3, a3, junk)
  lpost$post3_2_really_worse <- cbind(a1, a2+1.5, a3, junk)
  
  lpost$post3_2_3_worst <- cbind(a1, a2+1, a3+1, junk)
  
  
  lpost$post4_1_best <- cbind(a1-2.2, a2, a3, a4, junk)
  lpost$post4_2_best <- cbind(a1, a2-2.2, a3, a4, junk)
  
  lpost$post4_2_4_best <- cbind(a1, a2-2.2, a3, a4-2.2, junk)
  
  lpost$post4_1_worst <- cbind(a1+1.5, a2, a3, a4, junk)
  lpost$post4_3_worst <- cbind(a1, a2, a3+1.5, a4, junk)
  lpost$post4_4_worst <- cbind(a1, a2, a3, a4+1.5, junk)
  
  # prob_min(-lpost$post3_2_really_worse[, 1:2])
  # prob_lt(lpost$post3_2_really_worse[, 1:2])
  
  lpost <- lapply(lpost, function(x){
    colnames(x) <- c(paste0("b.arm_", 1:(ncol(x)-1)), "sigu")
    x
  })
  
  lpost
}


test_RAR_alloc <- function(){
  
  # Balanced randomisation interim 1
  message("----> TESTING RAR_alloc")
  
  
  message("----> TEST Randomisation during first interim")
  a_s <- list(K = 4, 
              active = c(T, T, T, F),
              arms_enabled_for_anly = c(1, 1, 1, 2),
              arms_in_post = c(T, T, T, F),
              p_rand = rep(0, 4),
              p_best = c(0.1, 0.2, 0.4, 0),
              p_beat_soc = c(0, 0.7, 0.5, DUMVAL),
              var_k = c(0.2, 0.3, 0.3, DUMVAL))
  
  d <- tibble(clustid = rep(1:300, each = 2))


  # Indexes are participant level.
  idxstart = 1; idxend = 100; interim_idx = 1
  
  
  res <- RAR_alloc(d$clustid, 
                   a_s, 
                   idxstart, idxend, interim_idx)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  # Shouldn't get updated in this case
  stopifnot(res$arm_status$p_best == a_s$p_best)
  # Cannot compare NA
  stopifnot(res$arm_status$p_beat_soc[!is.na(res$arm_status$p_beat_soc)] == 
              a_s$p_beat_soc[!is.na(res$arm_status$p_beat_soc)])
  stopifnot(res$arm_status$var_k[!is.na(res$arm_status$var_k)] == 
              a_s$var_k[!is.na(res$arm_status$var_k)])
  stopifnot(res$arm_status$p_rand == c(1/3, 1/3, 1/3, 0))
  
  stopifnot(1:3 %in% res$rand_arm)
  stopifnot(!(4 %in% res$rand_arm))

  # 
  message("----> TEST Randomisation during activation of arm 4 with all other arms active")
  nk <- c(38,30,32, 0)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0.4, 0),
                   p_beat_soc = c(0, 0.7, 0.5, DUMVAL),
                   var_k = c(0.2, 0.3, 0.3, DUMVAL),
                   nk = nk)

  # Indexes are participant level.
  idxstart = 101; idxend = 120; interim_idx = 2

  res <- RAR_alloc(clustid=d$clustid, 
                   a_s, 
                   idxstart, idxend, interim_idx)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  
  r <- sqrt(a_s$p_best[1:3] * a_s$var_k[1:3] / (nk[1:3]+1))
  p_rand = numeric(4)
  p_rand[1:3] <- r/sum(r) * 0.75
  p_rand[4] <- 1/4

  stopifnot(abs(res$arm_status$p_rand - p_rand) < 0.0001)
  
  stopifnot(res$arm_status$p_beat_soc[!is.na(res$arm_status$p_beat_soc)] == 
              a_s$p_beat_soc[!is.na(res$arm_status$p_beat_soc)])
  stopifnot(res$arm_status$var_k[!is.na(res$arm_status$var_k)] == 
              a_s$var_k[!is.na(res$arm_status$var_k)])
  
  stopifnot(1:4 %in% res$rand_arm)
  

  # 
  message("----> TEST Randomisation during activation of arm 4 with one of the other arms deactivated")
  nk <- c(38,30,32, 0)
  a_s <- list(K = 4, 
                   active = c(T, F, T, T),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, DUMVAL, 0.4, DUMVAL),
                   p_beat_soc = c(0, DUMVAL, 0.5, DUMVAL),
                   var_k = c(0.2, DUMVAL, 0.3, DUMVAL),
                   nk = nk)

  # Indexes are participant level.
  idxstart = 101; idxend = 120; interim_idx = 2
  
  res <- RAR_alloc(clustid=d$clustid, 
                   a_s=a_s, 
                   idxstart, idxend, interim_idx)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  
  r <- sqrt(a_s$p_best[c(1, 3)] * a_s$var_k[c(1, 3)] / (nk[c(1, 3)]+1))
  p_rand = numeric(4)
  p_rand[c(1, 3)] <- r/sum(r) * 2/3
  p_rand[2] <- 0
  p_rand[4] <- 1/3
  
  stopifnot(abs(res$arm_status$p_rand - p_rand) < 0.0001)
  
  stopifnot(res$arm_status$p_beat_soc[!is.na(res$arm_status$p_beat_soc)] == 
              a_s$p_beat_soc[!is.na(res$arm_status$p_beat_soc)])
  stopifnot(res$arm_status$var_k[!is.na(res$arm_status$var_k)] == 
              a_s$var_k[!is.na(res$arm_status$var_k)])
  
  stopifnot(c(1, 3, 4) %in% res$rand_arm)
  stopifnot(!(2 %in% res$rand_arm))
  
  
  
  
  # 
  message("----> TEST Randomisation after activation of arm 4 with all other arms active and arm 4 in post")
  nk <- c(42,30,38, 10); sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0.4, 0.2),
                   p_beat_soc = c(0, 0.7, 0.5, 0.2),
                   var_k = c(0.2, 0.3, 0.3, 0.4),
              nk = nk)

  # Indexes are participant level.
  idxstart = 121; idxend = 160; interim_idx = 3
  
  # 
  set.seed(1)
  res <- RAR_alloc(d$clustid, a_s,  idxstart, idxend, interim_idx)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  
  # The 60 is the previous number of clusters 
  r <- sqrt(a_s$p_best[1:4] * a_s$var_k[1:4] / (nk+1))
  p_rand = numeric(4)
  p_rand[1:4] <- r/sum(r) 
  
  stopifnot(abs(res$arm_status$p_rand - p_rand) < 0.0001)
  
  stopifnot(res$arm_status$p_beat_soc == a_s$p_beat_soc)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  # Occasionally this might fail because there is no guarantee that all arms
  # will be randomised during a given interim especially if pbest is tiny for
  # an arm. I have set the seed above so that it doesn't fail in the test 
  # environment.
  stopifnot(1:4 %in% res$rand_arm)
  # Now jumble the seed again.
  set.seed(as.numeric(Sys.time()))
  
  
  
  
  # 
  
  message("----> TEST Randomisation after activation of arm 4 with all other arms active and arm 4 ***NOT*** in post")
  nk <- c(42,30,38, 10); sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0.4, DUMVAL),
                   p_beat_soc = c(0, 0.7, 0.5, DUMVAL),
                   var_k = c(0.2, 0.3, 0.3, DUMVAL),
              nk = nk)

  # Indexes are participant level.
  idxstart = 121; idxend = 160; interim_idx = 3
  
  # 
  set.seed(1)
  res <- RAR_alloc(clustid = d$clustid, a_s,  idxstart, idxend, interim_idx)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  
  # The 60 is the previous number of clusters 
  r <- sqrt(a_s$p_best[1:3] * a_s$var_k[1:3] / (nk[1:3]+1))
  p_rand = numeric(4)
  p_rand[1:3] <- r/sum(r)  * 3/4
  p_rand[4] <- 1/4
  
  
  stopifnot(abs(res$arm_status$p_rand - p_rand) < 0.0001)
  
  stopifnot(res$arm_status$p_beat_soc[!is.na(res$arm_status$p_beat_soc)] == 
              a_s$p_beat_soc[!is.na(res$arm_status$p_beat_soc)])
  stopifnot(res$arm_status$var_k[!is.na(res$arm_status$var_k)] == 
              a_s$var_k[!is.na(res$arm_status$var_k)])
  
  # Occasionally this might fail because there is no guarantee that all arms
  # will be randomised during a given interim especially if pbest is tiny for
  # an arm. I have set the seed above so that it doesn't fail in the test 
  # environment.
  stopifnot(1:4 %in% res$rand_arm)
  # Now jumble the seed again.
  set.seed(as.numeric(Sys.time()))
  
  

  # 
  message("----> TEST Randomisation after activation of arm 4 with one of the other arms inactive and arm 4 ***NOT*** in post")
  nk <- c(42,30,38, 10); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, T),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.4, 0.5, DUMVAL, DUMVAL),
                   p_beat_soc = c(0, 0.6, DUMVAL, DUMVAL),
                   var_k = c(0.2, 0.3, DUMVAL, DUMVAL),
              nk = nk)

  # Indexes are participant level.
  idxstart = 121; idxend = 160; interim_idx = 3
  
  # 
  set.seed(1)
  res <- RAR_alloc(d$clustid, a_s,  idxstart, idxend, interim_idx)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  
  # The 60 is the previous number of clusters 
  r <- sqrt(a_s$p_best[1:2] * a_s$var_k[1:2] / (nk[1:2]+1))
  p_rand = numeric(4)
  p_rand[1:2] <- r/sum(r)  * 2/3
  p_rand[3] <- 0
  p_rand[4] <- 1/3
  
  
  stopifnot(abs(res$arm_status$p_rand - p_rand) < 0.0001)
  
  stopifnot(res$arm_status$p_beat_soc[!is.na(res$arm_status$p_beat_soc)] == 
              a_s$p_beat_soc[!is.na(res$arm_status$p_beat_soc)])
  stopifnot(res$arm_status$var_k[!is.na(res$arm_status$var_k)] == 
              a_s$var_k[!is.na(res$arm_status$var_k)])
  
  # Occasionally this might fail because there is no guarantee that all arms
  # will be randomised during a given interim especially if pbest is tiny for
  # an arm. I have set the seed above so that it doesn't fail in the test 
  # environment.
  stopifnot(c(1, 2, 4) %in% res$rand_arm)
  # Now jumble the seed again.
  set.seed(as.numeric(Sys.time()))
  
  
  # 
  message("----> TEST Randomisation during activation of arm 4 with all other trt arms inactive and only Soc and arm 2 have post")
  nk <- c(42,30,38, 10); # sum(nk)
  a_s <- list(K = 4, 
              active = c(T, F, F, T),
              enabled_for_anly = c(1, 1, 1, 3),
              arms_in_post = c(T, T, F, F),
              p_rand = rep(0, 4),
              p_best = c(0.4, DUMVAL, DUMVAL, DUMVAL),
              p_beat_soc = c(0, DUMVAL, DUMVAL, DUMVAL),
              var_k = c(0.2, DUMVAL, DUMVAL, DUMVAL),
              nk = nk)
  
  # Indexes are participant level.
  idxstart = 121; idxend = 160; interim_idx = 3
  
  # 
  set.seed(1)
  res <- RAR_alloc(d$clustid, a_s,  idxstart, idxend, interim_idx)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  
  p_rand <- c(1/2, 0, 0, 1/2)
  
  
  stopifnot(abs(res$arm_status$p_rand - p_rand) < 0.0001)
  
  stopifnot(res$arm_status$p_beat_soc[!is.na(res$arm_status$p_beat_soc)] == 
              a_s$p_beat_soc[!is.na(res$arm_status$p_beat_soc)])
  stopifnot(res$arm_status$var_k[!is.na(res$arm_status$var_k)] == 
              a_s$var_k[!is.na(res$arm_status$var_k)])
  
  # Occasionally this might fail because there is no guarantee that all arms
  # will be randomised during a given interim especially if pbest is tiny for
  # an arm. I have set the seed above so that it doesn't fail in the test 
  # environment.
  stopifnot(c(1, 4) %in% res$rand_arm)
  # Now jumble the seed again.
  set.seed(as.numeric(Sys.time()))
  
}





test_decision_intrm <- function(){
  
  message("----> TESTING decision_intrm")
  
  # post, a_s, t_s, lpar
  
  # posterior samples
  lpost <- post_samples()
  # decision_intrm <- function(post, cur_arm_status, lpar)
  lpar <- list()
  lpar$prob_symp = c(0.15, 0.15, 0.15, 0.15)
  lpar$thresh_sup <- 0.9
  lpar$thresh_fut <- 0.1

  message("\n\n TESTS ON 3 ARMS")
  
  # 
  message("----> TEST 1 Decision with 3 arms active prior to activating arm 4. All arms are equivalent.")
  nk <- c(36,30,34, 0); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, F),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0.4, DUMVAL),
                   p_beat_soc = c(0, 0.7, 0.5, DUMVAL),
                   is_sup = rep(F, 4),
                   is_fut = rep(F, 4),
                   var_k = c(0.2, 0.3, 0.3, DUMVAL),
              nk = nk)

  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 100  # current total number of clusters
                       )
  
  set.seed(1)
  a_s$p_best[1:3] <- prob_min(-lpost$post3_all_equal[, 1:3])
  set.seed(1)
  a_s$p_beat_soc[1:3] <- prob_lt(-lpost$post3_all_equal[, 1:3])
  a_s$var_k[1:3] <- apply(lpost$post3_all_equal[, 1:3], 2, var)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_all_equal, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
                 
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 

  stopifnot(res$arm_status$p_best[-4] == a_s$p_best[-4])
  stopifnot(abs(res$arm_status$p_beat_soc[-4] - a_s$p_beat_soc[-4])<0.1)
  stopifnot(res$arm_status$var_k[-4] == a_s$var_k[-4])
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(all(is.na(res$arm_status$fut_at)))

  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 100)
 
  # 
  message("----> TEST 2 Decision with 3 arms (arm 2 is inactive) prior to activating arm 4. All arms are equivalent.")
  nk <- c(36,30,34, 0); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, F, T, F),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, 10, NA, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, DUMVAL, 0.4, DUMVAL),
                   p_beat_soc = c(0, DUMVAL, 0.5, DUMVAL),
                   is_sup = rep(F, 4),
                   is_fut = c(F, T, F, F),
                   var_k = c(0.2, DUMVAL, 0.3, DUMVAL),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 100  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 3)] <- prob_min(-lpost$post3_all_equal[, c(1, 3)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 3)] <- prob_lt(-lpost$post3_all_equal[, c(1, 3)])
  a_s$var_k[c(1, 3)] <- apply(lpost$post3_all_equal[, c(1, 3)], 2, var)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_all_equal, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(res$arm_status$p_best[c(1, 3)] == a_s$p_best[c(1, 3)])
  stopifnot(abs(res$arm_status$p_beat_soc[c(1, 3)] - a_s$p_beat_soc[c(1, 3)])<0.1)
  stopifnot(res$arm_status$var_k[c(1, 3)] == a_s$var_k[c(1, 3)])
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[2]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 3, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 100) 
  
  
  # 
  message("----> TEST 3 Decision with 3 arms active prior to activating arm 4. SOC best")
  nk <- c(36,30,34, 0); # sum(nk)
  # The comparison is always to SOC so this should not update trial state.
  a_s <- list(K = 4, 
                   active = c(T, T, T, F),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0.4, DUMVAL),
                   p_beat_soc = c(0, 0.7, 0.5, DUMVAL),
                   is_sup = rep(F, 4),
                   is_fut = rep(F, 4),
                   var_k = c(0.2, 0.3, 0.3, DUMVAL),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 100  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:3] <- prob_min(-lpost$post3_1_best[, 1:3])
  set.seed(1)
  a_s$p_beat_soc[1:3] <- prob_lt(lpost$post3_1_best[, 1:3])
  a_s$var_k[1:3] <- apply(lpost$post3_1_best[, 1:3], 2, var)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_1_best, 
                        a_s = a_s, 
                        t_s = trial_status, 
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(res$arm_status$p_best[1:3] == a_s$p_best[1:3])
  stopifnot(abs(res$arm_status$p_beat_soc[1:3] - a_s$p_beat_soc[1:3])<0.1)
  stopifnot(res$arm_status$var_k[1:3] == a_s$var_k[1:3])
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(all(is.na(res$arm_status$fut_at)))
  
  # We don't stop if SOC is best arm because of the definition of the stopping rule.
  # The stopping rule tests whether an arm is better than SOC.
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 100) 
  
  
  
  # 
  message("----> TEST 4 Decision with 3 arms (all active) prior to activating arm 4. Arm 3 best")
  nk <- c(36,30,34, 0); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, F),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0.4, DUMVAL),
                   p_beat_soc = c(0, 0.7, 0.5, DUMVAL),
                   is_sup = rep(F, 4),
                   is_fut = rep(F, 4),
                   var_k = c(0.2, 0.3, 0.3, DUMVAL),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 100  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:3] <- prob_min(-lpost$post3_3_best[, 1:3])
  set.seed(1)
  a_s$p_beat_soc[1:3] <- prob_lt(lpost$post3_3_best[, 1:3])
  a_s$var_k[1:3] <- apply(lpost$post3_3_best[, 1:3], 2, var)
  
  lpar$thresh_sup <- 0.8
  lpar$prob_symp = c(0.15, 0.15, 0.1, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_3_best, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(res$arm_status$p_best[1:3] == a_s$p_best[1:3])
  stopifnot(abs(res$arm_status$p_beat_soc[1:3] - a_s$p_beat_soc[1:3])<0.1)
  stopifnot(res$arm_status$var_k[1:3] == a_s$var_k[1:3])
  
  stopifnot(res$arm_status$is_sup == c(F, F, T, F))
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(res$arm_status$eff_at[3] == nk[3])
  stopifnot(all(is.na(res$arm_status$eff_at[c(1, 2, 4)])))
  stopifnot(all(is.na(res$arm_status$fut_at)))

  stopifnot(res$trial_status$stop == T)
  stopifnot(res$trial_status$win == T)
  stopifnot(res$trial_status$eff == T)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == T)
  stopifnot(res$trial_status$nk == 100) 
  
  
  
  # 
  message("----> TEST 5 Decision with 3 arms (arm 2 inactive) prior to activating arm 4. Arm 3 best")
  nk <- c(36,30,34, 0); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, F, T, F),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, 10, NA, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, DUMVAL, 0.4, DUMVAL),
                   p_beat_soc = c(0, DUMVAL, 0.5, DUMVAL),
                   var_k = c(0.2, DUMVAL, 0.3, DUMVAL),
                   is_sup = c(F, F, T, F),
                   is_fut = c(F, T, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 100  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 3)] <- prob_min(-lpost$post3_3_best[, c(1, 3)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 3)] <- prob_lt(lpost$post3_3_best[, c(1, 3)])
  a_s$var_k[c(1, 3)] <- apply(lpost$post3_3_best[, c(1, 3)], 2, var)
  
  lpar$thresh_sup <- 0.8
  lpar$prob_symp = c(0.15, 0.15, 0.1, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_3_best, 
                        a_s = a_s, 
                        t_s = trial_status, 
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(res$arm_status$p_best[c(1, 3)] == a_s$p_best[c(1, 3)])
  stopifnot(abs(res$arm_status$p_beat_soc[c(1, 3)] - a_s$p_beat_soc[c(1, 3)])<0.1)
  stopifnot(res$arm_status$var_k[c(1, 3)] == a_s$var_k[c(1, 3)])
  
  stopifnot(res$arm_status$is_sup == c(F, F, T, F))
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(res$arm_status$eff_at[3]==nk[3])
  stopifnot(all(is.na(res$arm_status$eff_at[c(1, 2, 4)])))
  stopifnot(res$arm_status$fut_at[2]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 3, 4)])))

  stopifnot(res$trial_status$stop == T)
  stopifnot(res$trial_status$win == T)
  stopifnot(res$trial_status$eff == T)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == T)
  stopifnot(res$trial_status$nk == 100)  
  
  
  
  # 
  message("----> TEST 6 Decision with 3 arms (all active) prior to activating arm 4. Arm 1 worst")
  nk <- c(36,30,34, 0); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, F),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0.4, DUMVAL),
                   p_beat_soc = c(0, 0.2, 0.5, DUMVAL),
                   var_k = c(0.2, 0.2, 0.3, DUMVAL),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 100  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:3] <- prob_min(-lpost$post3_1_worst[, 1:3])
  set.seed(1)
  a_s$p_beat_soc[1:3] <- prob_lt(lpost$post3_1_worst[, 1:3])
  a_s$var_k[1:3] <- apply(lpost$post3_1_worst[, 1:3], 2, var)
  
  lpar$thresh_sup <- 0.9
  lpar$prob_symp = c(0.15, 0.11, 0.1, 0.11)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_1_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(res$arm_status$p_best[1:3] == a_s$p_best[1:3])
  stopifnot(abs(res$arm_status$p_beat_soc[1:3] - a_s$p_beat_soc[1:3])<0.01)
  stopifnot(res$arm_status$var_k[1:3] == a_s$var_k[1:3])
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(all(is.na(res$arm_status$fut_at)))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 100)  
  
  
  # 
  message("----> TEST 7 Decision with 3 arms (arm 3 inactive) prior to activating arm 4. Arm 1 worst")
  nk <- c(36,30,34, 0); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, F),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, NA, 10, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, DUMVAL, DUMVAL),
                   p_beat_soc = c(0, 0.2, DUMVAL, DUMVAL),
                   var_k = c(0.2, 0.2, DUMVAL, DUMVAL),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, T, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 100  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:2] <- prob_min(-lpost$post3_1_worst[, 1:2])
  set.seed(1)
  a_s$p_beat_soc[1:2] <- prob_lt(lpost$post3_1_worst[, 1:2])
  a_s$var_k[1:2] <- apply(lpost$post3_1_worst[, 1:2], 2, var)
  
  lpar$thresh_sup <- 0.9
  lpar$prob_symp = c(0.15, 0.11, 0.1, 0.11)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_1_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(res$arm_status$p_best[1:2] == a_s$p_best[1:2])
  stopifnot(abs(res$arm_status$p_beat_soc[1:2] - a_s$p_beat_soc[1:2])<0.1)
  stopifnot(res$arm_status$var_k[1:2] == a_s$var_k[1:2])
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[3]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 2, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 100)
  
  # 
  message("----> TEST 8 Decision with 3 arms (all active) prior to activating arm 4. Arm 2 worst")
  nk <- c(36,30,34, 0); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, F),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0.4, DUMVAL),
                   p_beat_soc = c(0, 0.2, 0.5, DUMVAL),
                   var_k = c(0.2, 0.2, 0.3, DUMVAL),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = T, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 100  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:3] <- prob_min(-lpost$post3_2_worst[, 1:3])
  set.seed(1)
  a_s$p_beat_soc[1:3] <- prob_lt(lpost$post3_2_worst[, 1:3])
  a_s$var_k[1:3] <- apply(lpost$post3_2_worst[, 1:3], 2, var)

  lpar$thresh_fut <- 0.2
  lpar$prob_symp = c(0.15, 0.2, 0.15, 0.15)

  set.seed(1)
  res <- decision_intrm(post = lpost$post3_2_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == c(T, F, T, F))
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == c(F, T, F, F))
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[2]==nk[2])
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 3, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 100)
  
  
  message("----> TEST 9 Decision with 3 arms (arm 3 inactive) prior to activating arm 4. Arm 2 worst")
  nk <- c(36,30,34, 0); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, F),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, NA, 10, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, DUMVAL, DUMVAL),
                   p_beat_soc = c(0, 0.2, DUMVAL, DUMVAL),
                   var_k = c(0.2, 0.2, DUMVAL, DUMVAL),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, T, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 100  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:2] <- prob_min(-lpost$post3_2_worst[, 1:2])
  set.seed(1)
  a_s$p_beat_soc[1:2] <- prob_lt(lpost$post3_2_worst[, 1:2])
  a_s$var_k[1:2] <- apply(lpost$post3_2_worst[, 1:2], 2, var)
  
  lpar$thresh_fut <- 0.2
  lpar$prob_symp = c(0.15, 0.2, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_2_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == c(T, F, F, F))
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == c(F, T, T, F))
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[2]==nk[2])
  stopifnot(res$arm_status$fut_at[3]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 4)])))

  stopifnot(res$trial_status$stop == T)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == T)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 100)
  
  
  message("\n\n TESTS ON 4 ARMS")
  
  message("----> 10 TEST Decision with 4 arms active. All arms are equivalent. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0, 0),
                   p_beat_soc = c(0, 0.2, 0, 0),
                   var_k = c(0.2, 0.2, 0, 0),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:4] <- prob_min(-lpost$post4_all_equal[, 1:4])
  set.seed(1)
  a_s$p_beat_soc[1:4] <- prob_lt(lpost$post4_all_equal[, 1:4])
  a_s$var_k[1:4] <- apply(lpost$post4_all_equal[, 1:4], 2, var)
  
  lpar$thresh_fut <- 0.1
  lpar$prob_symp = c(0.15, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_all_equal, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(all(is.na(res$arm_status$fut_at)))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  
  message("----> 11 TEST Decision with 4 arms active. All arms are equivalent. Arm 4 not in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0, DUMVAL),
                   p_beat_soc = c(0, 0.2, 0, DUMVAL),
                   var_k = c(0.2, 0.2, 0, DUMVAL),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  # NOTE It is correct to use post3_all_equal rather than post4_all_equal because
  # we only want to represent 3 arms in the posterior draws.
  set.seed(1)
  a_s$p_best[1:3] <- prob_min(-lpost$post3_all_equal[, 1:3])
  set.seed(1)
  a_s$p_beat_soc[1:3] <- prob_lt(lpost$post3_all_equal[, 1:3])
  a_s$var_k[1:3] <- apply(lpost$post3_all_equal[, 1:3], 2, var)
  
  lpar$thresh_fut <- 0.1
  lpar$prob_symp = c(0.15, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_all_equal, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(all(is.na(res$arm_status$fut_at)))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  

  
  message("----> 12 TEST Decision with 4 arms. All arms bar arm 3 active. All arms are equivalent. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, T),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, NA, 10, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, DUMVAL, 0),
                   p_beat_soc = c(0, 0.2, DUMVAL, 0),
                   var_k = c(0.2, 0.2, DUMVAL, 0),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, T, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 2, 4)] <- prob_min(-lpost$post4_all_equal[, c(1, 2, 4)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 2, 4)] <- prob_lt(lpost$post4_all_equal[, c(1, 2, 4)])
  a_s$var_k[c(1, 2, 4)] <- apply(lpost$post4_all_equal[, c(1, 2, 4)], 2, var)
  
  lpar$thresh_fut <- 0.1
  lpar$prob_symp = c(0.15, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_all_equal, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[3]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 2, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)

  
  message("----> 13 TEST Decision with 4 arms. All arms bar arm 3 active. All arms are equivalent. Arm 4 not in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, T),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, NA, 10, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, DUMVAL, DUMVAL),
                   p_beat_soc = c(0, 0.2, DUMVAL, DUMVAL),
                   var_k = c(0.2, 0.2, DUMVAL, DUMVAL),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, T, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  # NOTE It is correct to use post3_all_equal rather than post4_all_equal because
  # we only want to represent 3 arms in the posterior draws.
  set.seed(1)
  a_s$p_best[1:2] <- prob_min(-lpost$post3_all_equal[, 1:2])
  set.seed(1)
  a_s$p_beat_soc[1:2] <- prob_lt(lpost$post3_all_equal[, 1:2])
  a_s$var_k[1:2] <- apply(lpost$post3_all_equal[, 1:2], 2, var)
  
  lpar$thresh_fut <- 0.1
  lpar$prob_symp = c(0.15, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_all_equal, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[3]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 2, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  
  message("----> 14 TEST Decision with 4 arms active. SOC best. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0, 0),
                   p_beat_soc = c(0, 0.2, 0, 0),
                   var_k = c(0.2, 0.2, 0, 0),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:4] <- prob_min(-lpost$post4_1_best[, 1:4])
  set.seed(1)
  a_s$p_beat_soc[1:4] <- prob_lt(lpost$post4_1_best[, 1:4])
  a_s$var_k[1:4] <- apply(lpost$post4_1_best[, 1:4], 2, var)
  
  lpar$thresh_fut <- 0.01
  lpar$prob_symp = c(0.1, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_1_best, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(all(is.na(res$arm_status$fut_at)))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  
  
  message("----> 15 TEST Decision with 4 arms active. SOC best. All in post. All other arms futile.")
  nk <- c(36,30,34, 40); # sum(nk)
  # Identical to the above test but with an weaker futility threshold.
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0, 0),
                   p_beat_soc = c(0, 0.2, 0, 0),
                   var_k = c(0.2, 0.2, 0, 0),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim 
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:4] <- prob_min(-lpost$post4_1_best[, 1:4])
  set.seed(1)
  a_s$p_beat_soc[1:4] <- prob_lt(lpost$post4_1_best[, 1:4])
  a_s$var_k[1:4] <- apply(lpost$post4_1_best[, 1:4], 2, var)
  
  lpar$thresh_fut <- 0.2
  lpar$prob_symp = c(0.1, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_1_best, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == c(T, F, F, F))
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == c(F, T, T, T))
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[c(2, 3, 4)] == nk[2:4])
  stopifnot(is.na(res$arm_status$fut_at[1]))
  
  stopifnot(res$trial_status$stop == T)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == T)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)

  
  message("----> 16 TEST Decision with 4 arms active. SOC best. Arm 4 not in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0, DUMVAL),
                   p_beat_soc = c(0, 0.2, 0, DUMVAL),
                   var_k = c(0.2, 0.2, 0, DUMVAL),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  # Correct to be using 3 col post again.
  set.seed(1)
  a_s$p_best[1:3] <- prob_min(-lpost$post3_1_best[, 1:3])
  set.seed(1)
  a_s$p_beat_soc[1:3] <- prob_lt(lpost$post3_1_best[, 1:3])
  a_s$var_k[1:3] <- apply(lpost$post3_1_best[, 1:3], 2, var)
  
  lpar$thresh_fut <- 0.02
  lpar$prob_symp = c(0.1, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_1_best, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(all(is.na(res$arm_status$fut_at)))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  message("----> 17 TEST Decision with 4 arms. All arms bar arm 3 active. SOC best. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, T),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, NA, 10, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, DUMVAL, 0),
                   p_beat_soc = c(0, 0.2, DUMVAL, 0),
                   var_k = c(0.2, 0.2, DUMVAL, 0),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, T, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 2, 4)] <- prob_min(-lpost$post4_1_best[, c(1, 2, 4)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 2, 4)] <- prob_lt(lpost$post4_1_best[, c(1, 2, 4)])
  a_s$var_k[c(1, 2, 4)] <- apply(lpost$post4_1_best[, c(1, 2, 4)], 2, var)
  
  lpar$thresh_fut <- 0.02
  lpar$prob_symp = c(0.1, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_1_best, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[3]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 2, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  
  message("----> 18 TEST Decision with 4 arms. All arms bar arm 3 active. SOC best. Arm 4 not in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, T),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, NA, 10, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, DUMVAL, DUMVAL),
                   p_beat_soc = c(0, 0.2, DUMVAL, DUMVAL),
                   var_k = c(0.2, 0.2, DUMVAL, DUMVAL),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, T, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 2)] <- prob_min(-lpost$post3_1_best[, c(1, 2)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 2)] <- prob_lt(lpost$post3_1_best[, c(1, 2)])
  a_s$var_k[c(1, 2)] <- apply(lpost$post3_1_best[, c(1, 2)], 2, var)
  
  lpar$thresh_sup <- 0.9
  lpar$thresh_fut <- 0.02
  lpar$prob_symp = c(0.1, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_1_best, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[3]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 2, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  message("----> 19 TEST Decision with 4 arms active. Arm 2 best. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0, 0),
                   p_beat_soc = c(0, 0.2, 0, 0),
                   var_k = c(0.2, 0.2, 0, 0),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:4] <- prob_min(-lpost$post4_2_best[, 1:4])
  set.seed(1)
  a_s$p_beat_soc[1:4] <- prob_lt(lpost$post4_2_best[, 1:4])
  a_s$var_k[1:4] <- apply(lpost$post4_2_best[, 1:4], 2, var)
  
  lpar$thresh_sup <- 0.9
  lpar$thresh_fut <- 0.02
  lpar$prob_symp = c(0.15, 0.1, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_2_best, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)

  stopifnot(res$arm_status$is_sup == c(F, T, F, F))
  stopifnot(res$arm_status$is_fut == a_s$is_fut)

  stopifnot(res$arm_status$eff_at[2]==nk[2])
  stopifnot(all(is.na(res$arm_status$eff_at[c(1, 3, 4)])))
  stopifnot(all(is.na(res$arm_status$fut_at)))
  
  stopifnot(res$trial_status$stop == T)
  stopifnot(res$trial_status$win == T)
  stopifnot(res$trial_status$eff == T)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == T)
  stopifnot(res$trial_status$nk == 140)

  message("----> 20 TEST Decision with 4 arms active. Arm 3 best. Arm 4 not in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, 0, DUMVAL),
                   p_beat_soc = c(0, 0.2, 0, DUMVAL),
                   var_k = c(0.2, 0.2, 0, DUMVAL),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:3] <- prob_min(-lpost$post3_3_best[, 1:3])
  set.seed(1)
  a_s$p_beat_soc[1:3] <- prob_lt(lpost$post3_3_best[, 1:3])
  a_s$var_k[1:3] <- apply(lpost$post3_3_best[, 1:3], 2, var)
  
  lpar$thresh_sup <- 0.8
  lpar$thresh_fut <- 0.02
  lpar$prob_symp = c(0.15, 0.15, 0.1, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_3_best, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)

  stopifnot(res$arm_status$is_sup == c(F, F, T, F))
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(res$arm_status$eff_at[3]==nk[3])
  stopifnot(all(is.na(res$arm_status$eff_at[c(1, 2, 4)])))
  stopifnot(all(is.na(res$arm_status$fut_at)))
  
  stopifnot(res$trial_status$stop == T)
  stopifnot(res$trial_status$win == T)
  stopifnot(res$trial_status$eff == T)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == T)
  stopifnot(res$trial_status$nk == 140)
  
  
  message("----> 21 TEST Decision with 4 arms. All arms bar arm 3 active. Arm 2 best. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, T),
                   enabled_for_anly = c(1, 1, 1, 2),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, NA, 10, NA),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0.2, DUMVAL, 0),
                   p_beat_soc = c(0, 0.2, DUMVAL, 0),
                   var_k = c(0.2, 0.2, DUMVAL, 0),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, T, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 2, 4)] <- prob_min(-lpost$post4_2_best[, c(1, 2, 4)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 2, 4)] <- prob_lt(lpost$post4_2_best[, c(1, 2, 4)])
  a_s$var_k[c(1, 2, 4)] <- apply(lpost$post4_2_best[, c(1, 2, 4)], 2, var)
  
  lpar$thresh_sup <- 0.9
  lpar$thresh_fut <- 0.02
  lpar$prob_symp = c(0.15, 0.1, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_2_best, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == c(F, T, F, F))
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(res$arm_status$eff_at[2]==nk[2])
  stopifnot(all(is.na(res$arm_status$eff_at[c(1, 3, 4)])))
  stopifnot(res$arm_status$fut_at[3]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 2, 4)])))
  
  stopifnot(res$trial_status$stop == T)
  stopifnot(res$trial_status$win == T)
  stopifnot(res$trial_status$eff == T)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == T)
  stopifnot(res$trial_status$nk == 140)
  
  
  
  
  
  message("----> 22 TEST Decision with 4 arms. All arms bar arm 2 active. Arm 3 best. Arm 4 not in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, F, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, 10, NA, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1,DUMVAL, 0, DUMVAL),
                   p_beat_soc = c(0, DUMVAL, 0, DUMVAL),
                   var_k = c(0.2, DUMVAL, 0, DUMVAL),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, T, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 3)] <- prob_min(-lpost$post3_3_best[, c(1, 3)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 3)] <- prob_lt(lpost$post3_3_best[, c(1, 3)])
  a_s$var_k[c(1, 3)] <- apply(lpost$post3_3_best[, c(1, 3)], 2, var)
  
  lpar$thresh_sup <- 0.8
  lpar$thresh_fut <- 0.02
  lpar$prob_symp = c(0.15, 0.15, 0.1, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_3_best, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == c(F, F, T, F))
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(res$arm_status$eff_at[3]==nk[3])
  stopifnot(all(is.na(res$arm_status$eff_at[c(1, 2, 4)])))
  stopifnot(res$arm_status$fut_at[2]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 3, 4)])))
  
  stopifnot(res$trial_status$stop == T)
  stopifnot(res$trial_status$win == T)
  stopifnot(res$trial_status$eff == T)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == T)
  stopifnot(res$trial_status$nk == 140)

  message("----> 23 TEST Decision with 4 arms active. SOC worst. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1,0, 0, 0),
                   p_beat_soc = c(0, 0, 0, 0),
                   var_k = c(0.2, 0, 0, 0),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:4] <- prob_min(-lpost$post4_1_worst[, 1:4])
  set.seed(1)
  a_s$p_beat_soc[1:4] <- prob_lt(lpost$post4_1_worst[, 1:4])
  a_s$var_k[1:4] <- apply(lpost$post4_1_worst[, 1:4], 2, var)
  
  lpar$thresh_sup <- 0.95
  lpar$thresh_fut <- 0.02
  lpar$prob_symp = c(0.2, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_1_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  # Arm 2 should have been deactivated
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)

  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(all(is.na(res$arm_status$fut_at)))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  
  message("----> 24 TEST Decision with 4 arms active. SOC worst. Arm 4 not in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1,0, 0, DUMVAL),
                   p_beat_soc = c(0, 0, 0, DUMVAL),
                   var_k = c(0.2, 0, 0, DUMVAL),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:3] <- prob_min(-lpost$post3_1_worst[, 1:3])
  set.seed(1)
  a_s$p_beat_soc[1:3] <- prob_lt(lpost$post3_1_worst[, 1:3])
  a_s$var_k[1:3] <- apply(lpost$post3_1_worst[, 1:3], 2, var)
  
  lpar$thresh_sup <- 0.95
  lpar$thresh_fut <- 0.02
  lpar$prob_symp = c(0.2, 0.15, 0.15, 0.15)

  set.seed(1)
  res <- decision_intrm(post = lpost$post3_1_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(all(is.na(res$arm_status$fut_at)))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  message("----> 25 TEST Decision with 4 arms. All arms bar arm 3 active. SOC worst. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, T),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, NA, 10, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1,0, DUMVAL, 0),
                   p_beat_soc = c(0, 0, DUMVAL, 0),
                   var_k = c(0.2, 0, DUMVAL, 0 ),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, T, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 2, 4)] <- prob_min(-lpost$post4_1_worst[, c(1, 2, 4)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 2, 4)] <- prob_lt(lpost$post4_1_worst[, c(1, 2, 4)])
  a_s$var_k[c(1, 2, 4)] <- apply(lpost$post4_1_worst[, c(1, 2, 4)], 2, var)
  
  lpar$thresh_sup <- 0.95
  lpar$thresh_fut <- 0.02
  lpar$prob_symp = c(0.2, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_1_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[3]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 2, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  
  
  message("----> 26 TEST Decision with 4 arms. All arms bar arm 3 active. SOC worst. Arm 4 not in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, T),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, NA, 10, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1,0, DUMVAL, DUMVAL),
                   p_beat_soc = c(0, 0, DUMVAL, DUMVAL),
                   var_k = c(0.2, 0, DUMVAL, DUMVAL ),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, T, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 2)] <- prob_min(-lpost$post3_1_worst[, c(1, 2)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 2)] <- prob_lt(lpost$post3_1_worst[, c(1, 2)])
  a_s$var_k[c(1, 2)] <- apply(lpost$post3_1_worst[, c(1, 2)], 2, var)
  
  lpar$thresh_sup <- 0.95
  lpar$thresh_fut <- 0.02
  lpar$prob_symp = c(0.2, 0.15, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_1_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == a_s$active)
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == a_s$is_fut)
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[3]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 2, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == F)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  
  message("----> 27 TEST Decision with 4 arms active. Arm 3 worst. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1,0, 0, 0),
                   p_beat_soc = c(0, 0, 0, 0),
                   var_k = c(0.2, 0, 0, 0 ),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:4] <- prob_min(-lpost$post4_3_worst[, 1:4])
  set.seed(1)
  a_s$p_beat_soc[1:4] <- prob_lt(lpost$post4_3_worst[, 1:4])
  a_s$var_k[1:4] <- apply(lpost$post4_3_worst[, 1:4], 2, var)
  
  lpar$thresh_sup <- 0.95
  lpar$thresh_fut <- 0.2
  lpar$prob_symp = c(0.15, 0.15, 0.2, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_3_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == c(T, T, F, T))
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == c(F, F, T, F))
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[3]==nk[3])
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 2, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  

  message("----> 28 TEST Decision with 4 arms active. Arm 2 worst. Arm 4 not in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1,0, 0, DUMVAL),
                   p_beat_soc = c(0, 0, 0, DUMVAL),
                   var_k = c(0.2, 0, 0, DUMVAL ),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:3] <- prob_min(-lpost$post3_2_worst[, 1:3])
  set.seed(1)
  a_s$p_beat_soc[1:3] <- prob_lt(lpost$post3_2_worst[, 1:3])
  a_s$var_k[1:3] <- apply(lpost$post3_2_worst[, 1:3], 2, var)
  
  lpar$thresh_sup <- 0.95
  lpar$thresh_fut <- 0.2
  lpar$prob_symp = c(0.15, 0.2, 0.15, 0.15)
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_2_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == c(T, F, T, T))
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$arm_status$is_sup == a_s$is_sup)
  stopifnot(res$arm_status$is_fut == c(F, T, F, F))
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[2]==nk[2])
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 3, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)

  
  message("----> 29 TEST Decision with 4 arms. All arms bar arm 2 active. Arm 3 worst. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, F, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, 10, NA, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, DUMVAL, 0, DUMVAL),
                   p_beat_soc = c(0, DUMVAL, 0, DUMVAL),
                   var_k = c(0.2, DUMVAL, 0, DUMVAL ),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, T, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 3, 4)] <- prob_min(-lpost$post4_3_worst[, c(1, 3, 4)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 3, 4)] <- prob_lt(lpost$post4_3_worst[, c(1, 3, 4)])
  a_s$var_k[c(1, 3, 4)] <- apply(lpost$post4_3_worst[, c(1, 3, 4)], 2, var)
  
  lpar$thresh_sup <- 0.95
  lpar$thresh_fut <- 0.2
  lpar$prob_symp = c(0.15, 0.2, 0.15, 0.15)
  
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_3_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == c(T, F, F, T))
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$vret$is_sup == a_s$is_sup)
  stopifnot(res$vret$is_fut == c(F, T, T, F))
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[2]==10)
  stopifnot(res$arm_status$fut_at[3]==nk[3])
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  
  message("----> 30 TEST Decision with 4 arms. All arms bar arm 3 active. Arm 2 worst. Arm 4 not in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, T),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, NA, 10, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, F),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0, DUMVAL, DUMVAL),
                   p_beat_soc = c(0, 0, DUMVAL, DUMVAL),
                   var_k = c(0.2, 0, DUMVAL, DUMVAL ),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, T, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 2)] <- prob_min(-lpost$post3_2_really_worse[, c(1, 2)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 2)] <- prob_lt(lpost$post3_2_really_worse[, c(1, 2)])
  a_s$var_k[c(1, 2)] <- apply(lpost$post3_2_really_worse[, c(1, 2)], 2, var)
  
  lpar$thresh_sup <- 0.95
  lpar$thresh_fut <- 0.3
  lpar$prob_symp = c(0.15, 0.2, 0.15, 0.15)
  
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post3_2_really_worse, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == c(T, F, F, T))
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$vret$is_sup == a_s$is_sup)
  stopifnot(res$vret$is_fut == c(F, T, T, F))
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[2]==nk[2])
  stopifnot(res$arm_status$fut_at[3]==10)
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 4)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  
  
  
  message("----> 31 TEST Decision with 4 arms active. Arm 4 worst. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, T, T),
                   eff_at = rep(NA, 4),
                   fut_at = rep(NA, 4),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0, 0, 0),
                   p_beat_soc = c(0, 0, 0, 0),
                   var_k = c(0.2, 0, 0, 0 ),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, F, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[1:4] <- prob_min(-lpost$post4_4_worst[, 1:4])
  set.seed(1)
  a_s$p_beat_soc[1:4] <- prob_lt(lpost$post4_4_worst[, 1:4])
  a_s$var_k[1:4] <- apply(lpost$post4_4_worst[, 1:4], 2, var)
  
  lpar$thresh_sup <- 0.95
  lpar$thresh_fut <- 0.2
  lpar$prob_symp = c(0.15, 0.15, 0.15, 0.2)
  
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_4_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  # Arm 2 should have been deactivated
  stopifnot(res$arm_status$active == c(T, T, T, F))
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$vret$is_sup == a_s$is_sup)
  stopifnot(res$vret$is_fut == c(F, F, F, T))
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[4]==nk[4])
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 2, 3)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
  
  
  message("----> 32 TEST Decision with 4 arms. All arms bar arm 3 active. Arm 4 worst. All in post.")
  nk <- c(36,30,34, 40); # sum(nk)
  a_s <- list(K = 4, 
                   active = c(T, T, F, T),
                   eff_at = rep(NA, 4),
                   fut_at = c(NA, NA, 10, NA),
                   enabled_for_anly = c(1, 1, 1, 2),
                   arms_in_post = c(T, T, T, T),
                   p_rand = rep(0, 4),
                   p_best = c(0.1, 0, DUMVAL, 0),
                   p_beat_soc = c(0, 0, DUMVAL, 0),
                   var_k = c(0.2, 0, DUMVAL, 0 ),
                   is_sup = c(F, F, F, F),
                   is_fut = c(F, F, T, F),
              nk = nk)
  
  trial_status <- list(stop = F,# trial should stop (for eff or fut) determined at interim
                       win = F, # won regardless of whether at interim or final
                       eff = F, # stopped trial at interim since found effective treatment
                       fut_or_equ = F, # stopped trial at interim for futility - all arms futile
                       arm = F, # stopped an arm for fut (not necessarily stopping the trial)
                       pb = F,  # picked the best arm
                       nk = 140  # current total number of clusters
  )
  
  set.seed(1)
  a_s$p_best[c(1, 2, 4)] <- prob_min(-lpost$post4_4_worst[, c(1, 2, 4)])
  set.seed(1)
  a_s$p_beat_soc[c(1, 2, 4)] <- prob_lt(lpost$post4_4_worst[, c(1, 2, 4)])
  a_s$var_k[c(1, 2, 4)] <- apply(lpost$post4_4_worst[, c(1, 2, 4)], 2, var)
  
  lpar$thresh_sup <- 0.95
  lpar$thresh_fut <- 0.2
  lpar$prob_symp = c(0.15, 0.15, 0.15, 0.2)
  
  
  set.seed(1)
  res <- decision_intrm(post = lpost$post4_4_worst, 
                        a_s = a_s, 
                        t_s = trial_status,
                        lpar)
  
  stopifnot(res$arm_status$K == a_s$K)
  stopifnot(res$arm_status$active == c(T, T, F, F))
  stopifnot(res$arm_status$enabled_for_anly == a_s$enabled_for_anly)
  stopifnot(res$arm_status$arms_in_post == a_s$arms_in_post)
  stopifnot(res$arm_status$p_rand == a_s$p_rand) 
  
  stopifnot(abs(res$arm_status$p_best - a_s$p_best)<0.1)
  stopifnot(abs(res$arm_status$p_beat_soc - a_s$p_beat_soc)<0.1)
  stopifnot(res$arm_status$var_k == a_s$var_k)
  
  stopifnot(res$vret$is_sup == a_s$is_sup)
  stopifnot(res$vret$is_fut == c(F, F, T, T))
  
  stopifnot(all(is.na(res$arm_status$eff_at)))
  stopifnot(res$arm_status$fut_at[3]==10)
  stopifnot(res$arm_status$fut_at[4]==nk[4])
  stopifnot(all(is.na(res$arm_status$fut_at[c(1, 2)])))
  
  stopifnot(res$trial_status$stop == F)
  stopifnot(res$trial_status$win == F)
  stopifnot(res$trial_status$eff == F)
  stopifnot(res$trial_status$fut_or_equ == F)
  stopifnot(res$trial_status$arm == T)
  stopifnot(res$trial_status$pb == F)
  stopifnot(res$trial_status$nk == 140)
  
}



test_picked_best <- function(){
  
  
  # picked best when picked worse
  message("----> TEST Decision on best when no superior arms and all true probs equal.")
  a_s <- list(is_sup = c(F, F, F, F))
  lpar <- list(prob_symp = rep(0.15, 4))
  res <- picked_best_arm(a_s, lpar)
  stopifnot(res == F)

  message("----> TEST Decision on best when no superior arms and one true effective arm.")
  a_s <- list(is_sup = c(F, F, F, F))
  lpar <- list(prob_symp = c(0.15, 0.1, 0.15, 0.15))
  res <- picked_best_arm(a_s, lpar)
  stopifnot(res == F)
  
  message("----> TEST Decision on best when no superior arms and two true effective arm.")
  a_s <- list(is_sup = c(F, F, F, F))
  lpar <- list(prob_symp = c(0.15, 0.1, 0.08, 0.15))
  res <- picked_best_arm(a_s, lpar)
  stopifnot(res == F)
  
  message("----> TEST Decision on best when one superior arm and all true probs equal.")
  a_s <- list(is_sup = c(F, F, T, F))
  lpar <- list(prob_symp = c(0.15, 0.15, 0.15, 0.15))
  res <- picked_best_arm(a_s, lpar)
  stopifnot(res == F)
  
  message("----> TEST Decision on best when one superior (but not best) arm and one true effective arm.")
  a_s <- list(is_sup = c(F, F, T, F))
  lpar <- list(prob_symp = c(0.15, 0.1, 0.15, 0.15))
  res <- picked_best_arm(a_s, lpar)
  stopifnot(res == F)
  
  message("----> TEST Decision on best when one superior (but not best) arm and two true effective arm.")
  a_s <- list(is_sup = c(F, F, T, F))
  lpar <- list(prob_symp = c(0.15, 0.1, 0.15, 0.08))
  res <- picked_best_arm(a_s, lpar)
  stopifnot(res == F)
  
  message("----> TEST Decision on best when one superior (that is best) arm and one true effective arm.")
  a_s <- list(is_sup = c(F, T, F, F))
  lpar <- list(prob_symp = c(0.15, 0.1, 0.15, 0.15))
  res <- picked_best_arm(a_s, lpar)
  stopifnot(res == T)
  
  message("----> TEST Decision on best when one superior (that is best) arm and two true effective arm.")
  a_s <- list(is_sup = c(F, T, F, F))
  lpar <- list(prob_symp = c(0.15, 0.1, 0.1, 0.15))
  res <- picked_best_arm(a_s, lpar)
  stopifnot(res == T)  
  
  message("----> TEST Decision on best when two superior (one is best) arm and one true effective arm.")
  a_s <- list(is_sup = c(F, T, T, F))
  lpar <- list(prob_symp = c(0.15, 0.1, 0.15, 0.15))
  res <- picked_best_arm(a_s, lpar)
  stopifnot(res == T)  
  
  message("----> TEST Decision on best when two superior (one is best) arm and two true effective arm.")
  a_s <- list(is_sup = c(F, T, T, F))
  lpar <- list(prob_symp = c(0.15, 0.08, 0.15, 0.1))
  res <- picked_best_arm(a_s, lpar)
  stopifnot(res == T)  
  
  
}


test_prob_min <- function(){
  
  message("----> TESTING prob_min")
 
  lpost <- post_samples()
  
  message("----> TEST Minimum when first arm is best")
  res <- prob_min(lpost$post3_1_best[, -4])
  stopifnot(all(res[1] < res[2:3]))  
  
  message("----> TEST Minimum when last arm is best")
  res <- prob_min(lpost$post3_3_best[, -4])
  stopifnot(all(res[3] < res[1:2]))  
  
  message("----> TEST Minimum when arm 2 is best")
  res <- prob_min(lpost$post4_2_best[, -5])
  stopifnot(all(res[2] < res[c(1, 3, 4)]))  
  
  message("----> TEST Minimum for single arm")
  res <- prob_min(lpost$post4_2_best[, -c(2, 3, 4, 5)])
  stopifnot(res == 1)  
  

  
}


test_prob_lt <- function(){
  
  message("----> TESTING prob_min")
  
  lpost <- post_samples()
  
  message("----> TEST Prob LT SOC when first arm is best")
  res <- prob_lt(lpost$post3_1_best[, -4])
  stopifnot(res[1]==0)
  stopifnot(all(res[2:3]>res[1]))  
  
  message("----> TEST Prob LT SOC when last arm is best")
  res <- prob_lt(lpost$post3_3_best[, -4])
  stopifnot(res[1]==0)
  stopifnot(all(res[3]>res[1:2]))  
  
  message("----> TEST Prob LT SOC when arm 2 is best")
  res <- prob_lt(lpost$post4_2_best[, -5])
  stopifnot(res[1]==0)
  stopifnot(all(res[2]>res[c(1, 3, 4)]))  

  message("----> TEST Prob LT SOC for single arm")
  res <- prob_lt(lpost$post4_2_best[, -c(2, 3, 4, 5)])
  stopifnot(res == 0)  
  
  
  
}


