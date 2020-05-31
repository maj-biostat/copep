# !diagnostics off

library(tidyverse)
library(randomizr)
library(extraDistr)
library(parallel)
library(tictoc)
library(rjags)
library(coda)
library(digest)
library(rstanmod1)
library(rstanarm)

source("data.R")
source("util.R")
source("models.R")
source("randomise.R")
source("trials.R")
source("decision.R")

fitmodel <- function(cfg_interface = NULL){
 
  if(.Platform$OS.type == "unix") {
    ncores <- parallel::detectCores()-1
  } else {
    ncores <- 1
  }
  
  lpar <- list()
  lpar$Nmax <- 450
  # Just assume flat 10 per week
  lpar$enrl_lwr <- 3
  lpar$enrl_upr <- 10
  lpar$mu_n_household <- 1
  lpar$min_clust_size_1 <- 0
  lpar$max_clust_size_1 <- 5
  lpar$enrl_inc <- 1.01
  lpar$sig_u0 <- 0 # 0.3454296
  lpar$prob_symp <- c(0.2, 0.2, 0.2)
  lpar$arms_start_idx <- 8
  # rstan
  lpar$stanmodel <- 9 
  lpar$chains <- 1
  lpar$thin <- 1
  lpar$warmup <- 1000
  lpar$iter <- 5000
  lpar$refresh <- 0
  lpar$cores <- 1
  lpar$control <- list(adapt_delta=0.8)
  x=1
  
  
  ld <- mclapply(X=1:1000, mc.cores=ncores, FUN=function(x) {
    
    message(Sys.time(), " Starting x = ", x)
    
    # Response variable for all arms are generated regardless of whether
    # arm is active or not.
    dobs <- tryCatch(getdat1(lpar, x),
                  error=function(e) { 
                    message(" getdat1 " , e); 
                  })
    
    # my <- which.min(abs(d$entry_time - 100))
    # ggplot(d, aes(x = entry_time, y = id))+
    #   theme_bw()+
    #   geom_line()+
    #   geom_hline(yintercept = my)+
    #   geom_vline(xintercept = 100)+
    #   scale_x_continuous("Day")+
    #   scale_y_continuous("Cumulative patients enrolled", 
    #                      breaks = seq(0, max(d$id), by = 100))
    # csize <- as.numeric(table(d$clustid));    hist(csize)
                  
    arms <- seq_along(lpar$prob_symp)
    t <- as.numeric(table(dobs$clustid))
    
    rand <- c(sample(arms, size = length(arms), replace = F),
              sample(arms, size = lpar$Nmax - length(arms), replace = T))
    dobs$arm <- rep(rand, t)
    
    dobs$y[1:nrow(dobs)] <- unlist(lapply(1:nrow(dobs), function(j){
      y <- dobs[j, lpar$arms_start_idx + dobs$arm[j]]
      y
    }))
    
    lprior <- list()
    lprior$prior_only = 0
    lprior$prior_to_use = 2
    lprior$prior_soc_par1 = 7
    lprior$prior_soc_par2 = 0
    lprior$prior_soc_par3 = 2.5
    lprior$prior_trt_par1 = 10
    lprior$prior_trt_par2 = 0
    lprior$prior_trt_par3 = 1
    # student t nu
    lprior$pri_var_par1 = 5
    # student t scale
    lprior$pri_var_par2 = 2

    post <-  tryCatch(glmm_stan(dobs, lpar, lprior = lprior),
                      error=function(e) { 
                        message(" glmm_stan " , e); 
                      })
    
    # fit2 <- rstanarm::stan_glmer(y~as.factor(arm)+(1|clustid), 
    #                      data = dobs,
    #                      family = binomial(link = "logit"), 
    #                      prior = rstanarm::student_t(df = 7), 
    #                      prior_intercept = rstanarm::student_t(df = 7),
    #                      chains  = 1, 
    #                      thin = 1, 
    #                      warmup = 1000, 
    #                      iter = 5000,
    #                      refresh = 1000,
    #                      cores = 1, 
    #                      control = list(adapt_delta=0.8)
    #                      )
    # 
    # post2 <- as.matrix(fit2, pars = c("(Intercept)", "as.factor(arm)2", "as.factor(arm)3"))
    # ncols <- ncol(post2)
    # post2 <- cbind(post2[, 1], sweep(post2[,2:ncols], 1, post2[,1], "+"))
    
    # den1 <- density(post[, 3])
    # den2 <- density(post2[, 3])
    # plot(den1)
    # lines(den2, col = "red")
    # colMeans(post2[, 1:3])
    
    colMeans(post[, 1:3])
    
    

 
  })
  
  m <- do.call(rbind, ld)
  
  message(Sys.time(), " Done, saving output")
  
  saveRDS(m, "equivalence_quan_matrix.RDS")
  
  
}


assess <- function(){
  
  
  m <- readRDS("equivalence_quan_matrix.RDS")
  colMeans(m)
  
  m2 <- apply(m, 2, expit)
  colMeans(m2)
  
  den1 <- density(m2[,1])
  plot(den1)
    
}


fitmodel()