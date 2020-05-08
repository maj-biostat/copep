
source("data.R")
source("util.R")



#' GLMM implemented in jags 
#'
#' @param niter 
#'
#' @return
#' @export
#'
#' @examples
#' bequiet = F
#' dobs=d
glmm_jags <- function(dobs, lpar, bequiet = T){
  
  
  df <- list(n = nrow(dobs), 
             y = as.double(dobs$y),
             arm = as.double(dobs$arm),
             clust = as.double(dobs$clustid),
             n.arm = length(unique(dobs$arm)),
             n.clust = length(unique(dobs$clustid))
  )
  
  jags_inits <- function(){
    list(b.arm=rnorm(df$n.arm, 0, 1),
         tau.arm = rexp(1, 5),
         tau.clust = rexp(1, 5),
         .RNG.seed= runif(1, 1, 2000), 
         .RNG.name='base::Wichmann-Hill')    
  }
  
  jags_mod <- jags.model("BayesJAGSModel.txt",
                         data=df,
                         n.chains=1, 
                         n.adapt=1000, 
                         inits=jags_inits, 
                         quiet=bequiet )  #quiet=TRUE
  
  pbar <- ifelse(bequiet, "none", "text")
  # run a burn-in of 1000 iterations using update()
  update(jags_mod,
         n.iter=1000 , 
         progress.bar=pbar)
  
  # now run a further 5000 iterations monitoring the parameters of interest
  mcmc_samples <- coda.samples(jags_mod,
                               c("b.arm", "sigu"),
                               10000 , 
                               thin = 2,
                               progress.bar=pbar)
  
  # raftery.diag(mcmc_samples)
  # plot(mcmc_samples)
  
  post <- as.data.frame(mcmc_samples[[1]])
  
  names(post) <- gsub("[", "_", fixed = T, x = names(post))
  names(post) <- gsub("]", "", fixed = T, x = names(post))
  
  return(post)
}







glmm_stan <- function(dobs, lpar, bequiet = T){

  
  if (lpar$stanmodel == 9){
    
    message(get_hash(), " glmm_stan, model 9 ")
    
    fit <- logistic_mixed9(
      y = dobs$y, 
      K = length(unique(dobs$arm)), 
      arm = dobs$arm, 
      clustid = dobs$clustid, 
      prior_only = 0,
      pri_mu_int = 0,
      pri_sd_int = 100, # 1.81,
      pri_mu_trt = 0,
      pri_sd_trt = 1,
      pri_var_nu = 5,  # currently default at 3
      pri_var_sig = 2, # currently default at 1
      chains  = 1,
      thin    = 1,
      warmup = 1000, 
      iter    = 5000,
      refresh = 0,
      cores = 1,
      control = list(adapt_delta=0.8))
    
  }

  
  # names(as.data.frame(fit))
  
  res <- do.call(cbind, rstan::extract(fit, pars = c("bsoc", "btrt")))
  
  ncols <- ncol(res)
  post_mu <- cbind(res[, 1], sweep(res[,2:ncols], 1, res[,1], "+"))
  
  res <- do.call(cbind, rstan::extract(fit, pars = c("tau2")))
  
  post_mu <- cbind(post_mu, sqrt(res[, 1]))
  
  # b.arm_1   b.arm_2   b.arm_3      sigu
  colnames(post_mu) <- c(paste0("b.arm_", 1:ncols), "sigu")
  
  # MCMCvis::MCMCtrace(post_mu) 
  
  # post_trt <- res[, 2:ncols]
  # colnames(post_trt) <- paste0("b", 1:(ncols-1))

  # prob best
  # test <- function(post_mu){
  #   m <- matrix(0, ncol = ncol(post_mu), nrow = nrow(post_mu))
  #   coorx <- matrix(c(1:nrow(post_mu), max.col(post_mu)), ncol = 2)
  #   m[coorx] <- 1
  #   colMeans(m)
  # }
  
  # prob lt soc
  # test <- function(m){
  #   colMeans(m < 0)
  # }

  # microbenchmark::microbenchmark(test(-post_mu), prob_min(-post_mu), times=1000L)
  
  return(post_mu)
}







# posterior <- as.array(fit, pars = c("bsoc", "btrt", "sd_1"))
# lp_p <- log_posterior(fit)
# np_p <- nuts_params(fit)
# ratios_p <- neff_ratio(fit,  pars = c("b0", "b1", "lp__"))
# # mcmc_parcoord(posterior, np = np_p, pars = c("alpha", "b[1]"))
# mcmc_pairs(posterior, np = np_p, pars = c("b0", "b1"),
#            off_diag_args = list(size = 0.75))
# bayesplot::mcmc_trace(posterior, facet_args = list(nrow = dim(posterior)[3])) +
#   xlab("Post-warmup iteration")
# 
# bayesplot::mcmc_acf(posterior, lags = 10)
# raftery.diag(fit)
# plot(mcmc_samples)


#   sm <- stan_model(file = 'a.stan', save_dso = TRUE)
#   save('sm', file = 'sm.RData')
#   
#   Then, I submitted my R file to the SLURM and it works :)
# 
# The R file includes the following:
#   
#   load("sm.RData")
# 
# … simulated data…
# 
# fit <- sampling(sm, data=list(K, N, J, y, dir_alpha ), pars=c("pi", "mu", 
# "theta", "beta", "alpha", "prob"), warmup = 2000, iter = 5000, chains = 3)





