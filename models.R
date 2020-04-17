
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
glmm_jags <- function(dobs, bequiet = T){
  

  df <- list(n = nrow(dobs), 
             y = dobs$y,
             arm = as.numeric(dobs$arm),
             clust = dobs$clustid,
             n.arm = length(unique(dobs$arm)),
             n.clust = length(unique(dobs$clustid))
  )
  
  jags_inits <- function(){
    list(b.arm=rnorm(df$n.arm, 0, 1),
         tau.arm = rexp(1, 5),
         tau.clust = rexp(1, 5),
         .RNG.seed=1, 
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
                               5000 , 
                               thin = 3,
                               progress.bar=pbar)
  
  # plot(mcmc_samples)
  
  post <- as.data.frame(mcmc_samples[[1]])
  
  names(post) <- gsub("[", "_", fixed = T, x = names(post))
  names(post) <- gsub("]", "", fixed = T, x = names(post))
  
  return(post)
}
















