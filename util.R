# Don't use NA, it will feck up conditional logic
DUMVAL <- -9

MYHASH <- NULL

set_hash <- function(){
  MYHASH <<- digest(as.character(runif(1)), algo = "xxhash32")
}

get_hash <- function(){
  MYHASH
}



#' Probability a column is minimum
#'
#' @param mat Matrix of MC draws
#'
#' @return Numeric vector giving probability each column is the maximum
#' @export
prob_min <- function(m) {
  if(is.null(dim(m)) & is.vector(m)){
    message(get_hash(), " prob_min m passed as vector, converting to 1d matrix")
    dim(m) <- c(length(m), 1)
  }
  as.numeric(prop.table(table(factor(max.col(m), levels = 1:ncol(m)))))
}

#' Probability a column is minimum
#'
#' @param mat Matrix of MC draws
#'
#' @return Numeric vector giving probability each column is the maximum
#' @export
prob_lt <- function(soc, m) {
  if(is.null(dim(m)) & is.vector(m)){
    message(get_hash(), " prob_lt m passed as vector, converting to 1d matrix")
    return(c(0, mean(m < soc)))
  }
  unlist(lapply(1:(ncol(m)+1), function(x){
    if(x == 1){
      0
    } else{
      mean(m[, (x-1)] < soc)
    }
    }))
}

#' Probability column is equivalent to SoC
#'
#' @param mat Matrix of MC draws
#'
#' @return Numeric vector giving probability each column is equivalent to SoC
#' @export
prob_equiv <- function(soc, m, eq_delta) {
  if(is.null(dim(m)) & is.vector(m)){
    message(get_hash(), " prob_equiv m passed as vector, converting to 1d matrix")
    dim(m) <- c(length(m), 1)
  }
  unlist(lapply(1:(ncol(m)+1), function(x){

    if(x == 1){
      1
    } else{
      mean(abs(soc - m[, (x-1)]) < eq_delta)
    }
    
    }))
}

#' Calculate matrix of pairwise differences of MC draws
#'
#' @param mat The matrix of draws
#' @param trans A transformation function, default is identity \code{I()}.
#' @return A matrix of pairwise differences
#' @export
#' @importFrom utils combn
pairwise_diff <- function(mat, trans = I) {
  trans <- match.fun(trans)
  pair_comp <- combn(ncol(mat), 2)
  trans_mat <- trans(mat)
  pair_mat <- apply(pair_comp, 2, function(x) trans_mat[,x[1]] - trans_mat[,x[2]])
  colnames(pair_mat) <- apply(pair_comp, 2, paste, collapse = "-")
  return(pair_mat)
}


#' Make a pairwise comparison between parameter draws
#'
#' @param mat Matrix of parameter draws
#' @param eps Reference value, default is \code{0}.
#' @return Vector of probabilities that pairwise difference is greater than \code{eps}.
#' @export
pairwise_comp <- function(mat, eps = 0) {
  pmat <- pairwise_diff(mat)
  apply(pmat, 2, function(x) mean(x > eps))
}

#' Bayesian response adaptive randomisation
#'
#' @param pbest Probability arm is best
#' @param sampsize Current sample size allocated to arm
#' @param variance Current posterior variance
#' @param inactive Flag for if an arm should receive zero allocations
#' @param fix_ctrl Fix allocation to control by this amount
#'
#' @return A numeric vector giving the BRAR allocation probabilities
#' @export
brar <- function(pbest, sampsize, variance, inactive) {
  stopifnot(all(pbest >= 0))
  stopifnot(all(sampsize > 0))
  m <- length(pbest)
  r <- sqrt(pbest * variance / sampsize)
  r[inactive] <- 0
  w <- r / sum(r)
  
  return(w)
 
}








#' Parameter grid for scenarios
#'
#' @return
#' @export
#'
#' @examples
getdpars <- function(){
  
  # Figures for p1 p2 etc:
  # https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(20)30114-5/fulltext
  
  sig_u0_1 <- get_sd_for_icc(0.035)
  sig_u0_2 <- get_sd_for_icc(0.1)
  
  # Early stopping  trial
  # dpars <- as_tibble(expand.grid(Nc = c(250,300),
  #                                m_contacts = c(4),
  #                                sig_u0 = c(sig_u0_1),
  #                                p1 = c(0.15, 0.1, 0.075),
  #                                p2 = c(0.1, 0.075, 0.05),
  #                                p3 = c(0.12),
  #                                interim_at = c(20, 30, 50)))

  dpars <- tribble(~scen, ~Nc, ~m_contacts,  ~sig_u0, ~p1,    ~p2,   ~p3,  ~interim_at,  
                   "1-0",    300,     2,      sig_u0_1,  0.15, 0.15,  0.15,      30,     
                   "1-1",    300,     2,      sig_u0_1,  0.15, 0.1,   0.15,      30,     
                   "1-2a",   300,     2,      sig_u0_1,  0.15, 0.1,  0.125,      30,     
                   "1-2b",   300,     2,      sig_u0_1,  0.15, 0.1,    0.1,      30,     
                   "2-0",    300,     2,      sig_u0_2,  0.15, 0.15,  0.15,      30,     
                   "2-1",    300,     2,      sig_u0_2,  0.15, 0.1,   0.15,      30,     
                   "2-2a",   300,     2,      sig_u0_2,  0.15, 0.1,  0.125,      30,     
                   "2-2b",   300,     2,      sig_u0_2,  0.15, 0.1,    0.1,      30,     
                   "3-0",    300,     2,      sig_u0_1,  0.15, 0.15,  0.15,      20,     
                   "3-1",    300,     2,      sig_u0_1,  0.15, 0.1,   0.15,      20,     
                   "3-2a",   300,     2,      sig_u0_1,  0.15, 0.1,  0.125,      20,     
                   "3-2b",   300,     2,      sig_u0_1,  0.15, 0.1,    0.1,      20,     
                   "4-0",    300,     2,      sig_u0_2,  0.15, 0.15,  0.15,      20,     
                   "4-1",    300,     2,      sig_u0_2,  0.15, 0.1,   0.15,      20,     
                   "4-2a",   300,     2,      sig_u0_2,  0.15, 0.1,  0.125,      20,     
                   "4-2b",   300,     2,      sig_u0_2,  0.15, 0.1,    0.1,      20      
                   )

  dpars$id <- 1:nrow(dpars)
  
  dpars
}


checkprior1 <- function(a, b){
  
  # in jags parameterisation of gamma is
  # f(x|r, l) = (1/gamma(r)) * l^r * x^(r-1) * exp(-l x)
  
  # in r the parameterisation is in terms of shape l = 1/s
  # f(x|a, s) = (1/gamma(a)) * (1/s^a) * x^(a-1) * exp(- x/s)
  
  b <- 1/b
  
  tauarm <- rgamma(1000, a, b)
  sigarm <- sqrt(1/tauarm)
  barm <- rnorm(1000, 0, sigarm)
  
  hist(barm)
  
}



jags_setup <- function(){
  
  fileForJAGS <- "BayesJAGSModel.txt"
  cat("model {
  for (i in 1:n){
    y[i] ~ dbin (p.bound[i], 1)
    p.bound[i] <- max(0, min(1, p[i]))
    logit(p[i]) <- b.arm[arm[i]] + sigu*b.u[clust[i]]
  }
  
  for (j in 1:n.arm) {b.arm[j] ~ dnorm(-1.5, tau.arm)}
  for (j in 1:n.clust) {b.u[j] ~ dnorm(0, 1)}
  
  # Precision (inverse gamma on sig)
  tau.arm ~ dgamma(1,1)
  tau.clust ~ dgamma(1,1)
  
  sigu  <- pow(tau.clust, -0.5)
  
}", file= fileForJAGS)
  
}


jags_setup2 <- function(){
  
  fileForJAGS <- "BayesJAGSModel.txt"
  cat("model {
  for (i in 1:nsoc){
    y0[i] ~ dbin (p.bound[i], 1)
    p.bound[i] <- max(0, min(1, p[i]))
    logit(p[i]) <- b.soc + sigu*b.u[clust[i]]
  }
  for (i in 1:ntrt){
    y[i] ~ dbin (p.bound[i], 1)
    p.bound[i] <- max(0, min(1, p[i]))
    logit(p[i]) <- b.soc + b.trt[arm[i]] + sigu*b.u[clust[i]]
  }
  
  b.soc ~ dnorm(-1.7, tau.soc)
  for (j in 1:n.trt) {b.trt[j] ~ dnorm(0, tau.arm)}
  for (j in 1:n.clust) {b.u[j] ~ dnorm(0, 1)}
  
  # Precision (inverse gamma on sig)
  tau.soc ~ dgamma(1,1)
  tau.arm ~ dgamma(2,2)
  tau.clust ~ dgamma(1,1)
  
  sigu  <- pow(tau.clust, -0.5)
  
}", file= fileForJAGS)
  
}


#' Gets SD for a given icc assuming residual variance of (pi^2)/3
#' for the RE logistic model.
#'
#' @param icc 
#'
#' @return
#' @export
#'
#' @examples
get_sd_for_icc <- function(icc){
  var_re <- ((pi^2)/3)/(1/icc - 1)
  sqrt(var_re)
}

#' Returns the odds associated with probability \code{p}.
#'
#' @param p 
#'
#' @return
#' @export
#'
#' @examples
odd <- function(p){
  p[p==0] <- 0.001
  p[p==1] <- 0.999
  
  p/(1-p)
}


#' expit
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
expit <- function(x){
  1/(1+exp(-x))
}


#' Simple list of config items to be used in data simulation.
#'
#' @return
#' @export
#'
#' @examples
cfg <- function(trial_interface = NULL, outdir = NULL, nsim = 3){
  
  
  # Utility approach
  lpar <- list()
  
  lpar$save_all <- F
  
  lpar$stanmodel <- 9 # alternatively use 5
  
  lpar$outdir <- outdir
  lpar$trial_interface <- trial_interface
  lpar$nsim <- nsim
  lpar$scen_idx <- NA
  
  lpar$Nmax <- 200
  # Used to interpolate enrollment rate
  lpar$Nenrl_ramp <- 200
  
  # per week - but converted such that a subjects entry time 
  # is in days (not weeks) from start
  lpar$enrl_lwr <- 5
  lpar$enrl_upr <- 20
  lpar$enrl_inc <- 1.01
  
  # Interim every x clusters
  lpar$interim_start <- 150
  lpar$interim_at <- 40
  
  # Cluster size
  lpar$mu_n_household <- 1
  
  # implies minimum household of 1 other person
  # index cases with no other household members ineligible?
  lpar$min_clust_size_1 <- 0
  lpar$max_clust_size_1 <- 3
  
  # Also defines ordering of trt.
  lpar$trt_name <- c("pbo", "drg")
  
  # Prob of virus pos + sympt @ day 14 by trt arm (PBO, arm1, arm2, etc)
  # First slot always for PBO
  lpar$prob_symp <- c(0.15, 0.15, 0.15)
  
  # SD for cluster random effect
  lpar$sig_u0 <- 0.2
  
  lpar$test_at_day <- 10
  
  # Decision thresholds
  lpar$thresh_sup <- 0.975
  lpar$thresh_fut <- 0.15
  lpar$thresh_equ <- 0.92
  # this is on the log-odds scale
  lpar$eq_delta <- 0.45
  
  # Whether to keep all the interim posteriors or not
  lpar$keep_intrms <- F
  
  # RAR
  lpar$rand_par <- c(1, 1, 200, 5)
  
  lpar$arms_start_idx <- 8
  
  
  # Configuration to introduce new arms at a specified interim.
  # If no new arms are to be introduced set to rep(T, K) and rep(1, K)
  # where K is the total number of arms.
  lpar$arms_activ_at_start <- c(T, T, T)
  lpar$arms_enabled_for_anly <- c(1, 1, 1)

  return(lpar)
}




#' Simulation parameters and scenarios for fixed trial.
#' Method passed to sim function in main.
#' Contains the trials interface method along with the directory that 
#' output should be sent to. 
#'
#' @return
#' @export
#'
#' @examples
cfg_fixed <- function(trial_interface = NULL, outdir = NULL, nsim = 3){
  
  lpar <- cfg(trial_interface, outdir, nsim)

  # Make lpar contain the scenarios that you want to look at.  
  lpar$scenarios <- list()
  
  sig_u0_1 <- get_sd_for_icc(0.035)
  sig_u0_2 <- get_sd_for_icc(0.02)
  
  lpar$scenarios[[1]] <- list(scen = "1-0", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.15, 0.15),
                              arms_activ_at_start = c(T, T, T),
                              arms_enabled_for_anly = c(1, 1, 1),
                              thresh_sup = 0.95,
                              arms_start_idx=8,
                              interim_at = 400) # implies not to do interims
  
  lpar$scenarios[[2]] <- list(scen = "1-1", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.15),
                              arms_activ_at_start = c(T, T, T),
                              arms_enabled_for_anly = c(1, 1, 1),
                              thresh_sup = 0.95,
                              arms_start_idx=8,
                              interim_at = 400)
  
  lpar$scenarios[[3]] <- list(scen = "1-2", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.085, 0.15),
                              arms_activ_at_start = c(T, T, T),
                              arms_enabled_for_anly = c(1, 1, 1),
                              thresh_sup = 0.95,
                              arms_start_idx=8,
                              interim_at = 400)
  
  lpar$scenarios[[4]] <- list(scen = "1-3", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.125),
                              arms_activ_at_start = c(T, T, T),
                              arms_enabled_for_anly = c(1, 1, 1),
                              thresh_sup = 0.95,
                              arms_start_idx=8,
                              interim_at = 400)
  
  lpar$scenarios[[5]] <- list(scen = "1-4", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.1),
                              arms_activ_at_start = c(T, T, T),
                              arms_enabled_for_anly = c(1, 1, 1),
                              thresh_sup = 0.95,
                              arms_start_idx=8,
                              interim_at = 400)
  
  lpar$scenarios[[6]] <- list(scen = "2-0", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.15, 0.15),
                              arms_activ_at_start = c(T, T, T),
                              arms_enabled_for_anly = c(1, 1, 1),
                              thresh_sup = 0.95,
                              arms_start_idx=8,
                              interim_at = 400)
  
  lpar$scenarios[[7]] <- list(scen = "2-1", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.15),
                              arms_activ_at_start = c(T, T, T),
                              arms_enabled_for_anly = c(1, 1, 1),
                              thresh_sup = 0.95,
                              arms_start_idx=8,
                              interim_at = 400)
  
  lpar$scenarios[[8]] <- list(scen = "2-2", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.085, 0.15),
                              arms_activ_at_start = c(T, T, T),
                              arms_enabled_for_anly = c(1, 1, 1),
                              thresh_sup = 0.95,
                              arms_start_idx=8,
                              interim_at = 400)
  
  lpar$scenarios[[9]] <- list(scen = "2-3", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp =  c(0.15, 0.1, 0.125),
                              arms_activ_at_start = c(T, T, T),
                              arms_enabled_for_anly = c(1, 1, 1),
                              thresh_sup = 0.95,
                              arms_start_idx=8,
                              interim_at = 400)
  
  lpar$scenarios[[10]] <- list(scen = "2-4", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.1),
                              arms_activ_at_start = c(T, T, T),
                              arms_enabled_for_anly = c(1, 1, 1),
                              thresh_sup = 0.95,
                              arms_start_idx=8,
                              interim_at = 400)
  
  lpar
}


cfg_fixed4 <- function(trial_interface = NULL, outdir = NULL, nsim = 3){
  
  lpar <- cfg(trial_interface, outdir, nsim)
  
  # Make lpar contain the scenarios that you want to look at.  
  lpar$scenarios <- list()
  
  sig_u0_1 <- get_sd_for_icc(0.035)
  sig_u0_2 <- get_sd_for_icc(0.1)
  
  lpar$scenarios[[1]] <- list(scen = "1-0", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.15, 0.15, 0.15),
                              arms_activ_at_start = c(T, T, T, T),
                              arms_start_idx = 9,
                              arms_enabled_for_anly = c(1, 1, 1, 1),
                              interim_at = 400) # implies not to do interims

  
  lpar$scenarios[[2]] <- list(scen = "1-1", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.15, 0.15),
                              arms_activ_at_start = c(T, T, T, T),
                              arms_start_idx = 9,
                              arms_enabled_for_anly = c(1, 1, 1, 1),
                              interim_at = 400)
  
  lpar$scenarios[[3]] <- list(scen = "1-2", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.085, 0.15, 0.15),
                              arms_activ_at_start = c(T, T, T, T),
                              arms_start_idx = 9,
                              arms_enabled_for_anly = c(1, 1, 1, 1),
                              interim_at = 400)
  
  lpar$scenarios[[4]] <- list(scen = "1-3", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.125, 0.11),
                              arms_activ_at_start = c(T, T, T, T),
                              arms_start_idx = 9,
                              arms_enabled_for_anly = c(1, 1, 1, 1),
                              interim_at = 400)
  
  lpar$scenarios[[5]] <- list(scen = "1-4", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.1, 0.085),
                              arms_activ_at_start = c(T, T, T, T),
                              arms_start_idx = 9,
                              arms_enabled_for_anly = c(1, 1, 1, 1),
                              interim_at = 400)
  
  lpar$scenarios[[6]] <- list(scen = "2-0", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.15, 0.15, 0.15),
                              arms_activ_at_start = c(T, T, T, T),
                              arms_start_idx = 9,
                              arms_enabled_for_anly = c(1, 1, 1, 1),
                              interim_at = 400)
  
  lpar$scenarios[[7]] <- list(scen = "2-1", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.15, 0.15),
                              arms_activ_at_start = c(T, T, T, T),
                              arms_start_idx = 9,
                              arms_enabled_for_anly = c(1, 1, 1, 1),
                              interim_at = 400)
  
  lpar$scenarios[[8]] <- list(scen = "2-2", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.085, 0.15, 0.15),
                              arms_activ_at_start = c(T, T, T, T),
                              arms_start_idx = 9,
                              arms_enabled_for_anly = c(1, 1, 1, 1),
                              interim_at = 400)
  
  lpar$scenarios[[9]] <- list(scen = "2-3", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp =  c(0.15, 0.1, 0.125, 0.11),
                              arms_activ_at_start = c(T, T, T, T),
                              arms_start_idx = 9,
                              arms_enabled_for_anly = c(1, 1, 1, 1),
                              interim_at = 400)
  
  lpar$scenarios[[10]] <- list(scen = "2-4", 
                              Nc = 450, 
                              mu_n_household = 3, 
                              sig_u0_1 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.1, 0.085),
                              arms_activ_at_start = c(T, T, T, T),
                              arms_start_idx = 9,
                              arms_enabled_for_anly = c(1, 1, 1, 1),
                              interim_at = 400)
  
  lpar
}



cfg_rar_3 <- function(trial_interface = NULL, outdir = NULL, nsim = 3){
  
  lpar <- cfg(trial_interface, outdir, nsim)
  
  # Make lpar contain the scenarios that you want to look at.  
  lpar$scenarios <- list()
  
  sig_u0_1 <- get_sd_for_icc(0.035)
  sig_u0_2 <- get_sd_for_icc(0.1)
  
  lpar$scenarios[[1]] <- list(scen = "1-0", 
                              Nc = 450, 
                              mu_n_household = 1, 
                              sig_u0 = sig_u0_1, 
                              prob_symp = c(0.15, 0.15, 0.15),
                              arms_start_idx = 8,
                              interim_at = 20) # implies not to do interims
  
  lpar$scenarios[[2]] <- list(scen = "1-1", 
                              Nc = 450, 
                              mu_n_household = 1, 
                              sig_u0 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.15),
                              arms_start_idx = 8,
                              interim_at = 20)
  
  lpar$scenarios[[3]] <- list(scen = "1-2", 
                              Nc = 450, 
                              mu_n_household = 1, 
                              sig_u0 = sig_u0_1, 
                              prob_symp = c(0.15, 0.085, 0.15),
                              arms_start_idx = 8,
                              interim_at = 20)
  
  lpar$scenarios[[4]] <- list(scen = "1-3", 
                              Nc = 450, 
                              mu_n_household = 1, 
                              sig_u0 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.125),
                              arms_start_idx = 8,
                              interim_at = 20)
  
  lpar$scenarios[[5]] <- list(scen = "1-4", 
                              Nc = 450, 
                              mu_n_household = 1, 
                              sig_u0 = sig_u0_1, 
                              prob_symp = c(0.15, 0.1, 0.1),
                              arms_start_idx = 8,
                              interim_at = 20)
  
  # Different ICC
  
  lpar$scenarios[[6]] <- list(scen = "2-0", 
                              Nc = 450, 
                              mu_n_household = 1, 
                              sig_u0 = sig_u0_2, 
                              prob_symp = c(0.15, 0.15, 0.15),
                              arms_start_idx = 8,
                              interim_at = 20) # implies not to do interims
  
  lpar$scenarios[[7]] <- list(scen = "2-1", 
                              Nc = 450, 
                              mu_n_household = 1, 
                              sig_u0 = sig_u0_2, 
                              prob_symp = c(0.15, 0.1, 0.15),
                              arms_start_idx = 8,
                              interim_at = 20)
  
  lpar$scenarios[[8]] <- list(scen = "2-2", 
                              Nc = 450, 
                              mu_n_household = 1, 
                              sig_u0 = sig_u0_2, 
                              prob_symp = c(0.15, 0.085, 0.15),
                              arms_start_idx = 8,
                              interim_at = 20)
  
  lpar$scenarios[[9]] <- list(scen = "2-3", 
                              Nc = 450, 
                              mu_n_household = 1, 
                              sig_u0 = sig_u0_2, 
                              prob_symp = c(0.15, 0.1, 0.125),
                              arms_start_idx = 8,
                              interim_at = 20)
  
  lpar$scenarios[[10]] <- list(scen = "2-4", 
                              Nc = 450, 
                              mu_n_household = 1, 
                              sig_u0 = sig_u0_2, 
                              prob_symp = c(0.15, 0.1, 0.1),
                              arms_start_idx = 8,
                              interim_at = 20)
  
  
  
  
  lpar
  
}

cfg_rar_4 <- function(trial_interface = NULL, outdir = NULL, nsim = 3){
  
  lpar <- cfg_rar_3(trial_interface, outdir, nsim)
  
  # Scenarios 1-
  
  lpar$scenarios[[1]]$prob_symp <- c(0.15, 0.15, 0.15, 0.15)
  lpar$scenarios[[1]]$arms_activ_at_start <- c(T, T, T, F)
  lpar$scenarios[[1]]$arms_enabled_for_anly <- c(1, 1, 1, 2)
  lpar$scenarios[[1]]$arms_start_idx <- 9
  
  lpar$scenarios[[2]]$prob_symp <- c(0.15, 0.1, 0.15, 0.15)
  lpar$scenarios[[2]]$arms_activ_at_start <- c(T, T, T, F)
  lpar$scenarios[[2]]$arms_enabled_for_anly <- c(1, 1, 1, 2)
  lpar$scenarios[[2]]$arms_start_idx <- 9
  
  lpar$scenarios[[3]]$prob_symp <- c(0.15, 0.085, 0.15, 0.15)
  lpar$scenarios[[3]]$arms_activ_at_start <- c(T, T, T, F)
  lpar$scenarios[[3]]$arms_enabled_for_anly <- c(1, 1, 1, 2)
  lpar$scenarios[[3]]$arms_start_idx <- 9
  
  lpar$scenarios[[4]]$prob_symp <- c(0.15, 0.1, 0.125, 0.11)
  lpar$scenarios[[4]]$arms_activ_at_start <- c(T, T, T, F)
  lpar$scenarios[[4]]$arms_enabled_for_anly <- c(1, 1, 1, 2)
  lpar$scenarios[[4]]$arms_start_idx <- 9
  
  lpar$scenarios[[5]]$prob_symp <- c(0.15, 0.1, 0.1, 0.085)
  lpar$scenarios[[5]]$arms_activ_at_start <- c(T, T, T, F)
  lpar$scenarios[[5]]$arms_enabled_for_anly <- c(1, 1, 1, 2)
  lpar$scenarios[[5]]$arms_start_idx <- 9
  
  # Scenarios 2-
  
  lpar$scenarios[[6]]$prob_symp <- c(0.15, 0.15, 0.15, 0.15)
  lpar$scenarios[[6]]$arms_activ_at_start <- c(T, T, T, F)
  lpar$scenarios[[6]]$arms_enabled_for_anly <- c(1, 1, 1, 2)
  lpar$scenarios[[6]]$arms_start_idx <- 9
  
  lpar$scenarios[[7]]$prob_symp <- c(0.15, 0.1, 0.15, 0.15)
  lpar$scenarios[[7]]$arms_activ_at_start <- c(T, T, T, F)
  lpar$scenarios[[7]]$arms_enabled_for_anly <- c(1, 1, 1, 2)
  lpar$scenarios[[7]]$arms_start_idx <- 9
  
  lpar$scenarios[[8]]$prob_symp <- c(0.15, 0.085, 0.15, 0.15)
  lpar$scenarios[[8]]$arms_activ_at_start <- c(T, T, T, F)
  lpar$scenarios[[8]]$arms_enabled_for_anly <- c(1, 1, 1, 2)
  lpar$scenarios[[8]]$arms_start_idx <- 9
  
  lpar$scenarios[[9]]$prob_symp <- c(0.15, 0.1, 0.125, 0.11)
  lpar$scenarios[[9]]$arms_activ_at_start <- c(T, T, T, F)
  lpar$scenarios[[9]]$arms_enabled_for_anly <- c(1, 1, 1, 2)
  lpar$scenarios[[9]]$arms_start_idx <- 9
  
  lpar$scenarios[[10]]$prob_symp <- c(0.15, 0.1, 0.1, 0.085)
  lpar$scenarios[[10]]$arms_activ_at_start <- c(T, T, T, F)
  lpar$scenarios[[10]]$arms_enabled_for_anly <- c(1, 1, 1, 2)
  lpar$scenarios[[10]]$arms_start_idx <- 9
  
  
  lpar
  
}


