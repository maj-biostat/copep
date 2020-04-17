
# library(tidyverse)
# library(randomizr)
# library(GLMMadaptive)
# library(extraDistr)
# 
# ggplot2::theme_set(ggplot2::theme_bw())
# ggplot2::theme_update(text = element_text(size = 11))
# ggplot2::theme_update(legend.position = "top")
# ggplot2::theme_update(axis.text.x = element_text(size = 11))
# ggplot2::theme_update(axis.text.y = element_text(size = 11))
# 






#' PDF for symptoms score
#'
#' @param idx 
#'
#' @return
#' @export
#'
#' @examples
pdf_symptoms <- function(idx = NULL){
  
  # From
  # https://www.worldometers.info/coronavirus/coronavirus-symptoms/
  
  # Also refer to
  # Table 4. Possible scenarios for outcomes at day 15 from @who2020
  
  
  # 80.9% of infections are mild (with flu-like symptoms) and can recover at home.
  # 13.8% are severe, developing severe diseases including pneumonia and shortness of breath.
  # 4.7% as critical and can include: respiratory failure, septic shock, and multi-organ failure.
  # in about 2% of reported cases the virus is fatal.
  
  # 1. Not hospitalized, no limitations on activities
  # 2. Not hospitalized, limitation on activities;
  # 3. Hospitalized, not requiring supplemental oxygen;
  # 4. Hospitalized, requiring supplemental oxygen;
  # 5. Hospitalized, on non-invasive ventilation or high flow oxygen devices;
  # 6. Hospitalized, on invasive mechanical ventilation or ECMO;
  # 7. Death.
  
  y <- c(42,38,8,7,2,1,2)
  y <- y / sum(y)
  names(y) <- 1:7
  
  return(y)
}







#' Generates simulated data for design 1.
#' In this design we consider the proportion of those that were infected by the 
#' index individual. 
#' The number infected (which is our cluster of interest) is governed 
#' by a simple truncated poisson distribution. We randomise pbo to the cluster
#' but we also simulate the randomisation process of the index individual so that
#' we can adjust the analyses. The rationale here is that an index individual 
#' under therapy may be less infectious than an index individual under pbo.
#' For simplicity we assume that all those in close contact become infected.
#' 
#' For each individual we generate a binary indicator for which 1 represents
#' infected and showing symptoms at day 15 and 0 represents no symptoms at day
#' 15. The linear predictor includes a random effect for cluster group membership
#' with the thinking being that more severe index cases will result in more 
#' severe disease in those infected.
#' Nominally our assumptions are that 60% will still show symptoms at day 15 in 
#' the PBO group and we hope to reduce this to 40% in the drug group.
#'
#' @return
#' @export
#'
#' @examples
getdat1 <- function(lpar, datid = NULL){

  enrl_rate <- numeric(lpar$Nmax)
  enrl_rate <- lpar$enrl_lwr
  for(i in 2:lpar$Nmax){
    if(enrl_rate[i-1] > lpar$enrl_upr-1){
      enrl_rate[i] <- lpar$enrl_upr
    } else {
      enrl_rate[i] <- enrl_rate[i-1] * lpar$enrl_inc
    }
  }
  
  entry_time <- cumsum(c(0, rexp(n = (lpar$Nmax-1), rate=enrl_rate)))
  # entry_time
  
  
  # Cluster size - household/other
  clust_size_1 <- rtpois(lpar$Nmax,
                       lambda = lpar$mu_n_household,
                       a = lpar$min_clust_size_1,
                       b = lpar$max_clust_size_1)
  
  
  
  id <- cumsum(rep(1, sum(clust_size_1)))
  n <- length(id)
  
  clustid <- rep(seq_along(clust_size_1), clust_size_1)
  clustididx <- c(1, diff(clustid))
  
  entry_time <- entry_time[clustid] 
  entry_time[clustididx == 0] <- entry_time[clustididx == 0] 
  entry_time <- round(entry_time * 7, 2)


  
  
  # random effect
  clust_u0 <- rnorm(lpar$Nmax, 0, lpar$sig_u0)
  clust_u0 <- clust_u0[clustid]
  
  # Status at day 15 in each arm - i.e. we produce all counter factuals
  p <- do.call(cbind, lapply(1:length(lpar$prob_symp), function(i){
    # Probability of being symptomatic at day 15
    # Assume that this applies regardless of being infected or not
    p <- plogis(qlogis(lpar$prob_symp)[i] + clust_u0)
    p
  }))
  colnames(p) <- paste0("p", 1:ncol(p))
  
  y <- do.call(cbind, lapply(1:length(lpar$prob_symp), function(i){
    y_tmp <- rbinom(n, size = 1, prob = p[, i]) 
    y_tmp
  }))
  
  colnames(y) <- paste0("y", 1:ncol(y))
  
  d <- tibble(id = id, 
              clustid = clustid,
              arm = rep(NA, length(id)),
              entry_time = entry_time,
              # infc_trnsmt = infected,
              clust_u0 = clust_u0) %>%
    dplyr::bind_cols(as_tibble(p), 
                     as_tibble(y))
  
  d$y <- rep(NA, nrow(d))
  
  if(!is.null(datid)){
    attr(d, "datid") <- datid
  }
  
  # range(d$entry_time)
  
  d
}



