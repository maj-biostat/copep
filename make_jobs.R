library(tidyverse)

source("util.R")
source("trials.R")

# This is just temporary until I think of a better way to do it...
lpar <- cfg_rar_3(trial_interface = trial_GS_RAR, outdir = "outrar", nsim = 1)

n_machines <- 8

message("Schedulling ", length(lpar$scenarios), " scenarios across ", n_machines, " machines.")



# Load balance jobs across all available machines.
# Initialise job_indexes
job_indexes <- list()
for(i in 1:n_machines){
  job_indexes[[i]] <- numeric(0)
}

# While there are more scenarios, keep allocating to machines sequentially
idx <- 1
j <- 1
while(length(unlist(job_indexes)) < length(lpar$scenarios)){
  
  job_indexes[[idx]] <- c(job_indexes[[idx]], j)
  idx <- ifelse(idx == 8, 1, idx + 1)
  j <- j + 1
}

# For i in all machines available, write out its scenarios to cfg file.
for(i in 1:n_machines){
  
  message("Machine: ", i, " gets scenarios: ", 
          paste0(job_indexes[[i]], collapse = " "))
  cat(paste0(paste0(job_indexes[[i]], 
             collapse="\n"), "\n"), 
      file = paste0("scenarios", i, ".cfg"))
}







