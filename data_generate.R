#### Data Generation Code
#Author: Tugay Kacak

##---Load required library and functions##
library(mvtnorm) 
library(dplyr)
library(tidyr)
library(bifactor)
library(tictoc)
tic()
#fully-crossed design
cond <- expand.grid(
  NGEN = c(1,2,3),
  NFAC = c(3,4),
  NVAR = c(3,5),
  GLOAD = c("low","medium","high"),
  FLOAD = c("low","medium","high"),
  GRHO = c(0.00),
  FRHO = c(0.00,0.20,0.50),
  NOBS = c(250,500,750,1000),
  NCAT  = c(2,3,4,5,6,7))

current_condition <- cond

#Directory for saving datasets
output_folder <- "D:/Ranalysis/makale/small_bifac/datasets"
analyzes_folder <- "D:/Ranalysis/makale/small_bifac/analyzes"

n_conditions = nrow(current_condition)
n_rep = 100

#Multicore Analyzes
library(doParallel)
library(foreach)
#select cores
ncores = detectCores()
n.cluster = makeCluster(ncores - 4)
registerDoParallel(cl=n.cluster)

# Run the parallel simulation
foreach(i = 1:n_conditions, .packages = c("bifactor", "mvtnorm")) %dopar% {
  set.seed(1234 + i)
  
  # Generate the factor model
  model <- bifactor::sim_factor(
    n_generals = cond[i, "NGEN"],
    groups_per_general = cond[i, "NFAC"],
    items_per_group = cond[i, "NVAR"],
    loadings_g = cond[i, "GLOAD"],
    loadings_s = cond[i, "FLOAD"],
    crossloadings = 0.00,
    generals_rho = cond[i, "GRHO"],
    groups_rho = cond[i, "FRHO"],
    method = "minres"
  )
  
  # Generate datasets for each replication
  sim_data_i <- lapply(1:n_rep, function(j) {
    #dat <- MASS::mvrnorm(n = cond[i, "NOBS"], rep(0, nrow(model$R)), Sigma = model$R)
    dat <- data.frame(mvtnorm::rmvnorm(cond[i,"NOBS"], mean = rep(0,nrow(model$R)), sigma = model$R))
    dat <- data.frame(
             lapply(dat,
             FUN = function(x){cut(x, c(min(x)-0.01,
                                        qnorm((1:(cond[i,"NCAT"]-1))/cond[i,"NCAT"]),
                                        max(x) + 0.01), labels = FALSE,
                                        right =FALSE)}))
    return(dat)
  })
  
  # Construct the filename dynamically based on the condition number
  filename <- paste0(output_folder, "/sim_data_condition_", i, ".RDS")
  
  # Save the dataset for the current condition
  saveRDS(sim_data_i, file = filename)
  
  # Remove the dataset from memory and trigger garbage collection
  rm(sim_data_i)
}
stopCluster(n.cluster)

