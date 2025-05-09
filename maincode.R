#### Data Analyze Code
#Author: Tugay Kacak

##---Load required library and functions##
library(mvtnorm) 
library(dplyr)
library(tidyr)
library(bifactor)
library(tictoc)
library(fungible)
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

i = c #for loop

sim_data <- readRDS(paste0(output_folder, "/sim_data_condition_", c, ".RDS")) # get data

get_phi_target <- function(n_generals, specifics_per_general,       #get phi_target matrices 
                           generals_rho, specifics_rho) {
  
  nfactors <- n_generals + n_generals*specifics_per_general
  Phi_Target <- matrix(0, nfactors, nfactors)
  Phi_Target[1:n_generals, 1:n_generals] <- 1
  diag(Phi_Target) <- 1
  Phi_Weight <- 1-Phi_Target
  
  result <- list(Phi_Target = Phi_Target, Phi_Weight = Phi_Weight)
  return(result)
  
}
Phi_Target <- get_phi_target(cond[i,1],cond[i,2])$Phi_Target
Phi_Weight <- get_phi_target(cond[i,1],cond[i,2])$Phi_Weight


# Load necessary packages
library(doParallel)
library(foreach)
library(EGAnet)
library(bifactor)

# Set up parallel backend
ncores <- parallel::detectCores()
n.cluster <- parallel::makeCluster(ncores - 2)
doParallel::registerDoParallel(cl = n.cluster)

# Define necessary variables (assuming these exist)
# cond, i, n_rep, sim_data, Phi_Target, Phi_Weight should already be defined

# Run the parallel computation
bifa <- foreach(j = 1:n_rep, .packages = c("EGAnet", "bifactor")) %dopar% {
  tryCatch(
    {
      # Calculate the polychoric correlation matrix
      cor_matrix <- EGAnet::polychoric.matrix(sim_data[[j]])
      
      # Perform bifactor analysis
      bifactor_result <- bifactor::bifactor(
        cor_matrix,
        estimator = "uls",
        n_generals = cond[i, 1],
        n_groups   = cond[i, 1] * cond[i, 2],
        method     = "GSLiD",
        projection = "poblq",
        oblq_factors = cond[i, 1] * cond[i, 2],
        maxit = 100,
        PhiTarget = Phi_Target,
        PhiWeight = Phi_Weight,
        random_starts = 10,
        rot_control = list(maxit = 1e5, rotation = "oblimin"),
        verbose = FALSE
      )
      # Return only the 'bifactor' part
      bifactor_result$bifactor
    },
    error = function(e) {
      # On error, return NA
      NA
    }
  )
}

bifa_ml <- foreach(j = 1:n_rep, .packages = c("EGAnet", "bifactor")) %dopar% {
  tryCatch(
    {
      cor_matrix <- EGAnet::polychoric.matrix(sim_data[[j]])
      result <- bifactor::bifactor(
        cor_matrix,
        estimator = "ml",
        n_generals = cond[i, 1],
        n_groups   = cond[i, 1] * cond[i, 2],
        method     = "GSLiD",
        projection = "poblq",
        oblq_factors = cond[i, 1] * cond[i, 2],
        maxit      = 100,
        PhiTarget  = Phi_Target,
        PhiWeight  = Phi_Weight,
        random_starts = 10,
        rot_control   = list(maxit = 1e5,rotation = "oblimin"),
        verbose       = FALSE
      )$bifactor  # Extract and return the bifactor element
    },
    error = function(e) {
      # If an error occurs, return NA
      NA
    }
  )
}

# Stop the cluster after processing
parallel::stopCluster(n.cluster)

# SAVE ANALYZES
# Combine bifa and bifa_ml into one list
bifa_results <- list(
  ULS = bifa,
  ML = bifa_ml
)
# Save both into one RDS file
saveRDS(bifa_results, file = paste0(analyzes_folder, "/analysis_condition_", i, ".RDS"))


## Calculate Metrics ##

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
true_lambda <- model$lambda

g_index <- cond[i,1]
s_index <- cond[i,1]*cond[i,2]

## Accuracy via Congruence Coefficient
acc <- t(sapply(1:n_rep, function(x) {
  if (is.list(bifa[[x]]) && !is.null(bifa[[x]]$lambda)) {
    loadings <- bifa[[x]]$lambda
    out <- fungible::faAlign(F1 = true_lambda, F2 = loadings)
    c(
      ACC_general = mean(out$CC[1:g_index]),
      ACC_specific = mean(out$CC[(g_index + 1):s_index])
    )
  } else {
    # Return NA values if the replication didn't converge or has a different structure
    c(ACC_general = NA, ACC_specific = NA)
  }
}))
# Calculate the mean RMSE values ignoring any NAs from non-convergent replications
acc <- as.data.frame(acc) %>% colMeans(na.rm = TRUE)

acc_ml <- t(sapply(1:n_rep, function(x) {
  if (is.list(bifa_ml[[x]]) && !is.null(bifa_ml[[x]]$lambda)) {
    loadings <- bifa_ml[[x]]$lambda
    out <- fungible::faAlign(F1 = true_lambda, F2 = loadings)
    c(
      ACC_general = mean(out$CC[1:g_index]),
      ACC_specific = mean(out$CC[(g_index + 1):s_index])
    )
  } else {
    # Return NA values if the replication didn't converge or has a different structure
    c(ACC_general = NA, ACC_specific = NA)
  }
}))
# Calculate the mean RMSE values ignoring any NAs from non-convergent replications
acc_ml <- as.data.frame(acc_ml) %>% colMeans(na.rm = TRUE)


## RMSE for Loadings
rmse <- t(vapply(1:n_rep, function(x) {
  # Check if the replication result is a list and contains the "lambda" element
  if (is.list(bifa[[x]]) && !is.null(bifa[[x]]$lambda)) {
    loadings <- bifa[[x]]$lambda
    out <- fungible::faAlign(F1 = true_lambda, F2 = loadings)
    loadings <- out$F2
    c(
      RMSE_general = sqrt(mean((out$F2[, 1:g_index] - true_lambda[, 1:g_index])^2)),
      RMSE_specific = sqrt(mean((out$F2[, (g_index+1):s_index] - true_lambda[, (g_index+1):s_index])^2))
    )
  } else {
    # Return NA values for RMSE if the replication did not converge or the structure is different
    c(RMSE_general = NA, RMSE_specific = NA)
  }
}, numeric(2)))
# Calculate the mean RMSE values ignoring any NAs from non-convergent replications
rmse <- as.data.frame(rmse) %>% colMeans(na.rm = TRUE)

## RMSE for Loadings
rmse_ml <- t(vapply(1:n_rep, function(x) {
  # Check if the replication result is a list and contains the "lambda" element
  if (is.list(bifa_ml[[x]]) && !is.null(bifa_ml[[x]]$lambda)) {
    loadings <- bifa_ml[[x]]$lambda
    out <- fungible::faAlign(F1 = true_lambda, F2 = loadings)
    loadings <- out$F2
    c(
      RMSE_general = sqrt(mean((out$F2[, 1:g_index] - true_lambda[, 1:g_index])^2)),
      RMSE_specific = sqrt(mean((out$F2[, (g_index+1):s_index] - true_lambda[, (g_index+1):s_index])^2))
    )
  } else {
    # Return NA values for RMSE if the replication did not converge or the structure is different
    c(RMSE_general = NA, RMSE_specific = NA)
  }
}, numeric(2)))
# Calculate the mean RMSE values ignoring any NAs from non-convergent replications
rmse_ml <- as.data.frame(rmse_ml) %>% colMeans(na.rm = TRUE)

## Convergence for GSLiD estimations
heywood_bifa <- lapply(1:n_rep, function(x) {
  if (is.list(bifa[[x]]) && !is.null(bifa[[x]]$convergence)) {
    bifa[[x]]$convergence
  } else {
    NA
  }
}) %>% unlist()
# Calculate the average convergence flag across replications, ignoring NAs
heywood_bifa <- sum(heywood_bifa, na.rm = TRUE) / n_rep

## Convergence for Target rotation
heywood_target <- lapply(1:n_rep, function(x) {
  if (is.list(bifa[[x]]) && !is.null(bifa[[x]]$Target_convergence)) {
    bifa[[x]]$Target_convergence
  } else {
    NA  # Return NA if the element doesn't exist
  }
}) %>% unlist()
# Calculate the average convergence flag across replications, ignoring NAs
heywood_target <- sum(heywood_target, na.rm = TRUE) / n_rep

## Convergence for GSLiD estimations
heywood_bifa_ml <- lapply(1:n_rep, function(x) {
  if (is.list(bifa_ml[[x]]) && !is.null(bifa_ml[[x]]$convergence)) {
    bifa_ml[[x]]$convergence
  } else {
    NA
  }
}) %>% unlist()
# Calculate the average convergence flag across replications, ignoring NAs
heywood_bifa_ml <- sum(heywood_bifa_ml, na.rm = TRUE) / n_rep

## Convergence for Target rotation
heywood_target_ml <- lapply(1:n_rep, function(x) {
  if (is.list(bifa_ml[[x]]) && !is.null(bifa_ml[[x]]$Target_convergence)) {
    bifa_ml[[x]]$Target_convergence
  } else {
    NA  # Return NA if the element doesn't exist
  }
}) %>% unlist()
# Calculate the average convergence flag across replications, ignoring NAs
heywood_target_ml <- sum(heywood_target_ml, na.rm = TRUE) / n_rep

condition_results <- list(
  Method = c("ULS","ML"),
  ACC_general = c(acc["ACC_general"],acc_ml["ACC_general"]),
  ACC_specific = c(acc["ACC_specific"],acc_ml["ACC_specific"]),
  RMSE_general = c(rmse["RMSE_general"],rmse_ml["RMSE_general"]),
  RMSE_specific = c(rmse["RMSE_specific"],rmse_ml["RMSE_specific"]),
  Heywood_bifa = c(heywood_bifa, heywood_bifa_ml),
  Heywood_target = c(heywood_target, heywood_target_ml))

filename <- paste0(analyzes_folder, "/results/results_condition_", i, ".RDS")
saveRDS(condition_results, file = filename)

rm(list = ls())
toc()
# End of the script

