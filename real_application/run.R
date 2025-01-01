nCores = 8
myCluster <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(myCluster)

source("fun_model.R")
source("fun_data.R")
library(speff2trial)
library(foreach)
library(doParallel)
library(MASS)
library(dplyr)
library(rmutil)
library(epiDisplay)

set.seed(666)
real_app = function(group,group_method){ 
    path_name = paste(group, collapse = '_')
    log_file <- paste("log/", path_name, ".txt", sep = "")
    sink(log_file, append = TRUE)
    start_time = Sys.time()
    print("Start time is:")
    print(start_time)
    r = c(1/30,1/15,1/10,1/5,0.5)
    level = 0.2
    e = 4000   

    cat("correction term:", r, "\n")
    cat("confident level:", level, "\n")
    cat("privacy budget:", e, "\n")

    B = 200
    n_r = length(r)
    rep = 100
    BS_true = 0
    n_subgroup <- length(group)
    k_n <- 2^n_subgroup
    k_par = k_n

    data = OLS_data_real(group_method)
    D = data$D
    Y = data$Y

    tmp = fun_linear_regression(D,Y,B,level,e,r,BS_true,k_par)
    best_subgroup = tmp$best_subgroup
    l_naive = tmp$nonpri_naive
    k_fold = 3

    result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("dplyr","MASS","rmutil","epiDisplay","tmvtnorm","matrixcalc","corpcor")) %dopar% {
        unname(fun_linear_regression_parallel(D,Y,B,level,e,r,BS_true,k_par,k_fold))
    }
    parallel_result = rowMeans(result)
    final_result = c(parallel_result,l_naive,best_subgroup)
    sprintf("%.3f", final_result)
    final_result

    end_time <- Sys.time() 
    run_time <- end_time - start_time
    print("Time for the 100 epoch:")
    print(run_time)
    sink()
}

############################### run code ###############################
## Experiment with subgroup of race and age
# group = c('age','race')
# group_method = 'rage_age'
# real_app(group,group_method)

## Experiment with subgroup of symptom and offtrt
# group = c('symptom','offtrt')
# group_method = 'symptom_offtrt'
# real_app(group,group_method)