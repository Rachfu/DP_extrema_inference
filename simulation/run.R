source("fun_model.R")
source("utils.R")
source("fun_data.R")
library(corpcor)
library(MASS)

set.seed(888)

simulation_Gaussian = function(k,BS,e){
    path_name = "Gaussian"
    log_file <- paste("log/", path_name, ".txt", sep = "")

    sink(log_file, append = TRUE)
    n = 400
    rep = 1000
    B = 200
    level = 0.05
    r = c(-Inf,1/30,1/15,1/10,1/5,0.5)
    r_for_cv = c(1/30,1/15,1/10,1/5)

    theta_true = c(rep(0,k-1),BS)
    k_par = k/2
    cov_matrix = diag(length(theta_true))
    v=5 # fold
    seed_split = 666
    model = "Gaussian"
    cat("Parameter settings for", model, ": \n")
    cat("epsilon=",e,"\n")
    cat("theta_true=",theta_true,"\n")
    cat("k_par=",k_par,"\n")

    result <- simulation(n,rep,B,level,e,r,r_for_cv,theta_true,k_par,cov_matrix=diag(length(theta_true)),model,v,seed_split)
    print(result)
    sink()
}

simulation_OLS = function(k,BS,e){
    path_name = "OLS"
    log_file <- paste("log/", path_name, ".txt", sep = "")

    sink(log_file, append = TRUE)
    n = 400
    rep = 1000
    B = 200
    level = 0.05
    r = c(-Inf,1/30,1/15,1/10,1/5,0.5)
    r_for_cv = c(1/30,1/15,1/10,1/5)

    theta_true = c(rep(0,k-1),BS)
    k_par = k/2
    cov_matrix = diag(length(theta_true))
    v=5 # fold
    seed_split = 666
    model = "OLS"
    cat("Parameter settings for", model, ": \n")
    cat("epsilon=",e,"\n")
    cat("theta_true=",theta_true,"\n")
    cat("k_par=",k_par,"\n")

    result <- simulation(n,rep,B,level,e,r,r_for_cv,theta_true,k_par,cov_matrix=diag(length(theta_true)),model,v,seed_split)
    print(result)
    sink()
}

simulation_t_dis = function(){
    path_name = "OLS_t"
    log_file <- paste("log/", path_name, ".txt", sep = "")
    
    sink(log_file, append = TRUE)
    n = 400
    rep = 1000
    B = 200
    level = 0.05
    e = 1.5
    r = c(-Inf,1/30,1/15,1/10,1/5,0.5)
    r_for_cv = c(1/30,1/15,1/10,1/5)

    BS = 0
    k = 2
    theta_true = c(rep(0,k-1),BS)
    k_par = k/2
    cov_matrix = diag(length(theta_true))
    v=5 # fold
    model = "OLS_t"
    seed=666
    cat("Parameter settings for", model, ": \n")
    cat("epsilon=",e,"\n")
    cat("theta_true=",theta_true,"\n")
    cat("k_par=",k_par,"\n")

    result <- simulation(n,rep,B,level,e,r,r_for_cv,theta_true,k_par,cov_matrix=diag(length(theta_true)),model,v,seed_split)
    print(result)
    sink()
}


simulation_OLS_list_BS = function(BS){
    path_name = "OLS_list_BS"
    log_file <- paste("log/", path_name, ".txt", sep = "")
    
    sink(log_file, append = TRUE)
    n = 400
    rep = 1000
    B = 200
    level = 0.05
    e = 1.5
    r = c(-Inf,1/30,1/15,1/10,1/5,0.5)
    r_for_cv = c(1/30,1/15,1/10,1/5)

    k = 2
    theta_true = c(rep(0,k-1),BS)
    k_par = k/2
    cov_matrix = diag(length(theta_true))
    v=5 # fold
    model = "OLS"
    seed = 666

    cat("Parameter settings for", model, ": \n")
    cat("epsilon=",e,"\n")
    cat("theta_true=",theta_true,"\n")
    cat("k_par=",k_par,"\n")

    result <- simulation(n,rep,B,level,e,r,r_for_cv,theta_true,k_par,cov_matrix=diag(length(theta_true)),model,v,seed_split)
    print(result)
    sink()
}
############################### run code ###############################

## Experiment of Gaussian 
# k_list = c(2,8)
# BS_list = c(0,1)
# e_list = c(1.5,5)
# for (k in k_list){
#     for (BS in BS_list){
#         for (e in e_list){
#             simulation_Gaussian(k,BS,e)
#         }
#     }
# }

## Experiment of OLS 
# k_list = c(2,8)
# BS_list = c(0,1)
# e_list = c(1.5,5)
# for (k in k_list){
#     for (BS in BS_list){
#         for (e in e_list){
#             simulation_OLS(k,BS,e)
#         }
#     }
# }

## Experiment of t_distributed error in OLS 
# simulation_t_dis()

## Experiment of list of BSs in OLS 
# BS_list = c(0.05,0.1,0.3,0.5,0.7,0.9)
# for (BS in BS_list){
#     simulation_OLS_list_BS(BS)
# }