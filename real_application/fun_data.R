library(MASS)
library(speff2trial)
OLS_data_real = function(group_method){
  del_var <- c("pidnum", "arms", "cd40", "cd420", "cd496", "r", "cd80", "cens", "days")
  ACTG175 <- na.omit(ACTG175[,!colnames(ACTG175) %in% del_var])
  X <- ACTG175[, !colnames(ACTG175) %in% c("cd820","zprior","treat")]
  Y <- ACTG175[,"cd820"]
  Ti<- ACTG175$treat
  n = dim(ACTG175)[1]
  if (group_method == 'rage_age'){
    group = c('age','race')
    n_subgroup <- length(group)
    k_n <- 2^n_subgroup

    Z <- matrix(0,n,k_n)
    for (i in 1:n){ 
        indicator_2base = rep(0,n_subgroup)
        if (X$race[i] == 0) indicator_2base[1] = indicator_2base[1]+1
        if (X$age[i] > median(X$age)) indicator_2base[2] = indicator_2base[2]+1
        indicator_10base = indicator_2base[2]+indicator_2base[1]*2
        Z[i,indicator_10base+1] = 1
    }
    colnames(Z) <- c("G1","G2","G3","G4") 

    cat("X$age[i] > median(X$age)", median(X$age), "\n")
    cat("X$race[i] == 0 \n")


    DZ <- matrix(0,n,k_n)
    for (i in 1:n){
        DZ[i,] = Ti[i]*Z[i,]
    }
    colnames(DZ) <- c("TG1","TG2","TG3","TG4") 

    X = cbind(Z,X) 
    D = cbind(DZ,X) 

    D <- D[, !colnames(D) %in% c("race")]
  }

  if (group_method == 'symptom_offtrt'){
    group = c('symptom','offtrt')
    n_subgroup <- length(group)
    k_n <- 2^n_subgroup

    Z <- matrix(0,n,k_n)
    for (i in 1:n){ 
        indicator_2base = rep(0,n_subgroup)
        if (X$offtrt[i] == 0) indicator_2base[1] = indicator_2base[1]+1
        if (X$symptom[i] == 0) indicator_2base[2] = indicator_2base[2]+1
        indicator_10base = indicator_2base[2]+indicator_2base[1]*2
        Z[i,indicator_10base+1] = 1
    }
    colnames(Z) <- c("G1","G2","G3","G4") 

    cat("X$offtrt[i] == 0 \n")
    cat("X$symptom[i] == 0 \n")

    DZ <- matrix(0,n,k_n)
    for (i in 1:n){
        DZ[i,] = Ti[i]*Z[i,]
    }
    colnames(DZ) <- c("TG1","TG2","TG3","TG4") 

    X = cbind(Z,X) 
    D = cbind(DZ,X) 

    D <- D[, !colnames(D) %in% c("symptom")]
    D <- D[, !colnames(D) %in% c("offtrt")]
  }
  
  normalize_min_max <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
  }
  D$wtkg <- normalize_min_max(D$wtkg)
  D$age <- normalize_min_max(D$age)
  D$karnof <- normalize_min_max(D$karnof)
  D$preanti <- normalize_min_max(D$preanti)
  D = as.matrix(D)
  result = list(D,Y)
  names(result) = c("D","Y")
  return(result)
}