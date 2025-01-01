library(MASS)
Gaussian_data = function(n,theta_true,cov_matrix){
  D = mvrnorm(n,theta_true, cov_matrix)
  return(D)
}

OLS_data = function(n,theta_true){
  k = length(theta_true)
  D = matrix(0,n,k)
  for (i in 1:n){
    index = sample(1:k,1)
    D[i,index] = 1
  }
  Y = D %*% theta_true + rnorm(n)
  result = list(D,Y)
  names(result) = c("D","Y")
  return(result)
}

OLS_data_t = function(n,theta_true){
  k = length(theta_true)
  D = matrix(0,n,k)
  for (i in 1:n){
    index = sample(1:k,1)
    D[i,index] = 1
  }
  Y = D %*% theta_true + rt(n,df=10)
  result = list(D,Y)
  names(result) = c("D","Y")
  return(result)
}
