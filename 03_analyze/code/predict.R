source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "utils.R"))


compute_train_test_kernel <- function(X, X2, l, eta) {
  n1 <- nrow(X)
  n2 <- nrow(X2)
  
  K <- matrix(0, nrow = n1, ncol = n2)
  
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      sqdist <- sum((X[i, ] - X2[j, ])^2)
      exponent <- -1/2 * sqdist / (l^2)
      K[i, j] <- eta^2 * exp(exponent)
    }
  }
  return(K)
}


get_pred_dist_samples <- function(MCMC_samples,train_X, test_X){
  # extract samples about tau
  tau_samples <- MCMC_samples$tau
  l_samples   <- MCMC_samples$l_tau
  eta_samples <- MCMC_samples$eta_tau
  
  num_samples     <- nrow(tau_samples)
  num_test_points <- nrow(test_X)
  
  # store object
  pred_mean_samples <- matrix(NA, nrow=num_samples, ncol=num_test_points)
  
  # computation pred
  print("prediction start!")
  for(i in 1:num_samples){
    small_mat <- 1e-05 * diag(dim(train_X)[1])
    K_train       <- compute_kernel_mat(train_X, l_samples[i], eta_samples[i])
    K_train_test  <- compute_train_test_kernel(train_X, test_X, l_samples[i], eta_samples[i])
    K_test        <- compute_kernel_mat(test_X, l_samples[i], eta_samples[i])
    
    pred_mean_samples[i,] <- t(K_train_test) %*% chol_solve(K_train+small_mat) %*% tau_samples[i,]
    
    if(i%%50==0){
      message <- paste0(i, " iter has been done!")
      print(message)
    }
  }
  print("prediction end!")
  
  return(pred_mean_samples)
}


compute_pred_and_CI <- function(train_data, test_data, samples, method){
  if(method == "RCT"){
    idx = c(train_data$ID=="R")
  }else if(method == "observation"){
    idx = c(train_data$ID=="O")
  }else{
    idx = rep(TRUE, length(train_data$X))
  }
  
  train_X <- train_data$X[idx] |> as.matrix()
  test_X  <- test_data$X |> as.matrix()
  
  pred_dist_samples <- get_pred_dist_samples(samples, train_X, test_X)
  pred_mean <- apply(pred_dist_samples, 2, mean)
  
  CI_upper <- c()
  CI_lower <- c()
  
  for(i in 1:length(pred_mean)){
    CI_bound <- quantile(pred_dist_samples[,i], probs = c(0.025, 0.975))
    CI_lower[i] <- CI_bound[1]
    CI_upper[i] <- CI_bound[2]
  }
  
  result <- list("mean"=pred_mean, "CI_upper"=CI_upper, "CI_lower"=CI_lower)
  return(result)
}


# functions for make test data
squeared_CATE <- function(X){
  return(X^2)
}


Relu2_CATE <- function(X){
  
  Relu <- function(x){
    if(x > 0){
      return(x)
    }else{
      return(0)
    }}
  
  output <- sapply(X, Relu)^2
  return(output)
}


make_test_data <- function(data_name, data){
  test_X <- seq(-2, 2, by=0.01)|> as.matrix()
  
  # Dimitriou setting
  if(stringr::str_detect(data_name,  "Dimitriou")){
    true_HTE <- 1+test_X+test_X^2
  }

  # original setting
  else{
    if(stringr::str_detect(data_name,  "squared_n")){
      true_HTE <- squeared_CATE(test_X)
      
    }else if(stringr::str_detect(data_name,  "Relu2")){
      true_HTE <- Relu2_CATE(test_X)
      
    }else if(stringr::str_detect(data_name,  "squared_and_root")){
      true_HTE <- sq_and_root_CATE(test_X)
    }}
  
  test_data <- list("X"=test_X, "true_HTE"=true_HTE)
  return(test_data)
}