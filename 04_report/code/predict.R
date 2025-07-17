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
  pred_var_samples  <- matrix(NA, nrow=num_samples, ncol=num_test_points)
  
  # computation pred
  for(i in 1:num_samples){
    small_mat <- 1e-05 * diag(dim(train_X)[1])
    K_train       <- compute_kernel_mat(train_X, l_samples[i], eta_samples[i])
    K_train_test  <- compute_train_test_kernel(train_X, test_X, l_samples[i], eta_samples[i])
    K_test        <- compute_kernel_mat(test_X, l_samples[i], eta_samples[i])
    
    try(pred_mean_samples[i,] <- t(K_train_test) %*% chol_solve(K_train+small_mat) %*% tau_samples[i,])
    pred_cov_mat  <- K_test - t(K_train_test) %*% chol_solve(K_train+small_mat) %*% K_train_test
    pred_var_samples[i,] <- diag(pred_cov_mat)
    
    if(i%%50==0){
      message <- paste0(i, " iter has been done!")
      print(message)
    }
  }
  
  output <- list("mean"=pred_mean_samples, "var"=pred_var_samples)
  return(output)
}
