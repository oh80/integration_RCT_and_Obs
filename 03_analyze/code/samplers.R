#### utils ####
compute_kernel_mat <- function(X, l, eta){
  n <- nrow(X)
  base_mat <- matrix(0, nrow=n, ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      if(i > j){
        base <- -1/2 * sum((X[i,] - X[j,])^2)
        base_mat[i, j] <- base
        base_mat[j, i] <- base
      }}}
  log_K <- base_mat / l^2
  K <- eta^2 * exp(log_K)
  diag(K) <- eta^2
  
  return(K)
}


get_logdet <- function(K){
  return(msos::logdet(K))
}


chol_solve <- function(A) {
  L    <- chol(A)
  invA <- chol2inv(L)
  return(invA)
}


#### sampler functions ####

renew_eta <- function(GP, eta, X, l, sig_eta, alpha_eta, beta_eta){
  # sample from proposal distribution
  eta_new <- truncnorm::rtruncnorm(n=1, a=0, mean = eta, sd=sig_eta)
  
  small_mat <- 1e-05 * diag(length(GP))
  K     <- compute_kernel_mat(X, l, eta) + small_mat
  K_new <- compute_kernel_mat(X, l, eta_new) + small_mat
  
  log_like     <- -1/2 * (get_logdet(K) + t(GP) %*% chol_solve(K) %*% GP)
  log_like_new <- -1/2 * (get_logdet(K_new) + t(GP) %*% chol_solve(K_new) %*% GP)
  
  # compute log ratio
  log_r <- log_like_new - log_like                       # likelihood
  + log(dgamma(eta_new, alpha_eta, beta_eta))                  # prior
  - log(dgamma(eta, alpha_eta, beta_eta))
  - log(truncnorm::dtruncnorm(eta_new, a=0, mean = eta, sd=sig_eta)) 
  + log(truncnorm::dtruncnorm(eta, a=0, mean = eta, sd=sig_eta))     # proposal 
  
  log_U <- log(runif(n=1))
  if(log_r > log_U){
    eta <- eta_new
  }
  
  return(eta)
}


renew_l <- function(l, GP, X, eta, sig_hyper, alpha_l, beta_l){
  # sample from proposal distribution
  l_new <- truncnorm::rtruncnorm(n=1, a=0, mean = l, sd=sig_hyper)
  
  small_mat <- 1e-05 * diag(length(GP))
  
  K     <- compute_kernel_mat(X, l, eta) + small_mat
  K_new <- compute_kernel_mat(X, l_new, eta) + small_mat
  
  log_like     <- -1/2 * (get_logdet(K) + t(GP) %*% chol_solve(K) %*% GP)
  log_like_new <- -1/2 * (get_logdet(K_new) + t(GP) %*% chol_solve(K_new) %*% GP)
  
  # compute log ratio
  log_r <- log_like_new - log_like                           # likelihood
          + log(dgamma(l_new, alpha_l, beta_l))                  # prior
          - log(dgamma(l, alpha_l, beta_l))
          - log(truncnorm::dtruncnorm(l_new, a=0, mean = l, sd=sig_hyper)) 
          + log(truncnorm::dtruncnorm(l, a=0, mean = l, sd=sig_hyper))     # proposal 
  
  log_U <- log(runif(n=1))
  if(log_r > log_U){
    l <- l_new
    }
      
  return(l)
}



renew_g <- function(X, Y, Z, l, eta, sig, tau, b_O, Obs_flag){
  n <- nrow(X)
  
  small_mat <- 1e-05 * diag(n)
  K_g <- compute_kernel_mat(X, l, eta) + small_mat
  
  HTE <- tau
  HTE[Obs_flag] <- HTE[Obs_flag] + b_O
  
  A <- 1/sig * diag(n) + chol_solve(K_g)
  B <- (1/sig * diag(n)) %*% (Y - Z * HTE)
  
  g <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> 
    as.vector()
  return(g)
}



renew_tau <- function(X, Y, Z, l, eta, sig, g, b_O, Obs_flag){
  n <- nrow(X)
  
  small_mat <- 1e-05 * diag(n)
  K_tau <- compute_kernel_mat(X, l, eta) + small_mat
  
  b <- rep(0, length(Y))
  b[Obs_flag] <- b[Obs_flag] + b_O
  
  Z_left <- matrix(Z, nrow = n, ncol = n, byrow = TRUE)
  Z_right <- matrix(Z, nrow = n, ncol = n, byrow = FALSE)
  
  A <- Z_left * (1/sig * diag(n)) * Z_right + chol_solve(K_tau)
  B <- as.matrix(Z_left * (1/sig * diag(n))) %*% (Y - g - Z * b)
  
  tau <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> 
    as.vector()
  return(tau)
}


renew_l_b <- function(l, GP_O,  X_O, eta, sig_hyper, alpha_l, beta_l){
  # sample from proposal distribution
  l_new <- truncnorm::rtruncnorm(n=1, a=0, mean=l, sd=sig_hyper)
  
  small_mat_O <- 1e-05 * diag(length(GP_O))
  
  K_g_O     <- compute_kernel_mat(X_O, l, eta)     + small_mat_O
  K_g_O_new <- compute_kernel_mat(X_O, l_new, eta) + small_mat_O
  
  log_like <- -1/2 * (get_logdet(K_g_O) + t(GP_O) %*% chol_solve(K_g_O) %*% (GP_O))
  log_like_new <- -1/2 * (get_logdet(K_g_O_new) + t(GP_O) %*% chol_solve(K_g_O_new) %*% (GP_O))
  # compute log ratio
  log_r <- log_like_new - log_like                           # likelihood
  + log(dgamma(l_new, alpha_l, beta_l))                  # prior
  - log(dgamma(l, alpha_l, beta_l))
  - log(truncnorm::dtruncnorm(l_new, a=0, sd=sig_hyper)) 
  + log(truncnorm::dtruncnorm(l, a=0, sd=sig_hyper))     # proposal 
  
  log_U <- log(runif(n=1))
  if(log_r > log_U){
    l <- l_new
  }
  
  return(l)
}


renew_b <- function(X, Y, Z, l, eta, sig, g, tau){
  n <- nrow(X)
  
  small_mat <- 1e-05 * diag(n)
  K_tau     <- compute_kernel_mat(X, l, eta) + small_mat
  
  Z_left <- matrix(Z, nrow = n, ncol = n, byrow = TRUE)
  Z_right <- matrix(Z, nrow = n, ncol = n, byrow = FALSE)
  
  A <- Z_left * (1/sig * diag(n)) * Z_right + chol_solve(K_tau)
  B <- as.matrix(Z_left * (1/sig * diag(n))) %*% (Y - g - Z * tau)
  
  b <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> 
    as.vector()
  return(b)
}


renew_sig <- function(X, Y, Z, g, tau, b_O, Obs_flag, nu_0=10, sig_0=10){
  n <- length(Y)
  nu_new <- nu_0 + n
  
  HTE <- tau 
  HTE[Obs_flag] <-  HTE[Obs_flag] + b_O
  
  ss_n <- nu_0 * sig_0 + sum((Y - g - HTE*Z)^2)
  sig <- 1/rgamma(n = 1, nu_new/2, ss_n/2)
  return(sig)
}
