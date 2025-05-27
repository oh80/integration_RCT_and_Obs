#### utils ####
compute_kernel_mat <- function(X, l, eta=1){
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


#### sampler functions ####

# sampler of l_g and tau_g
# renew_l <- function(l, GP_O, GP_R, X_O, X_R, sig_hyper, alpha_l, beta_l){
#   # sample from proposal distribution
#   l_new <- truncnorm::rtruncnorm(n=1, a=0, sd=sig_hyper)
#   
#   small_mat_O <- 1e-05 * diag(length(GP_O))
#   small_mat_R <- 1e-05 * diag(length(GP_R))
#   
#   K_g_O     <- compute_kernel_mat(X_O, l)     + small_mat_O
#   K_g_R     <- compute_kernel_mat(X_R, l)     + small_mat_R
#   K_g_O_new <- compute_kernel_mat(X_O, l_new) + small_mat_O
#   K_g_R_new <- compute_kernel_mat(X_R, l_new) + small_mat_R
#   
#   log_like <- -1/2 * (get_logdet(K_g_O) + t(GP_O) %*% solve(K_g_O) %*% (GP_O)) +
#     -1/2 * (get_logdet(K_g_R) + t(GP_R) %*% solve(K_g_R) %*% (GP_R))
#   log_like_new <- -1/2 * (get_logdet(K_g_O_new) + t(GP_O) %*% solve(K_g_O_new) %*% (GP_O)) +
#     -1/2 * (get_logdet(K_g_R_new) + t(GP_R) %*% solve(K_g_R_new) %*% (GP_R))
# 
#   # compute log ratio
#   log_r <- log_like_new - log_like                           # likelihood
#       + log(dgamma(l_new, alpha_l, beta_l))                  # prior
#       - log(dgamma(l, alpha_l, beta_l))
#       - log(truncnorm::dtruncnorm(l_new, a=0, sd=sig_hyper)) 
#       + log(truncnorm::dtruncnorm(l, a=0, sd=sig_hyper))     # proposal 
# 
#   log_U <- log(runif(n=1))
#   if(log_r > log_U){
#     l <- l_new
#   }
#      
#   return(l)
# }

renew_l <- function(l, GP, X, sig_hyper, alpha_l, beta_l){
  # sample from proposal distribution
  l_new <- truncnorm::rtruncnorm(n=1, a=0, sd=sig_hyper)
  
  small_mat <- 1e-05 * diag(length(GP))
  
  K     <- compute_kernel_mat(X, l) + small_mat
  K_new <- compute_kernel_mat(X, l_new) + small_mat
  
  log_like     <- -1/2 * (get_logdet(K) + t(GP) %*% solve(K) %*% GP)
  log_like_new <- -1/2 * (get_logdet(K_new) + t(GP) %*% solve(K_new) %*% GP)
  
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


# 
# renew_g <- function(X, Y, Z, l, sig, HTE){
#   n <- nrow(X)
#   
#   small_mat <- 1e-05 * diag(n)
#   K_g     <- compute_kernel_mat(X, l) + small_mat
#   
#   A <- 1/sig * diag(n) + solve(K_g)
#   B <- (1/sig * diag(n)) %*% (Y - Z * HTE)
#   
#   g <- mvtnorm::rmvnorm(n=1, mean = solve(A+small_mat)%*%B, sigma=solve(A)) |> 
#     as.vector()
#   return(g)
# }


renew_g <- function(X, Y, Z, l, sig, tau, b_O, Obs_flag){
  n <- nrow(X)
  
  small_mat <- 1e-05 * diag(n)
  K_g <- compute_kernel_mat(X, l) + small_mat
  
  HTE <- tau
  tau[Obs_flag] <- tau[Obs_flag] + b_O
  
  A <- 1/sig * diag(n) + solve(K_g)
  B <- (1/sig * diag(n)) %*% (Y - Z * HTE)
  
  g <- mvtnorm::rmvnorm(n=1, mean = solve(A+small_mat)%*%B, sigma=solve(A)) |> 
    as.vector()
  return(g)
}


# renew_tau <- function(X, Y, Z, l, sig, g, b=0){
#   n <- nrow(X)
#   
#   small_mat <- 1e-05 * diag(n)
#   K_tau     <- compute_kernel_mat(X, l) + small_mat
#   
#   Z_left <- matrix(Z, nrow = n, ncol = n, byrow = TRUE)
#   Z_right <- matrix(Z, nrow = n, ncol = n, byrow = FALSE)
#   
#   A <- Z_left %*% (1/sig * diag(n)) %*% Z_right + solve(K_tau)
#   B <- Z_left %*% (1/sig * diag(n)) %*% (Y - g - Z * b)
#   
#   g <- mvtnorm::rmvnorm(n=1, mean = solve(A+small_mat)%*%B, sigma=solve(A)) |> 
#     as.vector()
#   return(g)
#}

renew_tau <- function(X, Y, Z, l, sig, g, b_O, Obs_flag){
  n <- nrow(X)
  
  small_mat <- 1e-05 * diag(n)
  K_tau <- compute_kernel_mat(X, l) + small_mat
  
  b <- rep(0, length(Y))
  b[Obs_flag] <- b[Obs_flag] + b_O
  
  Z_left <- matrix(Z, nrow = n, ncol = n, byrow = TRUE)
  Z_right <- matrix(Z, nrow = n, ncol = n, byrow = FALSE)
  
  A <- Z_left %*% (1/sig * diag(n)) %*% Z_right + solve(K_tau)
  B <- Z_left %*% (1/sig * diag(n)) %*% (Y - g - Z * b)
  
  g <- mvtnorm::rmvnorm(n=1, mean = solve(A+small_mat)%*%B, sigma=solve(A)) |> 
    as.vector()
  return(g)
}


renew_l_b <- function(l, GP_O,  X_O, sig_hyper, alpha_l, beta_l){
  # sample from proposal distribution
  l_new <- truncnorm::rtruncnorm(n=1, a=0, sd=sig_hyper)
  
  small_mat_O <- 1e-05 * diag(length(GP_O))
  
  K_g_O     <- compute_kernel_mat(X_O, l)     + small_mat_O
  K_g_O_new <- compute_kernel_mat(X_O, l_new) + small_mat_O
  
  log_like <- -1/2 * (get_logdet(K_g_O) + t(GP_O) %*% solve(K_g_O) %*% (GP_O))
  log_like_new <- -1/2 * (get_logdet(K_g_O_new) + t(GP_O) %*% solve(K_g_O_new) %*% (GP_O))
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


renew_b <- function(X, Y, Z, l, sig, g, tau){
  n <- nrow(X)
  
  small_mat <- 1e-05 * diag(n)
  K_tau     <- compute_kernel_mat(X, l) + small_mat
  
  Z_left <- matrix(Z, nrow = n, ncol = n, byrow = TRUE)
  Z_right <- matrix(Z, nrow = n, ncol = n, byrow = FALSE)
  
  A <- Z_left %*% (1/sig * diag(n)) %*% Z_right + solve(K_tau)
  B <- Z_left %*% (1/sig * diag(n)) %*% (Y - g - Z * tau)
  
  g <- mvtnorm::rmvnorm(n=1, mean = solve(A+small_mat)%*%B, sigma=solve(A)) |> 
    as.vector()
  return(g)
}
