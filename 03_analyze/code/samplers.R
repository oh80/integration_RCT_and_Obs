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

compute_cross_kernel_mat <- function(X1, X2, l, eta){
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  kernel_mat <- matrix(0, nrow = n1, ncol = n2)

  for(i in 1:n1){
    for(j in 1:n2){
      sq_dist <- sum((X1[i, ] - X2[j, ])^2)
      kernel_mat[i, j] <- eta^2 * exp(-1/2 * sq_dist / l^2)
    }
  }
  return(kernel_mat)
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
  eta_new <- truncnorm::rtruncnorm(n=1, a=0, b=2, mean = eta, sd=sig_eta)
  
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
  l_new <- truncnorm::rtruncnorm(n=1, a=0, b=2, mean = l, sd=sig_hyper)
  
  small_mat <- 1e-05 * diag(length(GP))
  
  K     <- compute_kernel_mat(X, l, eta) + small_mat
  K_new <- compute_kernel_mat(X, l_new, eta) + small_mat
  
  log_like     <- -1/2 * (get_logdet(K) + t(GP) %*% chol_solve(K) %*% GP)
  log_like_new <- -1/2 * (get_logdet(K_new) + t(GP) %*% chol_solve(K_new) %*% GP)
  
  # compute log ratio
  log_r <- log_like_new - log_like                                          # likelihood
          + log(dunif(l_new, min=0, max=2))                                 # prior
          - log(dunif(l, min=0, max=2))
          - log(truncnorm::dtruncnorm(l_new, a=0, mean = l, sd=sig_hyper))
          + log(truncnorm::dtruncnorm(l, a=0, mean = l, sd=sig_hyper))      # proposal
  
  log_U <- log(runif(n=1))
  if(log_r > log_U){
    l <- l_new
    }
      
  return(l)
}


renew_g <- function(X, Y, Z, l, eta, sig_R, sig_O, tau, b_O, Obs_flag){
  n <- nrow(X)
  
  small_mat <- 1e-05 * diag(n)
  K_g <- compute_kernel_mat(X, l, eta) + small_mat
  
  HTE <- tau
  HTE[Obs_flag] <- HTE[Obs_flag] + b_O
  
  precision_vec <- rep(1/sig_R, n)
  precision_vec[Obs_flag] <- 1/sig_O
  precision_mat <- diag(precision_vec)
  A <- precision_mat + chol_solve(K_g)
  B <- precision_mat %*% (Y - Z * HTE)
  
  g <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> 
    as.vector()
  return(g)
}


#renew_tau <- function(X, Y, Z, W,  l, eta, sig, g, b_O, Obs_flag){
renew_tau <- function(X, Y, Z, l, eta, sig_R, sig_O, g, b_O, Obs_flag){
  n <- nrow(X)
  
  small_mat <- 1e-05 * diag(n)
  K_tau <- compute_kernel_mat(X, l, eta) + small_mat
  
  b <- rep(0, length(Y))
  b[Obs_flag] <- b[Obs_flag] + b_O
  
  Z_left  <- matrix(Z, nrow = n, ncol = n, byrow = TRUE)
  Z_right <- matrix(Z, nrow = n, ncol = n, byrow = FALSE)
  
  precision_vec <- rep(1/sig_R, n)
  precision_vec[Obs_flag] <- 1/sig_O
  precision_mat <- diag(precision_vec)
  A <- Z_left * precision_mat * Z_right + chol_solve(K_tau)
  B <- as.matrix(Z_left * precision_mat) %*% (Y - g - Z * b)
  
  
  tau <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> 
    as.vector()
  return(tau)
}


renew_l_b <- function(l, GP_O,  X_O, eta, sig_hyper, alpha_l, beta_l){
  l_new <- truncnorm::rtruncnorm(n=1, a=0, b=2, mean = l, sd=sig_hyper)
  
  small_mat_O <- 1e-05 * diag(length(GP_O))
  
  K_g_O     <- compute_kernel_mat(X_O, l, eta)     + small_mat_O
  K_g_O_new <- compute_kernel_mat(X_O, l_new, eta) + small_mat_O
  
  log_like <- -1/2 * (get_logdet(K_g_O) + t(GP_O) %*% chol_solve(K_g_O) %*% (GP_O))
  log_like_new <- -1/2 * (get_logdet(K_g_O_new) + t(GP_O) %*% chol_solve(K_g_O_new) %*% (GP_O))
  
  # compute log ratio
  log_r <- log_like_new - log_like                                        # likelihood
    + log(dunif(l_new, min=0, max=2))                                     # prior
    - log(dunif(l, min=0, max=2))
    - log(truncnorm::dtruncnorm(l_new, a=0, b=2, mean = l, sd=sig_hyper)) 
    + log(truncnorm::dtruncnorm(l, a=0, b=2, mean = l, sd=sig_hyper))     # proposal 
  
  
  log_U <- log(runif(n=1))
  if(log_r > log_U){
    l <- l_new
  }
  
  return(l)
}


# renew_mu_b <- function(Y_O, Z_O, g_O, tau_O, b_tilde_O, sig_O, sig_mu_b = 100){
#   post_precision <- sum(Z_O^2) / sig_O + 1 / sig_mu_b
#   sig_sq_n <- 1 / post_precision
#   
#   resid <- Y_O - g_O - Z_O * (tau_O + b_tilde_O)
#   mu_n <- sig_sq_n * (sum(Z_O * resid) / sig_O)
# 
#   mu_b <- rnorm(n = 1, mean = mu_n, sd = sqrt(sig_sq_n))
#   
#   return(mu_b)
# }

renew_mu_b <- function(Y_O, Z_O, g_O, tau_O, b_tilde_O, sig_O, W_O, mu_b, sig_mu_b = 100){
  mu_b_new <- runif(n=1, min = mu_b - 0.1, max = mu_b + 0.1)
  
  n <- length(Y_O)
  log_like_new <- 0
  log_like_old <- 0
  
  # 後でpurrrに
  for(i in 1:n){
    if(Z_O[n] == 1){
      model_new <- g_O[n] + tau_O[n] + b_tilde_O[n] + mu_b_new
      model_old <- g_O[n] + tau_O[n] + b_tilde_O[n] + mu_b
      log_like_new <- log_like_new + dnorm(Y_O[n], mean = model_new, sig_O, log = TRUE) * W_O[n]
      log_like_old <- log_like_old + dnorm(Y_O[n], mean = model_old, sig_O, log = TRUE) * W_O[n]
    }
  }
  log_r <- log_like_new - log_like_old        # likelihood
          + dnorm(mu_b_new, mean = 0, sd = 0.1)
          - dnorm(mu_b, mean = 0, sd = 0.1)     # prior
  
  log_U <- log(runif(n=1))
  if(log_r > log_U){
    mu_b <- mu_b_new
  }
  return(mu_b)
}



renew_b_tilde <- function(X, Y, Z, l, eta, sig, g, tau, mu_b){
  n <- nrow(X)
  
  small_mat <- 1e-05 * diag(n)
  K_tau     <- compute_kernel_mat(X, l, eta) + small_mat
  
  Z_left <- matrix(Z, nrow = n, ncol = n, byrow = TRUE)
  Z_right <- matrix(Z, nrow = n, ncol = n, byrow = FALSE)
  
  A <- Z_left * (1/sig * diag(n)) * Z_right + chol_solve(K_tau)
  B <- as.matrix(Z_left * (1/sig * diag(n))) %*% (Y - g - Z * (tau + mu_b))
  
  b <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> 
    as.vector()
  return(b)
}


renew_b <- function(X, Y, Z, l, eta, sig, g, tau){
  n <- nrow(X)
  
  small_mat <- 1e-05 * diag(n)
  K_tau     <- compute_kernel_mat(X, l, eta) + small_mat
  
  Z_left <- matrix(Z, nrow = n, ncol = n, byrow = TRUE)
  Z_right <- matrix(Z, nrow = n, ncol = n, byrow = FALSE)
  
  shrink_penalty = 5
  A <- Z_left * (1/sig * diag(n)) * Z_right + chol_solve(K_tau) + diag(shrink_penalty, n)
  B <- as.matrix(Z_left * (1/sig * diag(n))) %*% (Y - g - Z * tau)
  
  b <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> 
    as.vector()
  return(b)
}


renew_sig <- function(X, Y, Z, g, tau, b_O, Obs_flag, nu_0=10, sig_0=10){
  n <- length(Y)
  HTE <- tau 
  HTE[Obs_flag] <-  HTE[Obs_flag] + b_O
  
  # RCT
  RCT_flag <- c(Obs_flag == 0) |> as.integer()
  n_R <- sum(RCT_flag)
  
  nu_new_R <- nu_0 + n_R
  ss_n_R <- nu_0 * sig_0 + sum((Y[RCT_flag] - g[RCT_flag] - HTE[RCT_flag]*Z[RCT_flag])^2)
  sig_R <- 1/rgamma(n = 1, nu_new_R/2, ss_n_R/2)

  # Obs
  n_O <- sum(Obs_flag)
  
  nu_new_O <- nu_0 + n_O
  ss_n_O <- nu_0 * sig_0 + sum((Y[Obs_flag] - g[Obs_flag] - HTE[Obs_flag]*Z[Obs_flag])^2)
  sig_O <- 1/rgamma(n = 1, nu_new_O/2, ss_n_O/2)
  
  return(list(sig_R, sig_O))
}


#### samplers for multi-task version ####

compute_block_mat_inv <- function(K_R, K_O, K_RO) {
  n_R <- nrow(K_R)
  n_O <- nrow(K_O)

  inv_K_R <- chol_solve(K_R + 1e-05 * diag(n_R))
  
  K_OR <- t(K_RO)
  S <- K_O - K_OR %*% inv_K_R %*% K_RO
  inv_S <- chol_solve(S + 1e-03 * diag(n_O))

  P <- inv_S
  N <- -inv_K_R %*% K_RO %*% inv_S
  M <- inv_K_R + inv_K_R %*% K_RO %*% inv_S %*% K_OR
  
  inv_K_g_top <- cbind(M, N)
  inv_K_g_bottom <- cbind(t(N), P)
  inv_K_g <- rbind(inv_K_g_top, inv_K_g_bottom)
  
  return(inv_K_g)
}

# K_ORの処理は後で確認
renew_g_MT <- function(X, Y, Z, l_O, l_R, eta_O, eta_R,
                       sig, tau, b_O, rho, Obs_flag, RCT_flag){
  X_O <- X[Obs_flag,] |> as.matrix()
  X_R <- X[RCT_flag,] |> as.matrix()
  
  # caliculate multi-task kernel matirx
  n <- nrow(X)
  n_O <- sum(Obs_flag)
  n_R <- sum(RCT_flag)
  
  K_g_O  <- compute_kernel_mat(X_O, l_O, eta_O)
  K_g_R  <- compute_kernel_mat(X_R, l_R, eta_R)
  K_g_RO <- 1/2 * (compute_cross_kernel_mat(X_R, X_O, l_R, eta_R) + 
                   compute_cross_kernel_mat(X_R, X_O, l_O, eta_O)) * rho
  
  K_g <- matrix(0, nrow = n, ncol = n)
  K_g[1:n_R, 1:n_R]         <- K_g_R
  K_g[(n_R+1):n, (n_R+1):n] <- K_g_O
  K_g[1:n_R, (1+n_R):n]     <- K_g_RO
  K_g[(1+n_R):n, 1:n_R]     <- t(K_g_RO)
  
  B <- matrix(rho, nrow = n, ncol = n)
  B[1:n_R, 1:n_R] <- 1
  B[(n_R+1):n, (n_R+1):n] <- 1
  
  small_mat <- 1e-05 * diag(n)
  #K_g <- (K_g * B) + small_mat
  K_g <- K_g + small_mat
  
  K_g_inv <- compute_block_mat_inv(K_g_R, K_g_O, K_g_RO)
  
  # sampling renewed g
  HTE <- tau
  HTE[Obs_flag] <- HTE[Obs_flag] + b_O
  
  A <- 1/sig * diag(n) + K_g_inv
  B <- (1/sig * diag(n)) %*% (Y - Z * HTE)
  
  g <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> 
    as.vector()
  return(g)
}
