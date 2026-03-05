source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "predict.R"))
source(here::here("03_analyze", "code", "utils.R"))

main <- function(){
  Date <- "1107"
  data_name <- "1d_constant_bias_n250_1.obj"
  location  <- "01_data"     # 01_data or 02_build
  
  seed    <- 42
  iter    <- 50
  burn_in <- 20
  desctiption <- ""
  
  # read data
  data_path <- here::here(location, "data", Date, data_name)
  data <- readRDS(data_path)
  
  X  <- data$X |> as.matrix()
  Y  <- data$Y |> as.vector()
  Z  <- data$Z |> as.vector()
  ID <- data$ID |> as.vector()
  
  # step-1 (insample Z=0)
  min_max_list <- get_RCT_sapport(X, ID)
  insample_flag <- extract_insample_id(X, min_max_list)
  X_in  <- X[insample_flag, , drop=FALSE]
  Y_in  <- Y[insample_flag]
  Z_in  <- Z[insample_flag]
  ID_in <- ID[insample_flag]
  
  set.seed(seed)
  step_1_samples <- run_MCMC(X_in[Z_in == 0, , drop=FALSE], Y_in[Z_in == 0], Z_in[Z_in == 0], ID_in[Z_in == 0], iter, burn_in, step = 1)
  
  # Step 1のプラグイン作成
  idx_in_O <- (Z_in == 0) & (ID_in == "O")
  b_g_O_plugin <- compute_b_g_O_pred(step_1_samples, X_in[idx_in_O, , drop=FALSE], X_in, ID_in)
  
  # step-2 (insample Z=0 or 1)
  step_2_samples <- run_MCMC(X_in, Y_in, Z_in, ID_in, iter, burn_in, step = 2, b_g_O_plugin)
  
  # Step 3のプラグイン作成（全領域への外挿）
  b_g_O_plugin <- compute_b_g_O_pred(step_1_samples, X_in[idx_in_O], X, ID)
  b_tau_O_plugin <- compute_b_tau_O_pred(step_2_samples, X_in, X, ID_in, ID)
  
  # step-3 (full-sample GP)
  step_3_samples <- run_MCMC(X, Y, Z, ID, iter, burn_in, step = 3, b_g_O_plugin, b_tau_O_plugin)
  
  # predict
  if(location == "01_data"){
    test_data <- make_test_data(data_name, data)
  }else{
    test_data <- data$test_X
  }
  pred_res <- compute_pred_and_CI(data, test_data, step_3_samples, method="proposal")
  
  # save result
  data_info    <- data$info
  analyze_info <- list(seed=seed, iter=iter, burn_in=burn_in, desctiption=desctiption)
  
  result <- list(samples=step_3_samples, pred_res=pred_res, b_g_O_plugin = b_g_O_plugin, b_tau_O_plugin = b_tau_O_plugin, 
                 data_info=data_info, analyze_info=analyze_info)
  #result |> save_result("3step-proposal", data_name)
}

get_RCT_sapport <- function(X, ID){
  X <- X |> as.matrix()
  min_list <- c()
  max_list <- c()
  n_col <- dim(X)[2]
  RCT_X <- X[c(ID == "R"), 1:n_col] |> as.matrix()
  n_RCT <- dim(RCT_X)[1]
  for(i in 1:n_col){
    min_list[i] <- min(RCT_X[1:n_RCT,i])
    max_list[i] <- max(RCT_X[1:n_RCT,i])
  }
  return(cbind(min_list, max_list))
}

extract_insample_id <- function(X, min_max_list){
  X <- X |> as.matrix()
  n     <- dim(X)[1]
  n_col <- dim(X)[2]
  insample_flag <- rep(TRUE, n)
  for(i in 1:n){
    for(d in 1:n_col){
      if(X[i,d] < min_max_list[d, 1] | min_max_list[d, 2] < X[i,d]){
        insample_flag[i] <- FALSE
      }
    }
  }
  return(insample_flag)
}

# --- Gibbs更新式を「行列演算」から「ベクトル演算」へ大幅最適化 ---

renew_g <- function(X, Y, Z, l, eta, sig_R, sig_O, tau, b_g_O, b_tau_O, Obs_flag){
  n <- nrow(X)
  small_mat <- 1e-05 * diag(n)
  K_g <- compute_kernel_mat(X, l, eta) + small_mat
  
  HTE <- tau
  HTE[Obs_flag] <- HTE[Obs_flag] + b_tau_O
  
  full_b_g <- rep(0, n)
  full_b_g[Obs_flag] <- b_g_O
  
  precision_vec <- rep(1/sig_R, n)
  precision_vec[Obs_flag] <- 1/sig_O
  
  A <- diag(precision_vec) + chol_solve(K_g)
  B <- precision_vec * (Y - full_b_g - Z * HTE)
  
  g <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> as.vector()
  return(g)
}

renew_tau <- function(X, Y, Z, l, eta, sig_R, sig_O, g, b_g_O, b_tau_O, Obs_flag){
  n <- nrow(X)
  small_mat <- 1e-05 * diag(n)
  K_tau <- compute_kernel_mat(X, l, eta) + small_mat
  
  full_b_tau <- rep(0, n)
  full_b_tau[Obs_flag] <- b_tau_O
  full_b_g <- rep(0, n)
  full_b_g[Obs_flag] <- b_g_O
  
  precision_vec <- rep(1/sig_R, n)
  precision_vec[Obs_flag] <- 1/sig_O
  
  # Z^2 = Z (since Z in {0,1})
  A <- diag(Z * precision_vec * Z) + chol_solve(K_tau)
  B <- Z * precision_vec * (Y - g - full_b_g - Z * full_b_tau)
  
  tau <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> as.vector()
  return(tau)
}

renew_b_g <- function(X_O, Y_O, Z_O, l, eta, sig_O, g_O, tau_O, b_tau_O){
  n <- nrow(X_O)
  small_mat <- 1e-05 * diag(n)
  K_bg <- compute_kernel_mat(X_O, l, eta) + small_mat
  
  A <- diag(1/sig_O, n) + chol_solve(K_bg)
  B <- (1/sig_O) * (Y_O - g_O - Z_O * (tau_O + b_tau_O))
  
  b <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> as.vector()
  return(b)
}

renew_b_tau <- function(X_O, Y_O, Z_O, l, eta, sig_O, g_O, tau_O, b_g_O){
  n <- nrow(X_O)
  small_mat <- 1e-05 * diag(n)
  K_btau <- compute_kernel_mat(X_O, l, eta) + small_mat
  
  A <- diag(Z_O * (1/sig_O) * Z_O) + chol_solve(K_btau)
  B <- Z_O * (1/sig_O) * (Y_O - g_O - b_g_O - Z_O * tau_O)
  
  b <- mvtnorm::rmvnorm(n=1, mean = chol_solve(A)%*%B, sigma=chol_solve(A)) |> as.vector()
  return(b)
}

renew_sig <- function(X, Y, Z, g, tau, b_g_O, b_tau_O, Obs_flag, nu_0=10, sig_0=10){
  HTE <- tau 
  HTE[Obs_flag] <- HTE[Obs_flag] + b_tau_O
  full_g <- g
  full_g[Obs_flag] <- full_g[Obs_flag] + b_g_O
  
  RCT_flag <- !Obs_flag
  n_R <- sum(RCT_flag)
  sig_R <- 1 # dummy
  if(n_R > 0){
    nu_new_R <- nu_0 + n_R
    ss_n_R <- nu_0 * sig_0 + sum((Y[RCT_flag] - full_g[RCT_flag] - HTE[RCT_flag]*Z[RCT_flag])^2)
    sig_R <- 1/rgamma(n = 1, nu_new_R/2, ss_n_R/2)
  }
  
  n_O <- sum(Obs_flag)
  nu_new_O <- nu_0 + n_O
  ss_n_O <- nu_0 * sig_0 + sum((Y[Obs_flag] - full_g[Obs_flag] - HTE[Obs_flag]*Z[Obs_flag])^2)
  sig_O <- 1/rgamma(n = 1, nu_new_O/2, ss_n_O/2)
  
  return(list(sig_R, sig_O))
}


run_MCMC <- function(X, Y, Z, ID, iter, burn_in, step, b_g_O_plugin = NA, b_tau_O_plugin = NA){
  X <- X |> as.matrix()
  Y <- Y |> as.vector()
  Z <- Z |> as.vector()
  
  X_O <- X[ID=="O", , drop=FALSE] |> as.matrix()
  X_R <- X[ID=="R", , drop=FALSE] |> as.matrix()
  Y_O <- Y[ID=="O"]
  Y_R <- Y[ID=="R"]
  Z_O <- Z[ID=="O"]
  Z_R <- Z[ID=="R"]
  
  n_O <- length(Y_O)
  n_R <- length(Y_R)
  Obs_flag <- (ID == "O")
  
  # store samples objects
  l_g_list   <- matrix(0, nrow=iter-burn_in)
  l_tau_list <- matrix(0, nrow=iter-burn_in)
  l_b_g_list   <- matrix(0, nrow=iter-burn_in)
  l_b_tau_list   <- matrix(0, nrow=iter-burn_in)
  
  eta_g_list   <- matrix(0, nrow=iter-burn_in)
  eta_tau_list <- matrix(0, nrow=iter-burn_in)
  eta_b_g_list   <- matrix(0, nrow=iter-burn_in)
  eta_b_tau_list   <- matrix(0, nrow=iter-burn_in)
  
  g_list    <- matrix(0, nrow=iter-burn_in, ncol=(n_O+n_R))
  tau_list  <- matrix(0, nrow=iter-burn_in, ncol=(n_O+n_R))
  b_g_O_list  <- matrix(0, nrow=iter-burn_in, ncol=n_O)
  b_tau_O_list  <- matrix(0, nrow=iter-burn_in, ncol=n_O)
  
  sig_R_list  <- matrix(0, nrow=iter-burn_in)
  sig_O_list  <- matrix(0, nrow=iter-burn_in)
  
  # initial value
  l_g   <- 0.5
  l_tau <- 0.5
  l_b_g   <- 0.5
  l_b_tau <- 0.5
  
  eta_g   <- 1
  eta_tau <- 1
  eta_b_g   <- 1
  eta_b_tau <- 1
  
  g       <- rep(1, n_O+n_R)
  tau     <- rep(1, n_O+n_R)
  b_g_O   <- rep(0, n_O)
  b_tau_O <- rep(0, n_O)
  
  sig_R <- 1
  sig_O <- 1
  
  if(step == 2 | step == 3){
    s <- sample(seq(1,iter - burn_in), size=1)
    mu_b <- b_g_O_plugin$mean[[s]]
    L_b  <- b_g_O_plugin$L[[s]]
    b_g_O <- mu_b + L_b %*% rnorm(n = length(mu_b))
  }
  
  if(step == 3){
    s <- sample(seq(1,iter - burn_in), size=1)
    mu_b <- b_tau_O_plugin$mean[[s]]
    L_b  <- b_tau_O_plugin$L[[s]]
    b_tau_O <- mu_b + L_b %*% rnorm(n = length(mu_b))
  }
  
  # hyper params
  sig_hyper <- 0.5
  alpha_l <- 2; beta_l  <- 2
  alpha_eta <- 2; beta_eta  <- 2; sig_eta   <- 1
  alpha_sig <- 2; beta_sig  <- 2
  
  # run mcmc
  for(t in 1:iter){
    # g 
    eta_g <- renew_eta(g, eta_g, X, l_g, sig_eta, alpha_eta, beta_eta)
    l_g   <- renew_l(l_g, g, X, eta_g, sig_hyper, alpha_l, beta_l)
    g     <- renew_g(X, Y, Z, l_g, eta_g, sig_R, sig_O, tau, b_g_O, b_tau_O, Obs_flag)
    
    # b_g
    if(step == 1){
      eta_b_g <- renew_eta(b_g_O, eta_b_g, X_O, l_b_g, sig_eta, alpha_eta, beta_eta)
      l_b_g   <- renew_l(l_b_g, b_g_O, X_O, eta_b_g, sig_hyper, alpha_l, beta_l)
      b_g_O   <- renew_b_g(X_O, Y_O, Z_O, l_b_g, eta_b_g, sig_O, g[Obs_flag], tau[Obs_flag], b_tau_O)
    }else{
      s <- sample(seq(1,iter - burn_in), size=1)
      mu_b <- b_g_O_plugin$mean[[s]]
      L_b  <- b_g_O_plugin$L[[s]]
      b_g_O <- mu_b + L_b %*% rnorm(n = length(mu_b))|> as.vector()
    }
    
    # tau
    if(step != 1){
      eta_tau <- renew_eta(tau, eta_tau, X, l_tau, sig_eta, alpha_eta, beta_eta)
      l_tau   <- renew_l(l_tau, tau, X, eta_tau, sig_hyper, alpha_l, beta_l)
      tau     <- renew_tau(X, Y, Z, l_tau, eta_tau, sig_R, sig_O, g, b_g_O, b_tau_O, Obs_flag)
    }
    
    # b_tau
    if(step == 2){
      eta_b_tau <- renew_eta(b_tau_O, eta_b_tau, X_O, l_b_tau, sig_eta, alpha_eta, beta_eta)
      l_b_tau   <- renew_l(l_b_tau, b_tau_O, X_O, eta_b_tau, sig_hyper, alpha_l, beta_l)
      b_tau_O   <- renew_b_tau(X_O, Y_O, Z_O, l_b_tau, eta_b_tau, sig_O, g[Obs_flag], tau[Obs_flag], b_g_O)
    }else if(step == 3){
      s <- sample(seq(1,iter - burn_in), size=1)
      mu_b <- b_tau_O_plugin$mean[[s]]
      L_b  <- b_tau_O_plugin$L[[s]]
      b_tau_O <- mu_b + L_b %*% rnorm(n = length(mu_b)) |> as.vector()
    }
    
    # sig
    sig <- renew_sig(X, Y, Z, g, tau, b_g_O, b_tau_O, Obs_flag)
    sig_R <- sig[[1]]
    sig_O <- sig[[2]]
    
    #save
    if(t > burn_in){
      j <- t - burn_in
      l_g_list[j]   <- l_g
      l_tau_list[j] <- l_tau
      l_b_g_list[j]   <- l_b_g
      l_b_tau_list[j] <- l_b_tau
      
      eta_g_list[j]   <- eta_g
      eta_tau_list[j] <- eta_tau
      eta_b_g_list[j]   <- eta_b_g
      eta_b_tau_list[j]   <- eta_b_tau
      
      g_list[j, ] <- g
      tau_list[j, ] <- tau
      b_g_O_list[j,]   <- as.vector(b_g_O)
      b_tau_O_list[j,] <- as.vector(b_tau_O)
      
      sig_R_list[j] <- sig_R
      sig_O_list[j] <- sig_O
    }
    
    if(t%%50 == 0){
      message = paste0("step ", step, ": ", t, " iter has been done!")
      print(message)
    }
  }
  
  samples <- list("l_g"=l_g_list, "l_tau"=l_tau_list, "l_b_g"=l_b_g_list, "l_b_tau"=l_b_tau_list,
                  "eta_g"=eta_g_list, "eta_tau"=eta_tau_list, "eta_b_g"=eta_b_g_list, "eta_b_tau"=eta_b_tau_list,
                  "tau"=tau_list, "g"=g_list,
                  "b_g_O"=b_g_O_list, "b_tau_O"=b_tau_O_list, "sig_R"=sig_R_list, "sig_O"=sig_O_list)
  return(samples)
}


compute_b_g_O_pred <- function(samples, X_train, X_test, ID_test){
  X_train_O <- X_train |> as.matrix()
  X_test_O  <- X_test[ID_test == "O", , drop=FALSE] |> as.matrix()
  
  b_O_in_sample <- samples$b_g_O  |> as.matrix()
  l_b_sample    <- samples$l_b_g 
  eta_b_sample  <- samples$eta_b_g
  
  small_mat <- 1e-05 * diag(nrow(X_train_O))
  n_sample <- length(l_b_sample)
  b_O_mean_sample <- list()
  b_O_L_sample    <- list()
  
  for(t in 1:n_sample){
    K_train       <- compute_kernel_mat(X_train_O, l_b_sample[t], eta_b_sample[t])
    K_train_test  <- compute_train_test_kernel(X_train_O, X_test_O, l_b_sample[t], eta_b_sample[t])
    K_test        <- compute_kernel_mat(X_test_O, l_b_sample[t], eta_b_sample[t])
    
    M <- t(K_train_test) %*% chol_solve(K_train+small_mat)
    b_O_mean_sample[[t]] <- M %*% b_O_in_sample[t,]
    
    b_O_Sigma <- K_test - M %*% K_train_test
    ev <- eigen(b_O_Sigma + diag(1e-05, nrow(b_O_Sigma)), symmetric = TRUE)
    ev$values[ev$values < 0] <- 0
    b_O_L_sample[[t]] <- ev$vectors %*% diag(sqrt(ev$values))
  }
  return(list("L"=b_O_L_sample, "mean"=b_O_mean_sample))
}


compute_b_tau_O_pred <- function(samples, X_train, X_test, ID_train, ID_test){
  X_train_O <- X_train[ID_train == "O", , drop=FALSE] |> as.matrix()
  X_test_O  <- X_test[ID_test == "O", , drop=FALSE] |> as.matrix()
  
  b_O_in_sample <- samples$b_tau_O  |> as.matrix()
  l_b_sample    <- samples$l_b_tau  
  eta_b_sample  <- samples$eta_b_tau
  
  small_mat <- 1e-05 * diag(nrow(X_train_O))
  n_sample <- length(l_b_sample)
  b_O_mean_sample <- list()
  b_O_L_sample    <- list()
  
  for(t in 1:n_sample){
    K_train       <- compute_kernel_mat(X_train_O, l_b_sample[t], eta_b_sample[t])
    K_train_test  <- compute_train_test_kernel(X_train_O, X_test_O, l_b_sample[t], eta_b_sample[t])
    K_test        <- compute_kernel_mat(X_test_O, l_b_sample[t], eta_b_sample[t])
    
    M <- t(K_train_test) %*% chol_solve(K_train+small_mat)
    b_O_mean_sample[[t]] <- M %*% b_O_in_sample[t,]
    
    b_O_Sigma <- K_test - M %*% K_train_test
    ev <- eigen(b_O_Sigma + diag(1e-05, nrow(b_O_Sigma)), symmetric = TRUE)
    ev$values[ev$values < 0] <- 0
    b_O_L_sample[[t]] <- ev$vectors %*% diag(sqrt(ev$values))
  }
  return(list("L"=b_O_L_sample, "mean"=b_O_mean_sample))
}

res <- main()

x <- seq(-2, 2, 0.01)
y_true <- x^2

# 1. 縦軸の最小値・最大値を計算（グラフがはみ出さないようにするため）
y_min <- min(res$pred_res$CI_lower, y_true)
y_max <- max(res$pred_res$CI_upper, y_true)

# 2. ベースとなるplot（平均値）を描画。type="l" で線グラフにする
plot(x, res$pred_res$mean, type="l", col="blue", lwd=2, 
     ylim=c(y_min, y_max), xlab="x", ylab="Treatment Effect")

# 3. linesで残りの線を上書きしていく
lines(x, res$pred_res$CI_upper, col="blue", lty=2) # 信頼区間上限（点線）
lines(x, res$pred_res$CI_lower, col="blue", lty=2) # 信頼区間下限（点線）
lines(x, y_true, col="red", lwd=2)                 # 真の関数（赤色）

mean((y_true-res$pred_res$mean)^2)|> sqrt()