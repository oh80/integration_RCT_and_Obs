source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "predict.R"))
source(here::here("03_analyze", "code", "utils.R"))

main <- function(){
  Date <- "1107"
  data_name <- "1d_constant_bias_n250_1.obj"
  location  <- "01_data"     # 01_data or 02_build
  
  seed    <- 42
  iter    <- 100
  burn_in <- 200
  desctiption <- ""
  
  # read data
  data_path <- here::here(location, "data", Date, data_name)
  data <- readRDS(data_path)
  
  X  <- data$X
  Y  <- data$Y
  Z  <- data$Z
  ID <- data$ID
    
  # step-1 (insample GP)
  min_max_list <- get_RCT_sapport(X, ID)
  insample_flag <- extract_insample_id(X, min_max_list)
  X_in  <- X[insample_flag]
  Y_in  <- Y[insample_flag]
  Z_in  <- Z[insample_flag]
  ID_in <- ID[insample_flag]
  
  set.seed(seed)
  step_1_samples <- run_MCMC(X_in, Y_in, Z_in, ID_in, iter, burn_in, step = 1)
  
  # step-2 (full-sample GP)
  b_O_plugin <- compute_b_O_pred(step_1_samples, X_in, X, ID_in, ID)
  step_2_samples <- run_MCMC(X, Y, Z, ID, iter, burn_in, step = 2, b_O_plugin)
  
  # predict
  if(location == "01_data"){
    test_data <- make_test_data(data_name, data)
  }else{
    test_data <- data$test_X
  }
  pred_res <- compute_pred_and_CI(data, test_data, step_2_samples, method="proposal")
  
  # save result
  data_info    <- data$info
  analyze_info <- list(seed=seed, iter=iter, burn_in=burn_in, desctiption=desctiption)
  
  result <- list(samples=step_2_samples, pred_res=pred_res, b_O_plugin = b_O_plugin, 
                 data_info=data_info, analyze_info=analyze_info)
  #result |> save_result("2step-proposal", data_name)
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
  
  output <- cbind(min_list, max_list)
  return(output)
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


run_MCMC <- function(X, Y, Z, ID, iter, burn_in, step, b_O_plugin = NA){
  # set data
  X <- X |> as.matrix()
  Y <- Y |> as.matrix()
  Z <- Z |> as.matrix()
  
  X_O <- X[ID=="O",] |> as.matrix()
  X_R <- X[ID=="R",] |> as.matrix()
  Y_O <- Y[ID=="O",] |> as.matrix()
  Y_R <- Y[ID=="R",] |> as.matrix()
  Z_O <- Z[ID=="O",] |> as.matrix()
  Z_R <- Z[ID=="R",] |> as.matrix()
  
  n_O <- length(Y_O)
  n_R <- length(Y_R)
  
  Obs_flag <- c(ID == "O")
  
  # store samples objects
  l_g_list   <- matrix(0, nrow=iter-burn_in)
  l_tau_list <- matrix(0, nrow=iter-burn_in)
  l_b_list   <- matrix(0, nrow=iter-burn_in)
  
  eta_g_list   <- matrix(0, nrow=iter-burn_in)
  eta_tau_list <- matrix(0, nrow=iter-burn_in)
  eta_b_list   <- matrix(0, nrow=iter-burn_in)
  
  g_list    <- matrix(0, nrow=iter-burn_in, ncol=(n_O+n_R))
  tau_list  <- matrix(0, nrow=iter-burn_in, ncol=(n_O+n_R))
  b_O_list  <- matrix(0, nrow=iter-burn_in, ncol=n_O)
  
  #sig_list  <- matrix(0, nrow=iter-burn_in)
  sig_R_list  <- matrix(0, nrow=iter-burn_in)
  sig_O_list  <- matrix(0, nrow=iter-burn_in)
  
  # initial value
  l_g   <- 0.5
  l_tau <- 0.5
  l_b   <- 0.5
  
  eta_g   <- 1
  eta_tau <- 1
  eta_b   <- 1
  
  g      <- rep(1, n_O+n_R)
  tau   <- rep(1, n_O+n_R)
  b_O   <- rep(0, n_O)
  b_tilde_O <- rep(0, n_O)
  
  sig_R <- 1
  sig_O <- 1
  
  if(step == 2){
    s <- sample(seq(1,iter - burn_in), size=1)
    mu_b <- b_O_plugin$mean[[s]]
    L_b  <- b_O_plugin$L[[s]]
    ep   <- ep <- rnorm(n = length(mu_b))
    b_O <- mu_b + L_b %*% ep
  }
  
  # hyper params
  sig_hyper <- 0.5
  
  alpha_l <- 2
  beta_l  <- 2
  
  alpha_eta <- 2
  beta_eta  <- 2
  sig_eta   <- 1
  
  alpha_sig <- 2
  beta_sig  <- 2
  
  # run mcmc
  for(t in 1:iter){
    # g
    eta_g <- renew_eta(g, eta_g, X, l_g, sig_eta, alpha_eta, beta_eta)
    l_g   <- renew_l(l_g, g, X, eta_g, sig_hyper, alpha_l, beta_l)
    g     <- renew_g(X, Y, Z, l_g, eta_g, sig_R, sig_O, tau, b_O, Obs_flag)
    
    # tau
    eta_tau <- renew_eta(tau, eta_tau, X, l_tau, sig_eta, alpha_eta, beta_eta)
    l_tau   <- renew_l(l_tau, tau, X, eta_tau, sig_hyper, alpha_l, beta_l)
    tau     <- renew_tau(X, Y, Z, l_tau, eta_tau, sig_R, sig_O, g, b_O, Obs_flag)
    
    # b
    if(step == 1){
      eta_b <- renew_eta(b_O, eta_b, X_O, l_b, sig_eta, alpha_eta, beta_eta)
      l_b   <- renew_l(l_b, b_O, X_O, eta_b, sig_hyper, alpha_l, beta_l)
      b_O   <- renew_b(X_O, Y_O, Z_O, l_b, eta_b, sig_O, g[Obs_flag], tau[Obs_flag])
    }else{
      s <- sample(seq(1,iter - burn_in), size=1)
      mu_b <- b_O_plugin$mean[[s]]
      L_b  <- b_O_plugin$L[[s]]
      ep   <- ep <- rnorm(n = length(mu_b))
      b_O <- mu_b + L_b %*% ep
    }
    
    # sig
    sig <- renew_sig(X, Y, Z, g, tau, b_O, Obs_flag)
    sig_R <- sig[[1]]
    sig_O <- sig[[2]]
    
    #save
    if(t > burn_in){
      j <- t - burn_in
      l_g_list[j]   <- l_g
      l_tau_list[j] <- l_tau
      l_b_list[j]   <- l_b
      
      eta_g_list[j]   <- eta_g
      eta_tau_list[j] <- eta_tau
      eta_b_list[j]   <- eta_b
      
      g_list[j, ] <- g
      tau_list[j, ] <- tau
      b_O_list[j,]   <- b_O
      
      sig_R_list[j] <- sig_R
      sig_O_list[j] <- sig_O
    }
    
    if(t%%50 == 0){
      message = paste0("step ", step, ": ", t, " iter has been done!")
      print(message)
    }
  }
  
  # output
  samples <- list("l_g"=l_g_list, "l_tau"=l_tau_list, "l_b"=l_b_list,
                  "eta_g"=eta_g_list, "eta_tau"=eta_tau_list, "eta_b"=eta_b_list,
                  "tau"=tau_list, "g"=g_list,
                  "b_O"=b_O_list, "sig_R"=sig_R_list, "sig_O"=sig_O_list)
  return(samples)
}


compute_b_O_pred <- function(samples, X_in, X, ID_in, ID){
  X_in <- X_in[ID_in == "O",,drop=FALSE] |> as.matrix()
  X    <- X[ID == "O", , drop=FALSE] |> as.matrix()
  
  b_O_in_sample <- samples$b_O  |> as.matrix()
  l_b_sample    <- samples$l_b  
  eta_b_sample  <- samples$eta_b 
  
  small_mat <- 1e-05 * diag(dim(X_in)[1])
  
  n_sample <- length(l_b_sample)
  n <- dim(X)[1]
  b_O_mean_sample <- list()
  b_O_L_sample    <- list()
  
  # extrapolate to observational region
  for(t in 1:n_sample){
    K_train       <- compute_kernel_mat(X_in, l_b_sample[t], eta_b_sample[t])
    K_train_test  <- compute_train_test_kernel(X_in, X, l_b_sample[t], eta_b_sample[t])
    K_test        <- compute_kernel_mat(X, l_b_sample[t], eta_b_sample[t])
    
    M <- t(K_train_test) %*% chol_solve(K_train+small_mat)
    b_O_mean_sample[[t]] <- M %*% b_O_in_sample[t,]
    
    b_O_Sigma <- K_test - M %*% K_train_test
    ev <- eigen(b_O_Sigma + diag(1e-05, nrow(b_O_Sigma)), symmetric = TRUE)
    ev$values[ev$values < 0] <- 0
    b_O_L_sample[[t]] <- ev$vectors %*% diag(sqrt(ev$values))
  }
  
  b_O_sample <- list("L"=b_O_L_sample, "mean"=b_O_mean_sample)
  return(b_O_sample)
}


#main()