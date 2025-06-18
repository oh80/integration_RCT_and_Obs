source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "utils.R"))

main <- function(){
  # settings
  Date <- "0602"
  data_name <- "1d_n550_1.obj"
  seed    <- 42
  iter    <- 5
  burn_in <- 2
  desctiption <- ""
  
  # read data
  data_path <- here::here("01_data", "data", Date, data_name)
  data <- readRDS(data_path)
  
  X  <- data$X
  Y  <- data$Y
  Z  <- data$Z
  ID <- data$ID
  
  # MCMC
  set.seed(seed)
  samples <- run_MCMC(X, Y, Z, ID, iter=iter, burn_in=burn_in)
  
  # save result
  data_info    <- data$info
  analyze_info <- list(seed=seed, iter=iter, burn_in=burn_in, desctiption=desctiption)
  
  result <- list(samples=samples, data_info=data_info, analyze_info=analyze_info)
  result |> save_result("proposal", data_path)
}


run_MCMC <- function(X, Y, Z, ID, iter=1000, burn_in=200){
  # set data
  X <- X |> matrix()
  Y <- Y |> matrix()
  Z <- Z |> matrix()
  
  X_O <- X[ID=="O"] |> matrix()
  X_R <- X[ID=="R"] |> matrix()
  Y_O <- Y[ID=="O"] |> matrix()
  Y_R <- Y[ID=="R"] |> matrix()
  Z_O <- Z[ID=="O"] |> matrix()
  Z_R <- Z[ID=="R"] |> matrix()
  
  n_O <- length(Y_O)
  n_R <- length(Y_R)
  
  Obs_flag <- c(ID == "O")
  RCT_flag <- c(ID == "R")
  
  # store samples objects
  l_g_O_list <- matrix(0, nrow=iter-burn_in)
  l_g_R_list <- matrix(0, nrow=iter-burn_in)
  l_tau_list <- matrix(0, nrow=iter-burn_in)
  l_b_list   <- matrix(0, nrow=iter-burn_in)
  
  eta_g_O_list <- matrix(0, nrow=iter-burn_in)
  eta_g_R_list <- matrix(0, nrow=iter-burn_in)
  eta_tau_list <- matrix(0, nrow=iter-burn_in)
  eta_b_list   <- matrix(0, nrow=iter-burn_in)
  
  g_list    <- matrix(0, nrow=iter-burn_in, ncol=(n_O+n_R))
  tau_list  <- matrix(0, nrow=iter-burn_in, ncol=(n_O+n_R))
  b_O_list  <- matrix(0, nrow=iter-burn_in, ncol=n_O)
  
  sig_list  <- matrix(0, nrow=iter-burn_in)
  
  # initial value
  l_g_O  <- 2
  l_g_R  <- 2
  l_tau  <- 2
  l_b    <- 2
  
  eta_g_O  <- 3
  eta_g_R  <- 3
  eta_tau  <- 3
  eta_b    <- 3
  
  g     <- rep(1, n_O+n_R)
  tau   <- rep(1, n_O+n_R)
  b_O   <- rep(0, n_O)
  
  sig <- 1
  
  rho <- 0.5
  
  # hyper params
  sig_hyper <- 1
  
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
    eta_g_O <- renew_eta(g[Obs_flag], eta_g_O, X_O, l_g_O, sig_eta, alpha_eta, beta_eta)
    eta_g_R <- renew_eta(g[RCT_flag], eta_g_R, X_R, l_g_R, sig_eta, alpha_eta, beta_eta)
    l_g_O   <- renew_l(l_g_O, g[Obs_flag], X_O, eta_g_O, sig_hyper, alpha_l, beta_l)
    l_g_R   <- renew_l(l_g_R, g[RCT_flag], X_R, eta_g_R, sig_hyper, alpha_l, beta_l)
    #rho   <- renew_rho()
    g     <- renew_g_MT(X, Y, Z, l_g_O, l_g_R, eta_g_O, eta_g_R, sig, tau, b_O, rho, Obs_flag, RCT_flag)
    
    # tau
    eta_tau <- renew_eta(tau, eta_tau, X, l_tau, sig_eta, alpha_eta, beta_eta)
    l_tau   <- renew_l(l_tau, tau, X, eta_tau, sig_hyper, alpha_l, beta_l)
    tau     <- renew_tau(X, Y, Z, l_tau, eta_tau, sig, g, b_O, Obs_flag)
    
    # b
    eta_b <- renew_eta(b_O, eta_b, X_O, l_b, sig_eta, alpha_eta, beta_eta)
    l_b <- renew_l_b(l_b, tau[Obs_flag], X_O, eta_b, sig_hyper, alpha_l, beta_l)
    b_O <- renew_b(X_O, Y_O, Z_O, l_b, eta_b, sig, g[Obs_flag], tau[Obs_flag])
    
    # sig
    sig <- renew_sig(X, Y, Z, g, tau, b_O, Obs_flag)
    
    #save
    if(t > burn_in){
      j <- t - burn_in
      l_g_O_list[j] <- l_g_O
      l_g_R_list[j] <- l_g_R
      l_tau_list[j] <- l_tau
      l_b_list[j]   <- l_b
      
      eta_g_O_list[j] <- eta_g_O
      eta_g_R_list[j] <- eta_g_R
      eta_tau_list[j] <- eta_tau
      eta_b_list[j]   <- eta_b
      
      g_list[j, ]   <- g
      tau_list[j, ] <- tau
      b_O_list[j,]  <- b_O
      
      sig_list[j] <- sig
    }
    
    if(t%%50 == 0){
      message = paste0(t, " iter has been done!")
      print(message)
    }
  }
  
  # output
  samples <- list("l_g_O"=l_g_O_list, "l_g_R"=l_g_R_list, "l_tau"=l_tau_list, 
                  "l_b"=l_b_list, "eta_g_O"=eta_g_O_list, "eta_g_R"=eta_g_R_list,
                  "eta_tau"=eta_tau_list, "eta_b"=eta_b_list, "tau"=tau_list,
                  "g"=g_list, "b_O"=b_O_list, "sig"=sig_list)
  return(samples)
}


main()