source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "predict.R"))
source(here::here("03_analyze", "code", "utils.R"))


main <- function(){
  # settings
  Date <- "0820"
  data_name <- "1d_inconstant_bias_n550_1.obj"
  location  <- "01_data"     # 01_data or 02_build
  
  seed    <- 42
  iter    <- 500
  burn_in <- 200
  desctiption <- ""
  
  # read data
  data_path <- here::here(location, "data", Date, data_name)
  data <- readRDS(data_path)
  
  X  <- data$X
  Y  <- data$Y
  Z  <- data$Z
  ID <- data$ID
  
  # MCMC
  set.seed(seed)
  samples <- run_MCMC(X, Y, Z, ID, iter=iter, burn_in=burn_in)
  
  # predict
  if(location == "01_data"){
    test_data <- make_test_data(data_name, data)
  }else{
    test_data <- data$test_X
  }
  pred_res <- compute_pred_and_CI(data, test_data, samples, method="proposal")
  
  # save result
  data_info    <- data$info
  analyze_info <- list(seed=seed, iter=iter, burn_in=burn_in, desctiption=desctiption)
  
  result <- list(samples=samples, pred_res=pred_res, 
                 data_info=data_info, analyze_info=analyze_info)
  result |> save_result("proposal", data_name)
}


run_MCMC <- function(X, Y, Z, ID, iter=1000, burn_in=200){
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
  
  sig_list  <- matrix(0, nrow=iter-burn_in)
  
  mu_b_list  <- matrix(0, nrow=iter-burn_in)
  
  # initial value
  l_g   <- 2
  l_tau <- 2
  l_b   <- 2
  
  eta_g   <- 3
  eta_tau <- 3
  eta_b   <- 3
  
  g      <- rep(1, n_O+n_R)
  tau   <- rep(1, n_O+n_R)
  b_O   <- rep(0, n_O)
  b_tilde_O <- rep(0, n_O)
  
  sig <- 1
  
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
    g     <- renew_g(X, Y, Z, l_g, eta_g, sig, tau, b_O, Obs_flag)
    
    # tau
    eta_tau <- renew_eta(tau, eta_tau, X, l_tau, sig_eta, alpha_eta, beta_eta)
    l_tau   <- renew_l(l_tau, tau, X, eta_tau, sig_hyper, alpha_l, beta_l)
    tau     <- renew_tau(X, Y, Z, l_tau, eta_tau, sig, g, b_O, Obs_flag)
    
    # b
    eta_b <- renew_eta(b_tilde_O, eta_b, X_O, l_b, sig_eta, alpha_eta, beta_eta)
    l_b   <- renew_l_b(l_b, b_tilde_O, X_O, eta_b, sig_hyper, alpha_l, beta_l)
    mu_b  <- renew_mu_b(Y_O, Z_O, g[Obs_flag], tau[Obs_flag], b_tilde_O, sig, sig_mu_b = 100)
    b_tilde_O  <- renew_b_tilde(X_O, Y_O, Z_O, l_b, eta_b, sig, g[Obs_flag], tau[Obs_flag], mu_b)
    
    b_O <- mu_b + b_tilde_O
    
    # sig
    sig <- renew_sig(X, Y, Z, g, tau, b_O, Obs_flag)
    
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
      
      sig_list[j] <- sig
      
      mu_b_list[j] <- mu_b
    }
    
    if(t%%50 == 0){
      message = paste0(t, " iter has been done!")
      print(message)
    }
  }
  
  # output
  samples <- list("l_g"=l_g_list, "l_tau"=l_tau_list, "l_b"=l_b_list,
                  "eta_g"=eta_g_list, "eta_tau"=eta_tau_list, "eta_b"=eta_b_list,
                  "tau"=tau_list, "g"=g_list,
                  "b_O"=b_O_list, "sig"=sig_list, "mu_b"=mu_b_list)
  return(samples)
}

main()