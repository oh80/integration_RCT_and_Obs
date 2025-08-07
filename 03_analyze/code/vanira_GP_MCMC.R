source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "utils.R"))

main <- function(){
  # settings
  Date     <- "0807"
  location <- "02_build" # 01_data or 02_build
  use_data <- "both" # RCT or observation or both
  data_name <- "lalonde_train_1.obj"
  
  seed    <- 42
  iter    <- 10
  burn_in <- 5
  desctiption <- ""
  
  # read data
  data_path <- here::here(location, "data", Date, data_name)
  data <- readRDS(data_path)
  data <- extract_use_data(data, use_data)
  
  # MCMC
  set.seed(seed)
  samples <- run_MCMC(data, iter=iter, burn_in=burn_in)
  
  # save result
  data_info    <- data$info
  analyze_info <- list(seed=seed, iter=iter, burn_in=burn_in, desctiption=desctiption)
  
  result <- list(samples=samples, data_info=data_info, analyze_info=analyze_info)
  result |> save_result(use_data, data_name)
}


extract_use_data <- function(data, use_data){
  X  <- data$X |> matrix()
  Y  <- data$Y |> matrix()
  Z  <- data$Z |> matrix()
  ID <- data$ID
  
  if(use_data == "RCT"){
    data = list("X" = X[ID=="R",], Y = Y[ID=="R"], Z = Z[ID=="R"])
    
  }else if(use_data == "observation"){
    data = list("X" = X[ID=="O",], Y = Y[ID=="O"], Z = Z[ID=="O"])
    
  }else if(use_data == "both"){
    data = data
    
  }else{
    print(use_data, " is not expected.")
  }
  
  return(data)
}



run_MCMC <- function(data, iter=1000, burn_in=200){
  # set data
  X <- data$X |> as.matrix()
  Y <- data$Y |> as.matrix()
  Z <- data$Z |> as.matrix()
  
  n <- length(Y)
  
  # store samples objects
  l_g_list   <- matrix(0, nrow=iter-burn_in)
  l_tau_list <- matrix(0, nrow=iter-burn_in)
  l_b_list   <- matrix(0, nrow=iter-burn_in)
  
  eta_g_list   <- matrix(0, nrow=iter-burn_in)
  eta_tau_list <- matrix(0, nrow=iter-burn_in)
  
  g_list    <- matrix(0, nrow=iter-burn_in, ncol=n)
  tau_list  <- matrix(0, nrow=iter-burn_in, ncol=n)
  
  sig_list  <- matrix(0, nrow=iter-burn_in)
  
  # initial value
  l_g   <- 1
  l_tau <- 1
  
  eta_g   <- 1
  eta_tau <- 1
  
  g      <- rep(1, n)
  tau    <- rep(1, n)
  
  sig <- 1
  
  # hyper params
  sig_hyper <- 0.1
  
  alpha_l <- 1
  beta_l  <- 1
  
  alpha_eta <- 1
  beta_eta  <- 1
  sig_eta   <- 0.1
  
  alpha_sig <- 1
  beta_sig  <- 1

  # run mcmc
  for(t in 1:iter){
    # g
    eta_g <- renew_eta(g, eta_g, X, l_g, sig_eta, alpha_eta, beta_eta)
    l_g   <- renew_l(l_g, g, X, eta_g, sig_hyper, alpha_l, beta_l)
    g     <- renew_g(X, Y, Z, l_g, eta_g, sig, tau, b_O=0, Obs_flag=rep(FALSE, n))
    
    # tau
    eta_tau <- renew_eta(tau, eta_tau, X, l_tau, sig_eta, alpha_eta, beta_eta)
    l_tau   <- renew_l(l_tau, tau, X, eta_tau, sig_hyper, alpha_l, beta_l)
    tau     <- renew_tau(X, Y, Z, l_tau, eta_tau, sig, g, b_O=0, Obs_flag=rep(FALSE, n))
    
    # sig
    sig <- renew_sig(X, Y, Z, g, tau, b_O=0, Obs_flag=rep(FALSE, n))
    
    #save
    if(t > burn_in){
      j <- t - burn_in
      l_g_list[j]   <- l_g
      l_tau_list[j] <- l_tau
      
      eta_g_list[j]   <- eta_g
      eta_tau_list[j] <- eta_tau
      
      g_list[j, ] <- g
      tau_list[j, ] <- tau
      
      sig_list[j] <- sig
    }
    
    if(t%%50 == 0){
      message = paste0(t, " iter has been done!")
      print(message)
    }
  }
  
  # output
  samples <- list("l_g"=l_g_list, "l_tau"=l_tau_list,
                  "eta_g"=eta_g_list, "eta_tau"=eta_tau_list, 
                  "tau"=tau_list, "g"=g_list, "sig"=sig_list)
  return(samples)
}


main()