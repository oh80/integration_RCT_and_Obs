source(here::here("03_analyze", "code", "samplers.R"))

main <- function(){
  # read data
  path <- here::here("01_data", "data", "1d_0531_1.obj")
  data <- readRDS(path)
  
  X  <- data$X
  Y  <- data$Y
  Z  <- data$Z
  ID <- data$ID
  
  # MCMC
  samples <- run_MCMC(X, Y, Z, ID, iter=100, burn_in=30)
  return(samples)
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
  
  sig <- 1
  
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
    # 
    eta_g <- renew_eta(g, eta_g, X, l_g, sig_eta, alpha_eta, beta_eta)
    l_g   <- renew_l(l_g, g, X, eta_g, sig_hyper, alpha_l, beta_l)
    g     <- renew_g(X, Y, Z, l_g, eta_g, sig, tau, b_O, Obs_flag)
    
    # tau
    eta_tau <- renew_eta(tau, eta_tau, X, l_tau, sig_eta, alpha_eta, beta_eta)
    l_tau   <- renew_l(l_tau, tau, X, eta_tau, sig_hyper, alpha_l, beta_l)
    tau     <- renew_tau(X, Y, Z, l_tau, eta_tau, sig, g, b_O, Obs_flag)
    
    # b
    #eta_b <- renew_eta(b_O, eta_b, X_O, l_b, sig_eta, alpha_eta, beta_eta)
    #l_b <- renew_l_b(l_b, tau[Obs_flag], X_O, eta_b, sig_hyper, alpha_l, beta_l)
    #b_O <- renew_b(X_O, Y_O, Z_O, l_b, eta_b, sig, g[Obs_flag], tau[Obs_flag])
    
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
    }
    
    if(t%%10 == 0){
      message = paste0(t, " iter has been done!")
      print(message)
    }
  }
  
  # output
  samples <- list("l_g"=l_g_list, "l_tau"=l_tau_list, "l_b"=l_b_list,
                  "tau"=tau_list, "g"=g_list,
                  "b_O"=b_O_list, "sig"=sig_list)
  return(samples)
}

res <- main()
# 
# 
# path <- here::here("01_data", "data", "1d_0531_1.obj")
# data <- readRDS(path)
# 
# tau_est <- apply(res$tau, 2, mean)
# b_est   <- apply(res$b_O, 2, mean)
# HTE_est <- tau_est
# HTE_est[data$ID=="O"] <- HTE_est[data$ID=="O"] +b_est
# 
# plot(x=data$Tau, y=HTE_est)
# 
# plot(x=data$Tau[data$ID=="O"], y=HTE_est[data$ID=="O"])
# plot(x=data$Tau[data$ID=="R"], y=HTE_est[data$ID=="R"])
# hist(res$l_g)
# hist(res$l_b)
# 
# hist(res$l_tau)
# hist(res$sig)
# 
# hist(apply(res$b_O, 2, mean))
# 
# plot(data$X, data$Tau)
# plot(data$X, HTE_est)
