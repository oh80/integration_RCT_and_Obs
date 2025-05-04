main <- function(){
  # read data
  path <- here::here("01_data", "data", "1d_0501_1.obj")
  data <- readRDS(path)
  
  X  <- data$X
  Y  <- data$Y
  Z  <- data$Z
  ID <- data$ID
  
  # MCMC
  samples <- run_MCMC(X, Y, Z, ID, iter=1000, burn_in=200)
}


run_MCMC <- function(X, Y, Z, ID, iter=1000, burn_in=200){
  # set data
  X_O <- X[ID=="O"]
  X_R <- X[ID=="R"]
  Y_O <- Y[ID=="O"]
  Y_R <- Y[ID=="R"]
  Z_O <- Z[ID=="O"]
  Z_R <- Z[ID=="R"]
  
  n_O <- length(Y_O)
  n_R <- length(Y_R)
  
  # store samples objects
  l_g_list   <- matrix(0, nrow=iter-burn_in)
  l_tau_list <- matrix(0, nrow=iter-burn_in)
  l_b_list   <- matrix(0, nrow=iter-burn_in)
  g_O_list   <- matrix(0, nrow=iter-burn_in, ncol=n_O)
  g_R_list   <- matrix(0, nrow=iter-burn_in, ncol=n_R)
  tau_O_list <- matrix(0, nrow=iter-burn_in, ncol=n_O)
  tau_R_list <- matrix(0, nrow=iter-burn_in, ncol=n_R)
  b_O_list   <- matrix(0, nrow=iter-burn_in, ncol=n_O)
  sig_list   <- matrix(0, nrow=iter-burn_in)
  
  # initial value
  l_g   <- 2
  l_tau <- 2
  l_b   <- 2
  
  g_O   <- rep(1, n_O)
  g_R   <- rep(1, n_R)
  tau_O <- rep(1, n_O)
  tau_R <- rep(1, n_R)
  b_O   <- rep(0, n_O)
  
  sig <- 1
  
  # hyper params
  sig_hyper <- 10
  
  alpha_l <- 5 
  beta_l  <- 5
  
  alpha_sig <- 5
  beta_sig  <- 5
  
  # run mcmc
  for(t in 1:iter){
    # g
    # tau
    # l
    # sig
    
    #save
    if(t > burn_in){
      j = t - burn_in
      l_g_list[j]   <- l_g
      l_tau_list[j] <- l_tau
      l_b_list[j]   <- l_b
      
      g_O_list[j,]   <- g_O
      g_R_list[j,]   <- g_R
      tau_O_list[j,] <- tau_O
      tau_R_list[j,] <- tau_R
      b_O_list[j,]   <- b_O
      
      sig_list[j] <- sig
    }
    
    if(t%%100 == 0){
      messagae = paste0(t, " iter has been done!")
      print(messagae)
    }
  }
}


main()