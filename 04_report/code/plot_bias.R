source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "predict.R"))
source(here::here("04_report", "code", "plot.R"))

main <- function(){
  #setting
  data_date    <- "0820"
  data_name    <- "1d_inconstant_bias_n550_1.obj"
  analyze_date <- "0926"
  analyze_name <- "proposal_1d_inconstant_bias_n550_1_1.obj"
  
  # read data and MCMC samples
  data_path <- here::here("01_data", "data", data_date, data_name)
  data <- readRDS(data_path)
  test_X  <- seq(-2, 2, by=0.01)|> as.matrix()
  
  analyze_path <- here::here("03_analyze", "result", analyze_date, analyze_name)
  MCMC_res <- readRDS(analyze_path)
  samples  <- MCMC_res$samples
  
  pred_bias_samples <- compute_bias_pred(data, samples, test_X)
  true_bias <- get_true_bias(test_X)
  plot <- plot_bias(pred_bias_samples, true_bias, test_X)
  
  save_plot(plot, paste0(analyze_name, "_bias"))
}


compute_bias_pred <- function(data, samples, test_X){
  train_X <- data$X[data$ID=="O"] |> as.matrix()
  
  # extract samples about tau
  b_samples <- samples$b_O
  l_samples   <- samples$l_b
  eta_samples <- samples$eta_b
  mu_b_sample <- samples$mu_b
  
  num_samples     <- nrow(b_samples)
  num_test_points <- nrow(test_X)
  
  # store object
  pred_mean_samples <- matrix(NA, nrow=num_samples, ncol=num_test_points)
  
  # computation pred
  for(i in 1:num_samples){
    small_mat <- 1e-05 * diag(dim(train_X)[1])
    K_train       <- compute_kernel_mat(train_X, l_samples[i], eta_samples[i])
    K_train_test  <- compute_train_test_kernel(train_X, test_X, l_samples[i], eta_samples[i])
    K_test        <- compute_kernel_mat(test_X, l_samples[i], eta_samples[i])
    
    M <- t(K_train_test) %*% chol_solve(K_train+small_mat)
    
    # pred mean
    b_tilde_O <- b_samples[i,] - mu_b_sample[i]
    pred_mean_samples[i,] <- M %*% b_tilde_O + mu_b_sample[i]
    if(i%%50==0){
      message <- paste0(i, " iter has been done!")
      print(message)
    }
  }

  return(pred_mean_samples)
}


get_true_bias <- function(test_X, alpha=0.5, S=10000){
  num_test_points <- nrow(test_X)
  
  Bias <- rep(NA, num_test_points)
  U_e <- runif(S, -1, 1)
  
  # Monte Carlo 
  for(l in 1:num_test_points){
    x <- test_X[l]
    d <- 1+alpha*exp(-x^2)
    U <- d*U_e
    
    # Z=1
    w1 <- 1/(1+exp(-2*U))
    w1 <- w1/sum(w1)
    
    # Z=0
    w0 <- 1/(1+exp(2*U))
    w0 <- w0/sum(w0)
    
    # bias 
    Bias[l] <- sum(w1*U) - sum(w0*U)
  }
  
  return(Bias)
}


plot_bias <- function(pred_samples, true_bias, test_X){
  pred_bias_mean <- apply(pred_samples, 2, mean)
  
  df_for_plot <- data.frame("X" = rep(test_X,2),
                        "bias"=c(pred_bias_mean, true_bias),
                        "label"=c(rep("pred", length(test_X)),
                                  rep("true", length(test_X))))
  bias_plot <- ggplot2::ggplot(data=df_for_plot,
                               mapping=ggplot2::aes(x=X, y=bias, color=label)) + 
    ggplot2::geom_point(size=0.5)
  return(bias_plot)
}


main()