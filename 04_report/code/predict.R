data_path <- here::here("01_data", "data", "0602", "1d_n550_1.obj")
data <- readRDS(data_path)

proposal_path <- here::here("03_analyze","result", "0606","proposal_1d_n550_1_1.obj")
proposal_res <- readRDS(proposal_path)

samples <- proposal_res$samples

dim(samples$tau)

test_X <- seq(-2, 2, by=0.01)

source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "utils.R"))


compute_train_test_kernel <- function(X, X2, l, eta) {
  n1 <- nrow(X)
  n2 <- nrow(X2)
  
  K <- matrix(0, nrow = n1, ncol = n2)
  
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      sqdist <- sum((X[i, ] - X2[j, ])^2)
      exponent <- -1/2 * sqdist / (l^2)
      K[i, j] <- eta^2 * exp(exponent)
    }
  }
  return(K)
}

get_pred_dist_samples <- function(MCMC_samples,train_X, test_X){
  # extract samples about tau
  tau_samples <- MCMC_samples$tau
  l_samples <- MCMC_samples$l_tau
  eta_samples <- MCMC_samples$eta_tau
  
  num_samples     <- nrow(tau_samples)
  num_test_points <- nrow(test_X)
  
  # store object
  pred_dist_samples <- matrix(NA, nrow=num_samples, ncol=num_test_points )
  
  # computation pred
  for(i in 1:num_samples){
    small_mat <- 1e-05 * diag(dim(train_X)[1])
    K_train       <- compute_kernel_mat(train_X, l_samples[i], eta_samples[i])
    K_train_test  <- compute_train_test_kernel(train_X, test_X, l_samples[i], eta_samples[i])
    
    try(pred_dist_samples[i,] <- t(K_train_test) %*% chol_solve(K_train+small_mat) %*% tau_samples[i,])
    
    if(i%%50==0){
      message <- paste0(i, " iter has been done!")
      print(message)
    }
  }
  
  return(pred_dist_samples)
}

pred_samples <- get_pred_dist_samples(samples, data$X|> as.matrix(), test_X|> as.matrix())

pred_mean <- apply(pred_samples,2, mean)

calc_CI_boound <- function(vec){
  return(quantile(vec, probs=c(0.025, 0.975)))
}
pred_CI_bound <- apply(pred_samples, 2, calc_CI_boound)

true_HTE <- test_X^2
num_test <- length(true_HTE)

df_for_plot <- data.frame("X"=rep(test_X, 4),
                          "HTE"=c(pred_CI_bound[1,],pred_CI_bound[2,],pred_mean,true_HTE),
                          "label"=c(rep("lower",num_test),rep("upper",num_test),rep("mean",num_test),rep("true",num_test)))
plot <- ggplot2::ggplot(data=df_for_plot,
                        mapping = ggplot2::aes(x=X, y=HTE, color=label)) +
  ggplot2::geom_point()
plot
