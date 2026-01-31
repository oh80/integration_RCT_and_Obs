source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "predict.R"))
source(here::here("03_analyze", "code", "GP_MCMC.R"))

proposal_inference <- function(train_data, test_data, iter, burn_in){
  X  <- train_data$X
  Y  <- train_data$Y
  Z  <- train_data$Z
  ID <- train_data$ID
  
  W <- compute_weights(train_data)
  
  samples  <- run_MCMC(X, Y, Z, W, ID, iter, burn_in)
  pred_res <- compute_pred_and_CI(train_data, test_data, samples, method="proposal")
  
  num_test <- dim(pred_res)[1]
  
  output <- pred_res |> as.data.frame() |> 
    dplyr::mutate(method   = "proposal",
                  test_X    = test_data$X, 
                  true_HTE =test_data$true_HTE)
  return(output)
}
