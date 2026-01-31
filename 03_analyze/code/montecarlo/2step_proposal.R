source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "predict.R"))
source(here::here("03_analyze", "code", "2step_GP_MCMC.R"))
source(here::here("03_analyze", "code", "montecarlo", "utils.R"))

step2_proposal_inference <- function(train_data, test_data, iter, burn_in){
  X  <- train_data$X
  Y  <- train_data$Y
  Z  <- train_data$Z
  ID <- train_data$ID
  
  # step-1 (insample GP)
  min_max_list <- get_RCT_sapport(X, ID)
  insample_flag <- extract_insample_id(X, min_max_list)
  X_in  <- X[insample_flag]
  Y_in  <- Y[insample_flag]
  Z_in  <- Z[insample_flag]
  ID_in <- ID[insample_flag]

  step_1_samples <- run_MCMC(X_in, Y_in, Z_in, ID_in, iter, burn_in, step = 1)
  
  # step-2 (full-sample GP)
  b_O_plugin <- compute_b_O_pred(step_1_samples, X_in, X, ID_in, ID)
  step_2_samples <- run_MCMC(X, Y, Z, ID, iter, burn_in, step = 2, b_O_plugin)
  
  # prediction
  pred_res <- compute_pred_and_CI(train_data, test_data, step_2_samples, method="proposal")
  num_test <- dim(pred_res)[1]
  
  output <- pred_res |> as.data.frame() |> 
    dplyr::mutate(method   = "2step_proposal",
                  test_X   = test_data$X, 
                  true_HTE = test_data$true_HTE)
  
  # save plugin estimator
  mcmc_sample_list <- list("b_O_plugin"=b_O_plugin, "step_2_samples"=step_2_samples,
                           "train_data"=train_data)
  save_result(mcmc_sample_list, "mcmc_sample_res_")
  
  return(output)
}
