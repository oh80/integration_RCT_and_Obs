source(here::here("04_report", "code", "predict.R"))
source(here::here("04_report", "code", "utils.R"))

main <- function(){
  #setting
  data_date    <- "0724"
  data_name    <- "1d_n550_1.obj"
  analyze_date <- "0724"
  analyze_name <- "both_1d_n550_1_1.obj"
  method <- "both"
  
  # read data and MCMC samples
  data_path <- here::here("01_data", "data", data_date, data_name)
  data <- readRDS(data_path)
  
  analyze_path <- here::here("03_analyze", "result", analyze_date, analyze_name)
  MCMC_res <- readRDS(analyze_path)
  samples  <- MCMC_res$samples
  
  # prepare test-data
  test_data <- prepare_test_data(data_name)
  
  # compute pred mean and 95%CI
  pred_result <- compute_pred_and_CI(data, test_data, samples, method)
  
  #compute metrics
  metrics_df <- coumpute_metrics(test_data, pred_result)
  save_metrics(metrics_df, analyze_name)
}


get_rmse <- function(y_true, y_pred){
  mse <- mean((y_true - y_pred)^2)
  rmse <- mse |> sqrt()
  return(rmse)
}

get_coverage_ratio <- function(y_true, upper, lower){
  num_pred <- length(y_true)
  in_intervel <- 0
  
  for(i in 1:num_pred){
    if(lower[i] < y_true[i] & y_true[i] < upper[i]){
      in_intervel <- in_intervel + 1
    }
  }
  coverage_ratio <- in_intervel/num_pred
  return(coverage_ratio)
}


coumpute_metrics <- function(test_data, pred_result){
  true_HTE <- test_data$true_HTE
  
  pred_mean <- pred_result$mean
  CI_lower  <- pred_result$CI_lower
  CI_upper  <- pred_result$CI_upper
  
  # metrics
  rmse <- get_rmse(true_HTE, pred_mean)
  coverage <- get_coverage_ratio(true_HTE, CI_upper, CI_lower)
  
  output <- data.frame("rmse"=rmse, "coverage"=coverage)
  return(output)
}


main()