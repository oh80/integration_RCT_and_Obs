source(here::here("04_report", "code", "predict.R"))
source(here::here("04_report", "code", "utils.R"))

main <- function(){
  #setting
  data_date    <- "0724"
  data_name    <- "1d_n550_1.obj"
  analyze_date <- "0724"
  analyze_name <- "RCT_1d_n550_1_1.obj"
  method <- "RCT"
  
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
  
  # plot
  plot <- plot_result(test_data, pred_result)
  save_plot(plot, analyze_name)
}


prepare_test_data <- function(data_name){
  test_X <- seq(-2, 2, by=0.01)|> as.matrix()
  
  # Dimitriou setting
  if(stringr::str_detect(data_name,  "Dimitriou")){
    true_HTE <- 1+test_X+test_X^2
  }
  
  # original setting
  else{
    true_HTE <- test_X^2
  }
  
  test_data <- list("X"=test_X, "true_HTE"=true_HTE)
  return(test_data)
}


plot_result <- function(test_data, pred_result){
  test_X   <- test_data$X
  true_HTE <- test_data$true_HTE
  
  pred_mean <- pred_result$mean
  CI_lower  <- pred_result$CI_lower
  CI_upper  <- pred_result$CI_upper
  
  num_test <- length(true_HTE)
  df_for_plot <- data.frame("X"=rep(test_X, 4),
                            "HTE"=c(CI_lower, CI_upper, pred_mean, true_HTE),
                            "label"=c(rep("lower",num_test),rep("upper",num_test),rep("mean",num_test),rep("true",num_test)))
  plot <- ggplot2::ggplot(data=df_for_plot,
                          mapping = ggplot2::aes(x=X, y=HTE, color=label)) +
    ggplot2::geom_point()
  return(plot)
}


main()