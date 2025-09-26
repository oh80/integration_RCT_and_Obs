source(here::here("04_report", "code", "utils.R"))

main <- function(){
  #setting
  data_date    <- "0820"
  data_name    <- "1d_inconstant_bias_n550_1.obj"
  analyze_date <- "0926"
  analyze_name <- "proposal_1d_inconstant_bias_n550_1_1.obj"
  method <- "proposal"
  
  # read data and MCMC samples
  data_path <- here::here("01_data", "data", data_date, data_name)
  data <- readRDS(data_path)
  
  analyze_path <- here::here("03_analyze", "result", analyze_date, analyze_name)
  analyze_res <- readRDS(analyze_path)
  pred_result <- analyze_res$pred_res
  
  # prepare test-data
  test_data <- prepare_test_data(data_name, data$info)
  
  # plot
  plot <- plot_result(test_data, pred_result)
  save_plot(plot, analyze_name)
}


prepare_test_data <- function(data_name, info){
  if(stringr::str_detect(data_name,  "narrow_1d")){
    test_X <- seq(-3, 3, by=0.01)|> as.matrix()
  }else{
    test_X <- seq(-2, 2, by=0.01)|> as.matrix()
  }
  
  # Dimitriou setting
  if(stringr::str_detect(data_name,  "Dimitriou")){
    true_HTE <- 1+test_X+test_X^2
  }
  
  # original setting
  if(stringr::str_detect(data_name, "squared_n")){
    true_HTE <- test_X^2
  }
  if(stringr::str_detect(data_name, "Relu")){
    true_HTE <- Relu2_CATE(test_X)
  }
  if(stringr::str_detect(data_name, "squared_and_root")){
    true_HTE <- sq_and_root_CATE(test_X)
  }
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
    ggplot2::geom_point(size=0.5)
  plot
  return(plot)
}


main()