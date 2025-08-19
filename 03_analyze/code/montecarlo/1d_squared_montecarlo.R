main <- function(){
  # settings
  num_repeat = 50
  
  ## data settings
  nO <- 200
  nR <- 50
  sig <- 1
  
  ## MCMC settings
  iter    <- 500
  burn_in <- 200
  
  desctiption = ""
  
  # experiment in Parallel 
  cores <- parallel::detectCores()-1
  clusters <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(clusters)
  
  res_list <- foreach(i=1:num_repeat) %dopar% {
    source(here::here("01_data", "code", "1d_square_generate.R"))
    
    # data generation
    seed = i
    train_data <- generate_data(seed, nO, nR, sig)
    
    test_X <- seq(-2, 2, by=0.01)|> as.matrix()
    test_data <- list("X"= test_X,
                      "true_HTE"=test_X^2)
    
    # inference 
    source(here::here("03_analyze", "code", "montecarlo", "proposal.R"))
    proposal_res <- proposal_inference(train_data, test_data, iter, burn_in)
    
    source(here::here("03_analyze", "code", "montecarlo", "pool.R"))
    pool_res     <- pool_inference(train_data, test_data, iter, burn_in)
    
    # combine
    result <- dplyr::bind_rows(proposal_res, pool_res)
  }
  
  stopCluster(clusters)
  
  # save
  source(here::here("03_analyze", "code", "montecarlo", "utils.R"))
  res_list |> bind_df_list() |> save_result("1d_squared")
}


bind_df_list <- function(df_list){
  df <- data.frame()
  for(i in 1:length(df_list)){
    df <- df |> dplyr::bind_rows(df_list[[i]])
  }
  return(df)
}


main()