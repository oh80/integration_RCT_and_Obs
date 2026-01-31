library(doParallel)


main <- function(){
  # settings
  num_repeat = 100
  
  ## data settings
  nO <- 1000
  nR <- 150
  sig <- 1
  alpha <- 0.5
  
  ## MCMC settings
  iter    <- 500
  burn_in <- 200
  
  desctiption = ""
  
  # experiment in Parallel 
  cores <- parallel::detectCores()-1
  clusters <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(clusters)
  
  res_list <- foreach(i=1:num_repeat) %dopar% {
    source(here::here("01_data", "code", "1d_inconstant_bias.R"))
    
    # data generation
    seed = i
    train_data <- generate_data(seed, nO, nR, sig, alpha)
    
    test_X <- seq(-2, 2, by=0.01)|> as.matrix()
    test_data <- list("X"= test_X,
                      "true_HTE"=test_X^2)
    
    data_list <- list("train_data" = train_data, "test_data" = test_data)
    data_list |> saveRDS(here::here("01_data", "data", "inconstant_bias", paste0(i,".obj")))
    
    # inference 
    source(here::here("03_analyze", "code", "montecarlo", "proposal.R"))
    proposal_res <- proposal_inference(train_data, test_data, iter, burn_in)
    
    source(here::here("03_analyze", "code", "montecarlo", "2step_proposal.R"))
    step2_proposal_res <- step2_proposal_inference(train_data, test_data, iter, burn_in)
    
    #source(here::here("03_analyze", "code", "montecarlo", "pool.R"))
    #pool_res     <- pool_inference(train_data, test_data, iter, burn_in)
    
    # combine
    result <- dplyr::bind_rows(step2_proposal_res, pool_res)
    #result <- step2_proposal_res
  }
  
  stopCluster(clusters)
  
  # save
  source(here::here("03_analyze", "code", "montecarlo", "utils.R"))
  res_list |> bind_df_list() |> save_result("1d_inconstant_bias")
}


bind_df_list <- function(df_list){
  df <- data.frame()
  for(i in 1:length(df_list)){
    df <- df |> dplyr::bind_rows(df_list[[i]])
  }
  return(df)
}


main()