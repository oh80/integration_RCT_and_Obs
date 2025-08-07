#install.packages('Matching')
library(Matching)
source(here::here("02_build", "code", "utils.R"))

main <- function(){
  # read raw data
  data("lalonde")
  raw_data <- lalonde
  set.seed(428)
  
  # sample size
  n_train <- 300
  n_test  <- 145
  n_exp   <- 100  # set n_exp + n_obs < n_train
  n_obs   <- 100
  
  # clean RCT
  splitted_data <- train_test_split(raw_data, n_train, n_test)
  train_data <- splitted_data[[1]] 
  test_data  <- splitted_data[[2]] 
  
  # build processed data
  exp_and_obs <- exp_obs_split(train_data, n_exp, n_obs)
  
  # save
  train_data  |> build(drop_cols=FALSE) |> save_data(name="lalonde_train")
  test_data   |> build(drop_cols=FALSE) |> save_data(name="lalonde_test")
  exp_and_obs |> build(drop_cols=TRUE)  |> save_data(name="lalonde_mix_train")
}


train_test_split <- function(data, n_train, n_test){
  data <- data |> dplyr::mutate(ID = rep("R", nrow(data)))
  train_data <- data |> dplyr::sample_n(size=n_train)
  train_idx <- rownames(train_data) |> as.integer()
  
  test_data <- data[-train_idx,] |> dplyr::sample_n(size=n_test)
  
  return(list(train_data, test_data))
}


build <- function(data, drop_cols){
  # drop some cols
  if(drop_cols == TRUE){
    data <- data |> dplyr::select(-educ)
  }
  
  # make X, Y, Z, ID list
  Y  <- data$re78 |> as.matrix()
  X  <- data |> dplyr::select(-re78, -treat, -ID) |> as.matrix()
  Z  <- data$treat |> as.matrix()
  ID <- data$ID
  
  output <- list("X"=X, "Y"=Y, "Z"=Z, "ID"=ID)
  return(output)
}


exp_obs_split <- function(data, n_exp, n_obs){
  # exp data contain 10~12 educ people
  exp_data <- data |> dplyr::filter(10<=educ & educ<=12) |> 
    dplyr::sample_n(size=n_exp)
  exp_idx <- rownames(exp_data) |> as.integer()
  
  # obs data 
  high_educ_treated_people <- data[-exp_idx,] |> dplyr::filter(10<=educ & educ<=12) |> 
    dplyr::filter(treat==1)
  high_educ_treated_idx <- rownames(high_educ_treated_people) |> as.integer()
  
  other_people <-  data[-c(exp_idx, high_educ_treated_idx),] 
  
  # sample
  high_educ_treated_people <- high_educ_treated_people |> dplyr::sample_n(n_obs*0.6)
  other_people             <- other_people |> dplyr::sample_n(n_obs*0.4)
  
  obs_data <- dplyr::bind_rows(high_educ_treated_people, other_people) |> 
    dplyr::mutate(ID=rep("O", n_obs))
  mix_data <- dplyr::bind_rows(exp_data, obs_data)
  
  return(mix_data)
}


main()