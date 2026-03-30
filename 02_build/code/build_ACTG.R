#install.packages("speff2trial")
library(speff2trial)
source(here::here("02_build", "code", "utils.R"))


main <- function(){
  set.seed(42)
  n_test <- 300
  n_RCT <- 20
  n_Obs <- 200
  data("ACTG175", package = "speff2trial")
  
  df <- ACTG175 |> dplyr::select(c(pidnum, cd496, treat, cd80, cd40, karnof)) |> 
    tidyr::drop_na() 
  
  # train test split
  df <- df[sample(nrow(df)),]
  train_df <- df[1:(nrow(df)-n_test),]
  test_df  <- df[(nrow(df)-n_test+1):nrow(df),]
  
  # prepare RCT and Obs
  RCT_df <- train_df |> dplyr::filter((cd80/cd40 < 1.5) & cd40 > 200)
  RCT_df <- RCT_df[1:n_RCT,] |> dplyr::mutate(ID = "R")
  RCT_df <- RCT_df|> dplyr::mutate(ID = "R")
  
  prior_model_df <- train_df |> dplyr::filter(!pidnum %in% RCT_df$pidnum)
  prior_model_df <- prior_model_df[1:50,]
  prior_model <- lm(cd496 ~ splines::ns(cd40, df=3) + splines::ns(cd80, df=3) + karnof, data = prior_model_df)
  
  Obs_df <- train_df |> dplyr::filter(!pidnum %in% c(RCT_df$pidnum, prior_model_df$pidnum))
  Obs_df <- Obs_df |> dplyr::mutate(y_pred = predict(prior_model, Obs_df)) 
  down_id_list <- c()
  for(i in 1:nrow(Obs_df)){
    if(Obs_df$treat[i] == 0){
      if(Obs_df$y_pred[i] > 200){
        p <- rbinom(1, 1, 0.8)
        if(p == 1){
          down_id_list <- c(down_id_list, Obs_df$pidnum[i])
        }
      }
    }
    else if(Obs_df$treat[i] == 1){
      if(Obs_df$y_pred[i] < 200){
        p <- rbinom(1, 1, 0.8)
        if(p == 1){
          down_id_list <- c(down_id_list, Obs_df$pidnum[i])
        }
      }
    }
  }
  Obs_df <- Obs_df |> 
    dplyr::filter(!pidnum %in% down_id_list) |> 
    dplyr::select(-y_pred) |> 
    dplyr::mutate(ID = "O")
  Obs_df <- Obs_df[1:n_Obs,]

  confound_train_df <- rbind(RCT_df, Obs_df) |> dplyr::select(-karnof)
  
  # save
  build_data <- list(train_df = train_df, test_data = test_df, 
                     confound_train_df = confound_train_df)
  save_data(build_data, "ACTG")
  
}

main()

  
# df <- ACTG175
# lm_df <- df |> dplyr::select(-c(cd420, pidnum, arms, r, treat, cd820, cens, preanti, offtrt, days))
# # #
# lm_fit <- lm(cd496~ ., data = lm_df)
# summary(lm_fit)

# hist(lm_df$cd80)
# hist(lm_df$cd40)
# hist(lm_df$karnof)
#
# cor(lm_df$cd80,lm_df$cd40)
# cor(lm_df$cd80, lm_df$karnof)

