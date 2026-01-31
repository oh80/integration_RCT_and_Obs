# read data
train_data_path <- here::here("02_build", "data", "0807", "lalonde_train_1.obj")
train_data <- readRDS(train_data_path)

mix_data_path <- here::here("02_build", "data", "0807", "lalonde_mix_train_1.obj")
mix_data <- readRDS(mix_data_path)

test_data_path <- here::here("02_build", "data", "0807", "lalonde_test_1.obj")
test_data <- readRDS(test_data_path)


# read res
prop_res_path <- here::here("03_analyze", "result", "0807", "proposal_lalonde_mix_train_1_1.obj")
prop_res <- readRDS(prop_res_path)

exp_res_path <- here::here("03_analyze", "result", "0807", "both_lalonde_train_1_2.obj")
exp_res <- readRDS(prop_res_path)



# predict
source(here::here("04_report", "code", "predict.R"))

train_X <- mix_data$X
train_X <- train_data$X
test_X  <- test_data$X
test_X <- test_X[,-2]

samples  <- exp_res$samples
tau_samples <- samples$tau
l_samples   <- samples$l_tau
eta_samples <- samples$eta_tau

num_samples     <- nrow(tau_samples)
num_test_points <- nrow(test_X)

pred_mean_samples <- matrix(NA, nrow=num_samples, ncol=num_test_points)
for(i in 1:num_samples){
  small_mat <- 1e-05 * diag(dim(train_X)[1])
  K_train       <- compute_kernel_mat(train_X, l_samples[i], eta_samples[i])
  K_train_test  <- compute_train_test_kernel(train_X, test_X, l_samples[i], eta_samples[i])
  K_test        <- compute_kernel_mat(test_X, l_samples[i], eta_samples[i])
  print(dim(t(K_train_test)))
  print(dim(chol_solve(K_train+small_mat)))
  
  pred_mean_samples[i,] <- t(K_train_test) %*% chol_solve(K_train+small_mat) %*% tau_samples[i,]
  
  if(i%%50==0){
    message <- paste0(i, " iter has been done!")
    print(message)
  }
}

dim(tau_samples)
dim(train_data$X)

pred_mean <- apply(pred_mean_samples, 2, mean) |> as.data.frame()
write.csv(pred_mean, file = here::here("04_report","result","0807","prop.csv"))

hist(test_data$Y)
plot(samples$tau)
plot(pred_mean)
