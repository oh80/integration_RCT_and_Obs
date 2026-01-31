result1 <- readRDS(here::here("03_analyze", "result", "0120", "1d_squared_montecarlo_1.obj"))
result2 <- readRDS(here::here("03_analyze", "result", "0113", "1d_squared_montecarlo_1.obj"))

result <- result2 |> dplyr::filter(method == "pool") |> 
  dplyr::bind_rows(result1)


ICM_res <- read.csv(here::here("03_analyze", "result", "0113", "python_montecarlo","ICM_RHO_08.csv"))
kallus_res <- read.csv(here::here("03_analyze", "result", "0113", "python_montecarlo","kallus.csv"))

# coverage
coverage_df <- result |>
  dplyr::mutate(isin_intervel = ifelse(CI_lower <= true_HTE & true_HTE <= CI_upper, 1, 0)) |>
  dplyr::group_by(method, test_X) |> dplyr::summarise(coverage = mean(isin_intervel))

ICM_coverage_df <- ICM_res |> 
  dplyr::mutate(isin_intervel = ifelse(lower <= true_HTE & true_HTE <= upper, 1, 0)) |>
  dplyr::group_by(test_X) |> dplyr::summarise(coverage = mean(isin_intervel)) |>
  dplyr::mutate(method = "ICM")

coverage_df <- coverage_df |> dplyr::bind_rows(ICM_coverage_df)

coverage_plot <- ggplot2::ggplot(data=coverage_df,
                                 mapping = ggplot2::aes(x=test_X, y=coverage, color=method)) + 
  ggplot2::geom_point()

coverage_plot

coverage_df |> dplyr::group_by(method) |> dplyr::summarise(mean(coverage))


# length
length_df <- result |> dplyr::mutate(length = CI_upper - CI_lower) |> 
  dplyr::group_by(method, test_X) |> dplyr::summarise(length_mean = mean(length))

ICM_length_df <- ICM_res |> dplyr::mutate(length = upper - lower) |> 
  dplyr::group_by(test_X) |> dplyr::summarise(length_mean = mean(length)) |>
  dplyr::mutate(method = "ICM")

length_df <- length_df |> dplyr::bind_rows(ICM_length_df)
  
length_plot <- ggplot2::ggplot(data = length_df,
                               mapping = ggplot2::aes(x=test_X, y=length_mean, color=method)) + 
  ggplot2::geom_point()

length_plot

length_df |> dplyr::group_by(method) |> dplyr::summarise(mean(length_mean))


# RMSE
rmse_df <- result |> dplyr::mutate(squared_error = (mean-true_HTE)^2) |> 
  dplyr::group_by(method, test_X) |> dplyr::summarise(rmse = mean(squared_error))

ICM_rmse_df <- ICM_res |> dplyr::mutate(squared_error = (pred_HTE-true_HTE)^2) |> 
  dplyr::group_by(test_X) |> dplyr::summarise(rmse = mean(squared_error)) |>
  dplyr::mutate(method = "ICM")

kallus_rmse_df <- kallus_res |> dplyr::mutate(squared_error = (pred_HTE-true_HTE)^2) |> 
  dplyr::group_by(test_X) |> dplyr::summarise(rmse = mean(squared_error)) |>
  dplyr::mutate(method = "2step_ML")

rmse_df <- rmse_df |> dplyr::bind_rows(ICM_rmse_df)

rmse_plot <- ggplot2::ggplot(data = rmse_df,
                            mapping = ggplot2::aes(x=test_X, y=rmse, color=method)) + 
  ggplot2::geom_point()

rmse_plot
rmse_df |> dplyr::group_by(method) |> dplyr::summarise(mean(rmse))
kallus_rmse_df$rmse |> mean()

# bias plot
library(purrr) 
library(dplyr)
library(here)

num_repeat <- 100
bias_df_all <- map_df(1:num_repeat, \(i) {

  file_path <- here("03_analyze", "result", "0120", paste0("mcmc_sample_res__montecarlo_", i ,".obj"))
  mcmc_res <- readRDS(file_path)

  X_val <- mcmc_res$train_data$X[mcmc_res$train_data$ID == "O"]
  b_mat <- matrix(unlist(mcmc_res$b_O_plugin$mean), byrow = TRUE, ncol = 400)

  bb <- colMeans(b_mat)

  data.frame(
    iter = i,
    x = X_val,
    bias = bb
  )
}, .progress = TRUE) 

true_bias_df <- true_bias_df |> dplyr::mutate(bias = b)


bias_plot <- bias_df_all |> ggplot2::ggplot(mapping = ggplot2::aes(x, bias)) +
  ggplot2::geom_point(alpha = 0.1) + 
  ggplot2::geom_hline(yintercept = 0.55, color="blue", linewidth=1)
  #ggplot2::geom_point(data = true_bias_df, color = "blue", size = 0.7)
bias_plot

bias_mse_df <- bias_df |> dplyr::filter(-1 < x & x < 1) |>
  dplyr::group_by(number) |>
  dplyr::summarise(bias_mse = mean((b - 0.55)^2))


# one-shot
m <- 87
idx <- seq((m-1)*401+1, m*401)
proposal_df <- result |> dplyr::filter(method=="2step_proposal")
proposal_oneshot <- proposal_df[idx,] 

proposal_oneshot |> ggplot2::ggplot(ggplot2::aes(test_X)) +
  ggplot2::geom_line(ggplot2::aes(y=mean, color="pred"), linewidth=1) + 
  ggplot2::geom_line(ggplot2::aes(y=true_HTE, color="true"), linewidth=1) +
  ggplot2::scale_color_manual(name = "Legend", 
                              values = c("pred" = "blue", 
                                         "true" = "red")) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin=CI_lower, ymax=CI_upper),
                       color="gray", alpha=0.2)

bias_res <- readRDS(here::here("03_analyze", "result", "0113", paste0("b_O_plugin_montecarlo_",m,".obj")))

bias_res |> ggplot2::ggplot(mapping = ggplot2::aes(x, b)) +
  ggplot2::geom_point(alpha = 0.1) + 
  ggplot2::geom_hline(yintercept = 0.55, color="blue", linewidth=1) + 
  ggplot2::geom_hline(yintercept = 0, color="black", linewidth=0.5) + 
  ggplot2::ylim(c(-1.2, 1.2))

# tau-mse vs bias-mse(insample)
mse_df <- data.frame()
for(i in 1:num_repeat){
  idx <- seq((i-1)*401+1, i*401)
  res_df_i <- proposal_df[idx,] 
  mse_df_i <- data.frame(number=i, tau_mse = mean((res_df_i$mean- res_df_i$true_HTE)^2))
  mse_df <- mse_df |> dplyr::bind_rows(mse_df_i)
}

mse_mse_df <- mse_df |> dplyr::left_join(bias_mse_df, by="number")

mse_mse_df |>
  ggplot2::ggplot(mapping = ggplot2::aes(x=bias_mse, y=tau_mse)) + 
  ggplot2::geom_point()
cor(mse_mse_df$tau_mse, mse_mse_df$bias_mse)
