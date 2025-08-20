
result <- readRDS(here::here("03_analyze", "result", "0820", "1d_squared_montecarlo_1.obj"))

# coverage
coverage_df <- result |>
  dplyr::mutate(isin_intervel = ifelse(CI_lower <= true_HTE & true_HTE <= CI_upper, 1, 0)) |>
  dplyr::group_by(method, test_X) |> dplyr::summarise(coverage = mean(isin_intervel))


coverage_plot <- ggplot2::ggplot(data=coverage_df,
                                 mapping = ggplot2::aes(x=test_X, y=coverage, color=method)) + 
  ggplot2::geom_point()

coverage_plot


# length
length_df <- result |> dplyr::mutate(length = CI_upper - CI_lower) |> 
  dplyr::group_by(method, test_X) |> dplyr::summarise(length_mean = mean(length))
  
length_plot <- ggplot2::ggplot(data = length_df,
                               mapping = ggplot2::aes(x=test_X, y=length_mean, color=method)) + 
  ggplot2::geom_point()

length_plot


# RMSE
rmse_df <- result |> dplyr::mutate(squared_error = (mean-true_HTE)^2) |> 
  dplyr::group_by(method, test_X) |> dplyr::summarise(rmse = mean(squared_error))

rmse_plot <- ggplot2::ggplot(data = rmse_df,
                            mapping = ggplot2::aes(x=test_X, y=rmse, color=method)) + 
  ggplot2::geom_point()

rmse_plot


# CATE (one-shot)
num_test <- 401
sim_number <- 2
one_shot_df <- result[(2*num_test*(sim_number-1)+1 ): (2*num_test*sim_number),]

prop_df <- one_shot_df |> dplyr::filter(method=="proposal")
pool_df <- one_shot_df |> dplyr::filter(method=="pool")

prop_df_for_plot <- data.frame("X" = rep(prop_df$test_X, 4),
                               "CATE" = c(prop_df$CI_lower, prop_df$CI_upper, prop_df$mean, prop_df$true_HTE),
                               "label" = c(rep("lower", num_test),rep("upper", num_test),rep("mean", num_test),rep("true", num_test)))

pool_df_for_plot <- data.frame("X" = rep(pool_df$test_X, 4),
                               "CATE" = c(pool_df$CI_lower, pool_df$CI_upper, pool_df$mean, pool_df$true_HTE),
                               "label" = c(rep("lower", num_test),rep("upper", num_test),rep("mean", num_test),rep("true", num_test)))

prop_cate_plot <- ggplot2::ggplot(data = prop_df_for_plot,
                             mapping = ggplot2::aes(x=X, y=CATE, color=label)) + 
  ggplot2::geom_point()

pool_cate_plot <- ggplot2::ggplot(data = pool_df_for_plot,
                                  mapping = ggplot2::aes(x=X, y=CATE, color=label)) + 
  ggplot2::geom_point()


prop_cate_plot
pool_cate_plot
