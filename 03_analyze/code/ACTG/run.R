# source(here::here("03_analyze", "code","ACTG", "pure_RCT.R"))
# source(here::here("03_analyze", "code","ACTG", "baseline_model.R"))
# source(here::here("03_analyze", "code","ACTG", "3step_proposal.R"))


gp_res <- readRDS(here::here("03_analyze", "result", "ACTG", "gp_pred"))
rf_res <- readRDS(here::here("03_analyze", "result", "ACTG", "RF_pred"))
lm_res <- readRDS(here::here("03_analyze", "result", "ACTG", "lm_pred"))
pool_res <- readRDS(here::here("03_analyze", "result", "ACTG", "pool_pred"))
RCT_only_res <- readRDS(here::here("03_analyze", "result", "ACTG", "RCT_only_gp"))
Obs_only_res <- readRDS(here::here("03_analyze", "result", "ACTG", "Obs_only_gp"))
prop3_res <- readRDS(here::here("03_analyze", "result", "ACTG", "proposal_3step_pred_v2"))
prop2_res <- readRDS(here::here("03_analyze", "result", "ACTG", "proposal_2step_pred_v2"))
ICM_res <- read.csv(here::here("03_analyze", "result", "ACTG", "ICM_pred.csv"), header = FALSE)
ICM_res <- ICM_res$V1 |> as.vector()
kallus_res <- read.csv(here::here("03_analyze", "result", "ACTG", "Kallus_pred.csv"), header = FALSE)
kallus_res <- kallus_res$V1 |> as.vector()


get_rmse <- function(y_true, y_pred){
  mse <- mean((y_true - y_pred)^2)
  rmse <- mse |> sqrt()
  return(rmse)
}

pure_RCT_list <- list(gp_res, rf_res, lm_res)
comp_list <- list(pool_res, RCT_only_res, Obs_only_res, prop2_res, prop3_res, ICM_res, kallus_res)

rmse_df <- data.frame(gp = rep(0, 7), rf = rep(0, 7), lm = rep(0, 7))
for(i in 1:3){
  for(j in 1:7){
    rmse <- get_rmse(pure_RCT_list[[i]], comp_list[[j]])
    rmse_df[j,i] <- rmse
  }
}
row.names(rmse_df) <- c("pool", "RCT_only", "Obs_only", "2step_proposal", "3step_proposal", "ICM_res", "2step-ML")

cor(gp_res, rf_res)

cor_df <- data.frame(gp = rep(0, 7), rf = rep(0, 7), lm = rep(0, 7))
for(i in 1:3){
  for(j in 1:7){
    cor <- cor(pure_RCT_list[[i]], comp_list[[j]])
    cor_df[j,i] <- cor
  }
}
row.names(cor_df) <- c("pool", "RCT_only", "Obs_only", "2step_proposal", "3step_proposal", "ICM_res", "2step-ML")


ATE_df <- data.frame(ATE = rep(0, 7))
for(j in 1:7){
    ATE <- mean(comp_list[[j]])
    ATE_df$ATE[j] <- ATE
  }

row.names(ATE_df) <- c("pool", "RCT_only", "Obs_only", "2step_proposal", "3step_proposal", "ICM_res", "2step-ML")


data <- readRDS(here::here("02_build", "data", "0321", "ACTG_15.obj"))
test_df <- data$test_data

ATE_true <- mean(test_df$cd496[test_df$treat==1]) - mean(test_df$cd496[test_df$treat==0])

# scatter
plot_df <- data.frame("pure_RCT" = gp_res, "pool_res" = pool_res,
                      "2step" = prop2_res, "3step" = prop3_res,
                      "ICM_res" = ICM_res, check.names = FALSE) |> 
  tidyr::pivot_longer(
    cols = c("pool_res", "ICM_res", "2step", "3step"),
    names_to = "model",
    values_to = "predicted"
  ) |>
  dplyr::mutate(model = factor(model, 
                        levels = c("pool_res", "ICM_res", "2step", "3step"),
                        labels = c("Pool Model", "ICM", "2-Step (Proposed)", "3-Step (Proposed)")))

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pure_RCT, y = predicted)) +
  ggplot2::geom_point(alpha = 0.5, color = "steelblue", size = 2) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  ggplot2::facet_wrap(~ model, nrow = 2, ncol = 2) +
  ggplot2::labs(
    x = "pure RCT (GP)",
    y = "Predicted HTE"
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "#f0f0f0", color = "black"),
    strip.text = ggplot2::element_text(face = "bold", size = 12),
    axis.title.x = ggplot2::element_text(face = "bold", margin = ggplot2::margin(t = 10)),
    axis.title.y = ggplot2::element_text(face = "bold", margin = ggplot2::margin(r = 10)),
    panel.grid.minor = ggplot2::element_blank()
  )
p
