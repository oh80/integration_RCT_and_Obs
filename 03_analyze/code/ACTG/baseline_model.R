rm(list = ls())

data <- readRDS(here::here("02_build", "data", "0321", "ACTG_15.obj"))
train_df <- data$confound_train_df
test_df <- data$test_data

# pool
source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "predict.R"))
source(here::here("03_analyze", "code", "vanira_GP_MCMC.R"))

X <- train_df |> dplyr::select(c(cd40, cd80))
X_mean <- apply(X, 2, mean)
X_sd <- apply(X, 2, sd)
X_train_scaled <- scale(X, center = X_mean, scale = X_sd)

Y <- train_df$cd496
Y_mean <- mean(Y)
Y_sd <- sd(Y)
Y_train_scaled <- scale(Y, Y_mean, Y_sd)

Z <- train_df$treat

train_data <- list(X=X_train_scaled, Y=Y_train_scaled, Z=Z)
MCMC_res <- run_MCMC(train_data, iter=500, burn_in=200)

X_test <- test_df|> dplyr::select(c(cd40, cd80))
X_test_scaled <- scale(X_test, center = X_mean, scale = X_sd) |> as.matrix()
test_data <- list(X=X_test_scaled)
gp_pred_res <- compute_pred_and_CI(train_data, test_data, MCMC_res, "pool")
gp_pred_mean <- gp_pred_res$mean * Y_sd

saveRDS(gp_pred_mean, here::here("03_analyze", "result", "ACTG", "pool_pred"))


# RCT-only
train_data$ID <- train_df$ID
RCT_only_train_data <- extract_use_data(train_data, "RCT")
RCT_only_MCMC_res <- run_MCMC(RCT_only_train_data, iter=500, burn_in=200)

RCT_only_gp_pred_res <- compute_pred_and_CI(RCT_only_train_data, test_data, RCT_only_MCMC_res, "pool")
RCT_only_gp_pred_mean <- RCT_only_gp_pred_res$mean * Y_sd

saveRDS(RCT_only_gp_pred_mean, here::here("03_analyze", "result", "ACTG", "RCT_only_gp"))

# Obs-only
Obs_only_train_data <- extract_use_data(train_data, "observation")
Obs_only_MCMC_res <- run_MCMC(Obs_only_train_data, iter=500, burn_in=200)

Obs_only_gp_pred_res <- compute_pred_and_CI(Obs_only_train_data, test_data, Obs_only_MCMC_res, "pool")
Obs_only_gp_pred_mean <- Obs_only_gp_pred_res$mean * Y_sd

saveRDS(Obs_only_gp_pred_mean, here::here("03_analyze", "result", "ACTG", "Obs_only_gp"))
