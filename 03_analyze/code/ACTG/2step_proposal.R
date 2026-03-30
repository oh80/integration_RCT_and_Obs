rm(list = ls())

data <- readRDS(here::here("02_build", "data", "0321", "ACTG_15.obj"))
train_df <- data$confound_train_df
test_df <- data$test_data

source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "predict.R"))
source(here::here("03_analyze", "code", "2step_GP_MCMC.R"))

seed    <- 42
iter    <- 500
burn_in <- 200

X <- train_df |> dplyr::select(c(cd40, cd80))
X_mean <- apply(X, 2, mean)
X_sd <- apply(X, 2, sd)
X_train_scaled <- scale(X, center = X_mean, scale = X_sd)

Y <- train_df$cd496
Y_mean <- mean(Y)
Y_sd <- sd(Y)
Y_train_scaled <- scale(Y, Y_mean, Y_sd)


Z  <- train_df$treat |> as.vector()
ID <- train_df$ID |> as.vector()

min_max_list <- get_RCT_sapport(X_train_scaled, ID)
insample_flag <- extract_insample_id(X_train_scaled, min_max_list)

X_in  <- X_train_scaled[insample_flag, , drop=FALSE]
Y_in  <- Y_train_scaled[insample_flag]
Z_in  <- Z[insample_flag]
ID_in <- ID[insample_flag]

# step-1 (insample GP)
set.seed(seed)
step_1_samples <- run_MCMC(X_in, Y_in, Z_in, ID_in, iter, burn_in, step = 1)

# step-2 (full-sample GP)
b_O_plugin <- compute_b_O_pred(step_1_samples, X_in, X_train_scaled, ID_in, ID)
step_2_samples <- run_MCMC(X_train_scaled, Y_train_scaled, Z, ID, iter, burn_in, step = 2, b_O_plugin)

# predict
X_test <- test_df|> dplyr::select(c(cd40, cd80))
X_test_scaled <- scale(X_test, center = X_mean, scale = X_sd) |> as.matrix()
pred_res <- compute_pred_and_CI(list(X=X_train_scaled), list(X = X_test_scaled), step_2_samples, method="proposal")

proposal_pred_mean <- pred_res$mean * Y_sd
saveRDS(proposal_pred_mean, here::here("03_analyze", "result", "ACTG", "proposal_2step_pred_v2"))
