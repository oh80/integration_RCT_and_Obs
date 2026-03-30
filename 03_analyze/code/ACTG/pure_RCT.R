rm(list = ls())

data <- readRDS(here::here("02_build", "data", "0321", "ACTG_3.obj"))
train_df <- data$train_df
test_df <- data$test_data

# PL-GP
source(here::here("03_analyze", "code", "samplers.R"))
source(here::here("03_analyze", "code", "predict.R"))
source(here::here("03_analyze", "code", "vanira_GP_MCMC.R"))

X <- train_df |> dplyr::select(c(cd40, cd80, karnof))
X_mean <- apply(X, 2, mean)
X_sd <- apply(X, 2, sd)
X_train_scaled <- scale(X, center = X_mean, scale = X_sd)

Y <- train_df$cd496
Y_mean <- mean(Y)
Y_sd <- sd(Y)
Y_train_scaled <- scale(Y, Y_mean, Y_sd)

Z <- train_df$treat

train_data <- list(X=X_train_scaled, Y=Y_train_scaled, Z=Z)
MCMC_res <- run_MCMC(train_data, iter=50, burn_in=20)

X_test <- test_df|> dplyr::select(c(cd40, cd80, karnof))
X_test_scaled <- scale(X_test, center = X_mean, scale = X_sd) |> as.matrix()
test_data <- list(X=X_test_scaled)
gp_pred_res <- compute_pred_and_CI(train_data, test_data, MCMC_res, "pool")
gp_pred_mean <- gp_pred_res$mean * Y_sd
saveRDS(gp_pred_mean, here::here("03_analyze", "result", "ACTG", "gp_pred"))


# lm
lm_fit1 <- lm(cd496 ~ cd40  + cd80 + karnof, data = train_df[train_df$treat==1,])
lm_fit0 <- lm(cd496 ~ cd40  + cd80 + karnof, data = train_df[train_df$treat==0,])

lm_pred <- predict(lm_fit1, test_df|> dplyr::select(c(cd40, cd80, karnof))) - 
  predict(lm_fit0, test_df|> dplyr::select(c(cd40, cd80, karnof)))

saveRDS(lm_pred, here::here("03_analyze", "result", "ACTG", "lm_pred"))


# RF
# install.packages("randomForest")
set.seed(42)
rf_fit1 <- randomForest::randomForest(cd496 ~ cd40 + cd80 + karnof, 
                        data = train_df[train_df$treat == 1, ])

rf_fit0 <- randomForest::randomForest(cd496 ~ cd40 + cd80 + karnof, 
                        data = train_df[train_df$treat == 0, ])

rf_pred1 <- predict(rf_fit1, newdata = test_df |> dplyr::select(cd40, cd80, karnof))
rf_pred0 <- predict(rf_fit0, newdata = test_df |> dplyr::select(cd40, cd80, karnof))

rf_hte_pred <- rf_pred1 - rf_pred0
saveRDS(rf_hte_pred, here::here("03_analyze", "result", "ACTG", "RF_pred"))
