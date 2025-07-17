#### read data ####
data_path <- here::here("01_data", "data", "0708", "1d_Dimitriou_n500_3.obj")
data <- readRDS(data_path)

res_path <- here::here("03_analyze","result", "0708","proposal_1d_Dimitriou_n500_3_1.obj")
res <- readRDS(res_path)

samples <- res$samples

test_X  <- seq(-2, 2, by=0.01)|> as.matrix()
train_X <- data$X|> as.matrix()
#train_X <- data$X[data$ID=="R"]|> as.matrix()

#true_HTE <- test_X^2
true_HTE <- 1+test_X+test_X^2

#### predict ####
source(here::here("04_report", "code", "predict.R"))
pred_dist_samples <- get_pred_dist_samples(samples, train_X, test_X)

pred_mean <- apply(pred_dist_samples$mean, 2, mean)
pred_var  <- apply(pred_dist_samples$var, 2, mean)

calc_IC_bound <- function(mean, var){
  if(var < 0){
    var <- 1e-5
  }
  upper <- qnorm(p=0.975, mean=mean, sd=sqrt(var))
  lower <- qnorm(p=0.025, mean=mean, sd=sqrt(var))
  
  return(c(upper, lower))
}

# calc_CI_boound <- function(vec){
#   return(quantile(vec, probs=c(0.025, 0.975)))
# }
# pred_CI_bound <- apply(pred_samples, 2, calc_CI_boound)

CI_upper <- c()
CI_lower <- c()

for(i in 1:length(pred_mean)){
  CI_bound <- calc_IC_bound(pred_mean[i], pred_var[i])
  CI_upper[i] <- CI_bound[1]
  CI_lower[i] <- CI_bound[2]
}


#### plot ####
num_test <- length(true_HTE)

df_for_plot <- data.frame("X"=rep(test_X, 4),
                          "HTE"=c(CI_lower, CI_upper, pred_mean, true_HTE),
                          "label"=c(rep("lower",num_test),rep("upper",num_test),rep("mean",num_test),rep("true",num_test)))
plot <- ggplot2::ggplot(data=df_for_plot,
                        mapping = ggplot2::aes(x=X, y=HTE, color=label)) +
  ggplot2::geom_point()
plot

plot(test_X, pred_mean)

plot(train_X, apply(samples$tau,2,mean))
hist(train_X)



b_samples <- samples$b_O
b_train <- apply(b_samples, 2, mean)
train_obs_X <- data$X[data$ID=="R"]
plot(train_obs_X, b_train)
hist(train_obs_X)
