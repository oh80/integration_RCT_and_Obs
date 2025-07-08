get_rmse <- function(y_true, y_pred){
  mse <- mean((y_true - y_pred)^2)
  rmse <- mse |> sqrt()
  return(rmse)
}