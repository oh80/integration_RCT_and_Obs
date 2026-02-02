N <- 500
X <- runif(n = N, min = -2, max = 2) |> as.matrix()
y <- X^2 + rnorm(n=N, sd = 0.5)
plot(X, y)


compute_kernel_mat <- function(X, l, eta){
  n <- nrow(X)
  base_mat <- matrix(0, nrow=n, ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      if(i > j){
        base <- -1/2 * sum((X[i,] - X[j,])^2)
        base_mat[i, j] <- base
        base_mat[j, i] <- base
      }}}
  log_K <- base_mat / l^2
  K <- eta^2 * exp(log_K)
  diag(K) <- eta^2
  
  return(K)
}

compute_train_test_kernel <- function(X, X2, l, eta) {
  n1 <- nrow(X)
  n2 <- nrow(X2)
  
  K <- matrix(0, nrow = n1, ncol = n2)
  
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      sqdist <- sum((X[i, ] - X2[j, ])^2)
      exponent <- -1/2 * sqdist / (l^2)
      K[i, j] <- eta^2 * exp(exponent)
    }
  }
  return(K)
}

# prepare NNGP
compute_distance <- function(x1, x2){
  distance <- (x1 - x2)^2 |> sum()
  return(distance)
}


get_neighbor_sample <- function(x, X_candidate, m){
  if(nrow(X_candidate) <= m){
    return(seq(1,nrow(X_candidate)))
    
  }else{
    distance_vec <- apply(X_candidate, 1, compute_distance, x1 = x)
    neighbor_idx <- order(distance_vec)[1:m]
    
    return(neighbor_idx)
  }
}


get_neighbor_sample_list <- function(X, m){
  N <- dim(X)[1]
  
  output_list      <- list()
  output_list[[1]] <- 1
  
  output_list[2:N] <- purrr::map(
    .x = seq(2, N),
     ~ get_neighbor_sample(x = X[.x, ] |> as.matrix(),
                             X_candidate = X[1:(.x-1), ] |> as.matrix(),
                             m = m))
  return(output_list)
}


get_children <- function(i, idx_list, N){
  output <- c()
  cnt <- 1
  for(t in 1:N){
    idx_vec <- idx_list[[t]]
    if(i %in% idx_vec){
      output[cnt] <- t
      cnt <- cnt + 1
    }
  }
  return(output)
}


get_children_list <- function(neighbor_sample_list){
  N <- length(neighbor_sample_list)
  output_list <- list()
  
  output_list[1:N] <- purrr::map(
    .x = seq(1, N),
    ~ get_children(i = .x,
                   idx_list = neighbor_sample_list,
                   N = N))
  return(output_list)
}


neighbor_sample_list <- get_neighbor_sample_list(X, m=30)
children_list <- get_children_list(neighbor_sample_list)


# MCMC
l <- 1
eta <- 0.5
epsilon <- 2

B_mat_list <- list()
F_mat_list <- list()

for(i in 1:N){
  if(i == 1){
    B_mat_list[[1]] <- 1
    F_mat_list[[1]] <- 1
  } else {
    
    idx <- neighbor_sample_list[[i]]
    
    X_i <- X[i, , drop=FALSE]
    X_n <- X[idx, , drop=FALSE]

    K_n   <- compute_kernel_mat(X_n, l, eta) + 1e-6 * diag(length(idx)) 
    K_i_n <- compute_train_test_kernel(X_i, X_n, l, eta)
    K_i   <- compute_kernel_mat(X_i, l, eta)
    
    B_mat <- K_i_n %*% solve(K_n)

    F_mat <- K_i - B_mat %*% t(K_i_n)
    
    B_mat_list[[i]] <- B_mat
    F_mat_list[[i]] <- F_mat
  }
}

gp <- y

for(i in 1:N){
  n_idx <- neighbor_sample_list[[i]]
  c_idx <- children_list[[i]]
  
  # posterior variance
  V_inv <- 1/epsilon + 1/as.numeric(F_mat_list[[i]])
  for(t in c_idx){
    p_of_t <- neighbor_sample_list[[t]]
    pos <- which(p_of_t == i)
    b_ti <- as.numeric(B_mat_list[[t]][pos])
    V_inv <- V_inv + (b_ti^2) / as.numeric(F_mat_list[[t]])
  }
  V <- 1 / V_inv
  
  # posterior mean
  p_pred <- 0
  if(i > 1) p_pred <- sum(as.numeric(B_mat_list[[i]]) * gp[n_idx])
  
  mu <- y[i]/epsilon + p_pred / as.numeric(F_mat_list[[i]])
  
  # information from children
  for(t in c_idx){
    p_of_t <- neighbor_sample_list[[t]]
    pos <- which(p_of_t == i)
    b_ti <- as.numeric(B_mat_list[[t]][pos])
    
    o_idx <- p_of_t[-pos]
    b_o   <- B_mat_list[[t]][-pos]
    
    # residuals from children
    resid_t <- gp[t] - sum(as.numeric(b_o) * gp[o_idx])
    mu <- mu + (b_ti * resid_t) / as.numeric(F_mat_list[[t]])
  }
  
  gp[i] <- rnorm(1, mean = V * mu, sd = sqrt(V))
}

plot(X, gp)


# predict
test_X <- seq(-2,2,0.01)

test_X_neighbor_list <- list()
test_X_neighbor_list[1:length(test_X)] <- purrr::map(
  .x = seq(1, length(test_X)),
  ~ get_neighbor_sample(x = test_X[.x] |> as.matrix(),
                        X_candidate = X,
                        m = 30))

N_sample <- 200
gp_test <- matrix(NA, nrow=N_sample, ncol=length(test_X))

for(i in 1:length(test_X)){
    idx <- test_X_neighbor_list[[i]]
    
    X_i <- test_X[i] |> as.matrix()
    X_n <- X[idx, , drop=FALSE]
    
    K_n   <- compute_kernel_mat(X_n, l, eta) + 1e-6 * diag(length(idx)) 
    K_i_n <- compute_train_test_kernel(X_i, X_n, l, eta)
    K_i   <- compute_kernel_mat(X_i, l, eta)
    
    B_mat <- K_i_n %*% solve(K_n)
    
    F_mat <- K_i - B_mat %*% t(K_i_n)
    
    gp_test[,i] <- rnorm(n=N_sample,
                         mean = sum(B_mat * gp[idx]),
                         sd = sqrt(F_mat))
}


gp_test_mean <- apply(gp_test, 2, mean)
gp_test_upper <- apply(gp_test, 2, quantile, probs=0.975)
gp_test_lower <- apply(gp_test, 2, quantile, probs=0.025)

plot(X, y, col = "gray", pch = 16, main = "NNGP Fit and Predict")
lines(test_X, gp_test_mean, col = "blue", lwd = 2)
lines(test_X, gp_test_upper, col = "red", lty = 2)
lines(test_X, gp_test_lower, col = "red", lty = 2)

gc()