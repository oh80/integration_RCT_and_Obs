rm(list=ls())

# settings
x_set <- seq(-2, 2, by=0.01)
L <- length(x_set)

S <- 10000
alpha <- 0.5

# Monte Carlo 
Bias <- rep(NA, L)
U_e <- runif(S, -1, 1)

for(l in 1:L){
  x <- x_set[l]
  d <- 1+alpha*exp(-x^2)
  U <- d*U_e
  
  # Z=1
  w1 <- 1/(1+exp(-2*U))
  w1 <- w1/sum(w1)
  
  # Z=0
  w0 <- 1/(1+exp(2*U))
  w0 <- w0/sum(w0)
  
  # bias 
  Bias[l] <- sum(w1*U) - sum(w0*U)
}


# plot 
plot(x_set, Bias, type="l", xlab="x")
true_bias_df <- data.frame(x = x_set, b = Bias)
