compute_E <- function(iter = 10000, Z){
  u_sample <- c()
  Max <- 1/(1+exp(-2))

  
  for(i in 1:iter){
    ui <- runif(n=1, min=-1, max = 1)

    if(Z==1){
      zi <- 1/(1+exp(-2*ui))
    }else if(Z==0){
      zi <- 1/(1+exp(2*ui))
    }

    ri <- zi/Max
    
    v <- runif(n=1)
    if(v < ri){
      u_sample <- append(u_sample, ui)
    }
  }
  
  E <- mean(u_sample)
  return(E)
}


e1 <- compute_E(Z=1)
e0 <- compute_E(Z=0)

e1 - e0
