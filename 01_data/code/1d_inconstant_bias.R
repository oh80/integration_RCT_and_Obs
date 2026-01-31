functions_path <- here::here("01_data", "code", "utils.R")
source(functions_path)

main <- function(){
  # settings
  seed = 1
  nO <- 200
  nR <- 50
  
  sig   <- 1
  alpha <- 0.5
  
  desctiption = ""
  
  # generate data
  data <- generate_data(seed, nO, nR, sig, alpha)
  
  # add meta information
  info <- list(desctiption=desctiption, seed=seed, nO=nO, nR=nR,sig=sig)
  data$info = info
  
  
  # save
  data |> save_data("1d_inconstant_bias", sample_size=nO+nR)
}


squeared_CATE <- function(X){
  return(X^2)
}


generate_data <- function(seed, nO, nR, sig, alpha){
  set.seed(seed)
  n <- nR + nO 
  ID <- c(rep("R", nR), rep("O", nO))
  
  # confounder
  X <- c()
  X[ID=="R"] <- runif(nR, -1, 1)
  X[ID=="O"] <- runif(nO, -2, 2)
  
  # unobserved confounder
  d <- 1+alpha*exp(-X^2)
  U <- c()
  for(i in 1:n){
    U[i] <- runif(1, min=-d[i], max=d[i])
  }
  
  # assignment
  LGS <- function(x){ 1/(1+exp(-x)) }
  Pi <- c(rep(0.5, nR), LGS(2*U[ID=="O"]))
  Z <- rbinom(n, 1, Pi)
  
  # tau
  Tau <- squeared_CATE(X)
  
  # outcome
  Base <- 1 + X + U
  Y <- Base + Z*Tau + sig*rnorm(n)
  
  data <- list("X"=X, "Y"=Y, "Tau"=Tau, "Base"=Base, "U"=U, "Z"=Z, "ID"=ID)
  return(data)
}


#main()