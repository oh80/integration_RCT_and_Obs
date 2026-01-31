functions_path <- here::here("01_data", "code", "utils.R")
source(functions_path)

main <- function(){
  # settings
  seed = 42
  nO <- 400
  nR <- 150
  sig <- 1
  beta <- 2
  
  desctiption = ""
  
  # generate data
  data <- generate_data(seed, nO, nR, sig, beta)
  
  # add meta information
  info <- list(desctiption=desctiption, seed=seed, nO=nO, nR=nR,sig=sig, beta=beta)
  data$info = info
  
  
  # save
  data |> save_data("1d_logexp", sample_size=nO+nR)
}


logexp_CATE <- function(X, beta){
  cate <- log(1+exp(beta*X))
  return(cate)
}


generate_data <- function(seed, nO, nR, sig, beta){
  set.seed(seed)
  n <- nR + nO 
  ID <- c(rep("R", nR), rep("O", nO))
  
  # confounder
  X <- c()
  X[ID=="R"] <- runif(nR, -1, 1)
  X[ID=="O"] <- runif(nO, -2, 2)
  U <- runif(n, -1, 1)
  
  # assignment
  LGS <- function(x){ 1/(1+exp(-x)) }
  Pi <- c(rep(0.5, nR), LGS(2*U[ID=="O"]))
  Z <- rbinom(n, 1, Pi)
  
  # tau
  Tau <- logexp_CATE(X,beta)
  
  # outcome
  Base <- 1 + X + U
  Y <- Base + Z*Tau + sig*rnorm(n)
  
  data <- list("X"=X, "Y"=Y, "Tau"=Tau, "Base"=Base, "U"=U, "Z"=Z, "ID"=ID)
  return(data)
}


#main()