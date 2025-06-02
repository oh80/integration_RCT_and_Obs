functions_path <- here::here("01_data", "code", "utils.R")
source(functions_path)

main <- function(){
  # settings
  seed = 42
  nO <- 100
  nR <- 50
  sig <- 1
  true_CATE = "squared"
  
  desctiption = ""
  
  # generate data
  data <- generate_data(seed, nO, nR, sig, true_CATE)
  
  # add meta information
  info <- list(desctiption=desctiption, seed=seed, nO=nO, nR=nR,sig=sig, true_CATE=true_CATE)
  data$info = info
  # data <- data |> add_meta_info(desctiption=desctiption,
  #                               seed=seed, nO=nO, nR=nR,
  #                               sig=sig, true_CATE=true_CATE)
  
  # save
  data |> save_data("1d", sample_size=nO+nR)
}


generate_data <- function(seed, nO, nR, sig, true_CATE){
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
  if(true_CATE=="squared"){
    CATE <- function(x){ x^2 }
  }
  Tau <- CATE(X)
  
  # outcome
  Base <- 1 + X + U
  Y <- Base + Z*Tau + sig*rnorm(n)
  
  data <- list("X"=X, "Y"=Y, "Tau"=Tau, "Base"=Base, "U"=U, "Z"=Z, "ID"=ID)
  return(data)
}



main()