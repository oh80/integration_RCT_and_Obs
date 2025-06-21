functions_path <- here::here("01_data", "code", "utils.R")
source(functions_path)

main <- function(){
  # settings
  seed = 42
  n <- 500
  
  desctiption = ""
  
  # generate data
  data <- generate_data(seed, n)
  
  # add meta information
  info <- list(desctiption=desctiption, seed=seed, n=n)
  data$info = info
  
  # save
  data |> save_data("1d_Dimitriou", sample_size=n)
}


participant_prob <- function(x){
  p_s <- exp(-3-3*x)/(1+exp(-3-3*x))
  return(p_s)
}

decision_participant <- function(X){
  n <- nrow(X)
  P_s <- apply(X, 1, participant_prob)
  S <- rbinom(n, size=1, prob = P_s)
  
  ID <- rep("O", n)
  ID[S==1] <- rep("R", sum(S))
  return(ID)
}

decision_assingment <- function(X, ID){
  n <- nrow(X)
  Z <- rep(0, n)
  for(i in 1:n){
    # Observation
    if(ID[i]=="O"){
      e_x <- 1/(1+exp(X[i]))
      Z[i] <- rbinom(n=1, size=1, prob=e_x)
      
      # RCT 
    }else if(ID[i]=="R"){
      Z[i] <- rbinom(n=1, size=1, prob=0.5)}}
  return(Z)
}

true_CATE <- function(x){
  CATE <- 1+X+X^2
  return(CATE)
}

geneate_outcome <- function(X, ID, Z, Tau){
  n <- nrow(X)
  Y <- rep(0, n)
  for(i in 1:n){
    epsilon <- rnorm(n=1, mean=0, sd=1)
    # Observation
    if(ID[i]=="O"){
      U <- rnorm(n=1, mean=(2*Z-1)*sin(X-1), sd=1)
      Y[i] <- Z[i] * Tau[i] * X[i]^2 - 1 + U + epsilon
      
      # RCT   
    }else if(ID[i]=="R"){
      Y[i] <- Z[i] * Tau[i] * X[i]^2 - 1 + epsilon
    }}
  return(Y)
}


generate_data <- function(seed, n){
  set.seed(42)
  
  # confounder
  X <- runif(n=n, min=-2, max=2) |> as.matrix()
  
  # experiment participant
  ID <- decision_participant(X)
  
  # assingment
  Z <- decision_assingment(X, ID)
  
  # CATE and outcome
  Tau <- apply(X, 2, true_CATE)
  Y <- geneate_outcome(X, ID, Z, Tau)
  
  data <- list("X"=X, "Y"=Y, "Tau"=Tau, "Z"=Z, "ID"=ID)
  return(data)
}


main()