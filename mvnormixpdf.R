mvnormixpdf <- function(X, pi,mu,Sigma){
  # X  = n.p array of n values of p-dim MV normal mixture 
  # w = column k vector probs
  # mu = p.k matrix of component means 
  # Sigma = p.p.k array of variance matrice 
  # pdf = n vector of pdf values
  #
  p <- nrow(mu); k <- ncol(mu); n <- nrow(X); pdf <- rep(0, n)
  for (j in 1:k){
    pdf <- pdf + pi[j]*dmvnorm(X, mu[,j], Sigma[,,j])
  }
  return (pdf)
}




