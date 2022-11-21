# set path under the folder Ethan Lynn Qiang/Code_fittingGMM"
#setwd("/Users/ethanfangxy/Dropbox/Ethan_Lynn_Qiang/Code_fittingGMM")
library(MASS)
library(MCMCpack)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
sourceCpp("Code_fittingGMM/express.cc", rebuild = TRUE, verbose = FALSE)
sourceCpp("Code_fittingGMM/express_scalar.cc", rebuild = TRUE, verbose = FALSE)
source("Code_fittingGMM/mvnormixpdf.R")
source("Code_fittingGMM/mvnormixrnd.R")
source("Code_fittingGMM/modesearch.R")
source("Code_fittingGMM/modesearchmerge.R")
### C code for sample Z (latent indicator)
set.seed(12345)

p <- ncol(Xm)
n <- nrow(Xm)
X = Xm
#signature(X="matrix")
X<-as.big.matrix(X)


#### Initialization ####
lambda <- 5
delta <- p-1+5
Phi <- diag(10,p)
m <- rep(0, p)
J = 120; # upper bound for the number of normal components
e = 20; f = 1; #alpha

TT = 3000
#################################
######Initialize parameters for MCMC
alpha <- rep(0, TT); alpha[1] <- 50
nu <- array(0, dim = c(J-1, TT))
nu[,1] <- rbeta(J-1,1,alpha[1])
pp <- array(0, dim = c(J, TT))
pp[1,1] <- nu[1,1]
tmp <- 1-nu[1,1]
for(jj in 2:(J-1)){
  pp[jj,1] = tmp * nu[jj,1]
  tmp <- tmp * (1-nu[jj,1])
}
pp[J,1] <- tmp
mu <- array(0, dim = c(p,J,TT))
Sigma <- vector("list", TT) 
Sigma[[1]]<- array(0, dim=c(p,p,J))
for(j in 1:J){
  Sigma[[1]][,,j] <- riwish(delta, Phi)
  mu[, j, 1] <- mvrnorm(1, m, lambda*Sigma[[1]][,,j])
}

source("Code_fittingGMM/Gibbs_DP.R")

?save.image("DPGMM_simu1_Gibbs.RData")
LpEM = -6.1110e+005; 
for(rr in 1:10){ # Number of EM algorithms
  inter <- TT-(rr-1)*1000 # Change the starting point for EM algorithms by using different MCMC samples. Make sure inter > 0.
  W_Gibbs <- pp[,inter]
  M_Gibbs <- array(unlist(mu[,,inter]),dim = c(p,J))
  Sig_Gibbs <- array(unlist(Sigma[[inter]]), dim = c(p, p, J))
  alpha_Gibbs <- alpha[inter]
  source("Code_fittingGMM/EM_DP.R")
  if(Lp > LpEM){
    W_final <- W # proportion
    M_final <- M # mean
    V_final <-V
    VV_final <- VV
    Sigma_final <- Sig # Covariance matrix
    alpha_final <- alpha_EM
    LpEM <- Lp
  }
}
save.image("DPGMM_simu1_EM.RData")
## Mode search 
tol <- 0.0001
modes <- modesearch(p, W_final, M_final, Sigma_final,0,tol)
m <- modes$m
pm <- modes$pm
imodes <- modes$imodes

tmp <- modesearchmerge(m, pm, imodes,W_final,tol)

save.image("DPGMM_simu1.RData")