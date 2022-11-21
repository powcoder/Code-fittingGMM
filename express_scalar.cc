https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <iostream>

using namespace Rcpp ;


double find_max(arma::vec& z, int J){
  double t = z[0];
  int i;
  for (i = 1; i < J; i++)
  {
    if (z[i] > t)
      t = z[i];
  }
  return t;
}

double find_max_arma(arma::rowvec& z, int J){
  double t = z[0];
  int i;
  for (i = 1; i < J; i++)
  {
    if (z[i] > t)
      t = z[i];
  }
  return t;
}


// [[Rcpp::export(".loopexp_scalar")]]
Rcpp::NumericVector loopexp_out(const Rcpp::NumericVector& helpa, const Rcpp::NumericVector& helpb, const arma::vec& helpc,
                                const Rcpp::NumericVector& X, int J, int n)
{
  int i, j;
  Rcpp::NumericVector Z(n);
  RNGScope scope;
  for (i = 1; i <= n; i++)
  {
    arma::vec prob(J);
    
    for (j = 1; j <= J; j++)
    {
      double dnor = Rf_dnorm4(X[i-1], helpb[j-1], helpc[j-1], 1);
      prob[j-1] = helpa[j-1] + dnor;
    }
    
    
    double mmax = find_max(prob, J);
    double deno = mmax + log(sum(exp(prob - mmax)));
    Rcpp::IntegerVector JJ = seq_len(J) ; 
    
    Z[i-1] = RcppArmadillo::sample(JJ, 1, false, (exp(prob - deno)))[0];
  }
  return Z;
}

// [[Rcpp::export(".loopexp_scalar_prob")]]
arma::mat loopexp_out2(const Rcpp::NumericVector& helpa, const Rcpp::NumericVector& helpb, 
                       const arma::vec& helpc, const Rcpp::NumericVector& X, int J, int n)
{
  int i, j;
  arma::mat prob(n,J);
  RNGScope scope;
  for (i = 1; i <= n; i++)
  {
    
    for (j = 1; j <= J; j++)
    {
      double dnor = Rf_dnorm4(X[i-1], helpb[j-1], helpc[j-1], 1);
      prob(i-1, j-1) = helpa[j-1] + dnor;
    }
    
    arma::rowvec test = prob.row(i-1);
    double mmax = find_max_arma(test, J);
    double deno = mmax + log(sum(exp(prob.row(i-1) - mmax)));
    prob.row(i-1) = exp(prob.row(i-1) - deno);
  }
  return prob;
}

// [[Rcpp::export(".loopexp_scalar_prob_n")]]
NumericVector loopexp_out3(const Rcpp::NumericVector& helpa, const Rcpp::NumericVector& helpb, 
                           const arma::vec& helpc, const Rcpp::NumericVector& X, int J, int n,
                           const arma::mat& alpmcmc)
{
  NumericVector pq(J*J);
  int ct = 0;
  int k, j;
  arma::mat alpbem= loopexp_out2(helpa, helpb, helpc, X, J, n);
  for(k = 0; k < J; k++){
    for(j = 0; j < J; j++){
      pq[ct] = sum(abs((alpbem.col(k) - alpmcmc.col(j))));
      ct +=1;
    }
  }
  return pq;
}

