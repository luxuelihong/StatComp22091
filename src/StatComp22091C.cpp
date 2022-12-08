#include <Rcpp.h>
using namespace Rcpp;

//' @title generate gibbs using Rcpp
//' @description generate gibbs using Rcpp
//' @param N the number of samples
//' @param mu1 Mean value of the first parameter
//' @param mu2 Mean value of the second parameter
//' @param sigma1 Standard deviation of the first parameter
//' @param sigma2 Standard deviation of the second parameter
//' @param rho correlation coefficient of the two parameter
//' @return Bivariate normal distribution x
//' @export
//[[Rcpp::export]]
NumericMatrix gibbs_cpp(int N,double mu1,double mu2,double sigma1,double sigma2,double rho){
  NumericMatrix X(N,2);
  double s1,s2;
  s1=sqrt(1-pow(rho,2))*sigma1;
  s2=sqrt(1-pow(rho,2))*sigma2;
  X(0,0)=mu1;
  X(0,1)=mu2;
  double x1,x2,m1,m2;
  for(int i=1;i<N;i++){
    x2=X(i-1,1);
    m1=mu1+rho*(x2-mu2)*sigma1/sigma2;
    X(i,0)=rnorm(1,m1,s1)[0];
    x1=X(i,0);
    m2=mu2+rho*(x1-mu1)*sigma2/sigma1;
    X(i,1)=rnorm(1,m2,s2)[0];
  }
  return(X);
} 

