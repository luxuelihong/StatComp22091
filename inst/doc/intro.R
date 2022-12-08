## ----eval=FALSE---------------------------------------------------------------
#  delta_i = function(data_i, big_delta_symbolic, params_first,  params_second, param_vals
#  ){
#  
#    # Iterate through the symbolic expression for Delta and
#    # evaluate each element numerically.
#    big_delta = matrix(NA, nrow = length(params_second), ncol = length(params_first))
#    env = as.list(c(data_i, param_vals))
#    for (r in 1:nrow(big_delta)) {
#      for (c in 1:ncol(big_delta)) {
#        big_delta[r,c] = eval(
#          big_delta_symbolic[r,c][[1]],
#          env = env)
#      }
#    }
#  
#    return(big_delta)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix gibbs_cpp(int N,double mu1,double mu2,double sigma1,double sigma2,double rho){
#    NumericMatrix X(N,2);
#    double s1,s2;
#    s1=sqrt(1-pow(rho,2))*sigma1;
#    s2=sqrt(1-pow(rho,2))*sigma2;
#    X(0,0)=mu1;
#    X(0,1)=mu2;
#    double x1,x2,m1,m2;
#    for(int i=1;i<N;i++){
#      x2=X(i-1,1);
#      m1=mu1+rho*(x2-mu2)*sigma1/sigma2;
#      X(i,0)=rnorm(1,m1,s1)[0];
#      x1=X(i,0);
#      m2=mu2+rho*(x1-mu1)*sigma2/sigma1;
#      X(i,1)=rnorm(1,m2,s2)[0];
#    }
#    return(X);
#  }

