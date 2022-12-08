## -----------------------------------------------------------------------------
plot(pressure)

## -----------------------------------------------------------------------------
n<-1000
u<-runif(n)
x<-2/(1-u)^{1/3}#F(x)=1-(2/x)^2
hist(x,prob=TRUE,main=expression(f(x)==8*x^(-3)))
y<-seq(0,20,.01)
lines(y,8*y^(-3))

## -----------------------------------------------------------------------------
n<-1000;k<-0;y<-numeric(n)
while (k<n) {
  u<-runif(1)
  x<-runif(1)#random variate from g(x)=1 and c=12
  if(x^2*(1-x)>u){
    #we accept x
    k<-k+1
    y[k]<-x
  }
}
date<-y
hist(date,prob=TRUE,main=expression(f(y)==12*y^2*(1-y)))
z<-seq(0,1,.01)
lines(z,12*z^2*(1-z))

## -----------------------------------------------------------------------------
n<-1000;r<-4;beta<-2
lambda<-rgamma(n,r,beta)
x<-rpois(n,lambda)#the length of lambda is n

## -----------------------------------------------------------------------------
n<-1000;r<-4;beta<-2
lambda<-rgamma(n,r,beta)
y<-rpois(n,lambda)#the length of lambda is n
hist(y,prob=TRUE,main=expression(f(y)==r*beta^r/(beta+y)^(r+1)))
z<-seq(0,12,.01)
lines(z,r*beta^r/(beta+z)^(r+1))

## -----------------------------------------------------------------------------
quick_sort<-function(x){
      num<-length(x)
      if(num==0||num==1){return(x)
        }else{
          a<-x[1]
          y<-x[-1]
          lower<-y[y<a]
          upper<-y[y>=a]
          return(c(quick_sort(lower),a,quick_sort(upper)))}
    }#define the fast sorting algorithm
n<-c(1e4,2e4,4e4,6e4,8e4)
m<-100;y<-numeric(m);average<-numeric(5)#use average to save the time of all the simulation 
for(k in 1:5){
  #simulate with different n ,5 times
  for(j in 1:m){
    test<-sample(1:n[k])
    y[j]<-system.time(quick_sort(test))[1]
  }
  average[k]<-mean(y)
  #get the average time of siumlation
} 
average

## -----------------------------------------------------------------------------
t<-n*log(n)
lm(average~t)
plot(t,average)
lines(t,average)#scatter plot and regression line

## -----------------------------------------------------------------------------
m<-2e4;x<-runif(m)
theta.hat1<-mean(exp(x))
## estimate θ by simple Monte Carlo method
n<-1e4;y<-runif(n)
theta.hat2<-(mean(exp(y))+mean(exp(1-y)))/2
## estimate θ by antithetic variate approach
theta.hat1
theta.hat2
var1<-var(exp(x))
var2<-var(exp(y)+exp(1-y))
var1
var2
percent<-(var1-var2)/var1
percent
##the percent reduction in variance using the antithetic variate.

## -----------------------------------------------------------------------------
n=1e4;
est <- sd <- numeric(2)
g <- function(x) {
  x^2*exp(-x^2/2)/sqrt(2*pi) * (x > 1)
  }
x<-rnorm(n)#using f1
fg <- g(x)/(exp(-x^2/2)/sqrt(2*pi))
est[1] <- mean(fg)
sd[1] <- sd(fg)
x <- rexp(n, 1) #using f2
fg <- g(x) / exp(-x)
est[2] <- mean(fg)
sd[2] <- sd(fg)
est
sd

## -----------------------------------------------------------------------------
M <- 10000; k <- 10 
r <- M/k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T <- numeric(k)
est <- numeric(N)
g <- function(x) {
  x^2*exp(-x^2/2)/sqrt(2*pi) * (x > 1)
}
for (i in 1:N) {
x<-rexp(M)
for(j in 1:k-1){
  T[j]<-mean(k*g(x)*(x>j)*(x<j+1)/exp(-x))
  }# cut the section into 10 subintervals
T[k]<-mean(k*g(x)*(x>k)/exp(-x))
est[i] <- mean(T)
}
est3<-mean(est)
var3<-sd(est)
est3
var3

## -----------------------------------------------------------------------------
set.seed(22091)
n<-1000
alpha<-0.05
muhat<-numeric(2)
Z<-numeric(1000)#the time of mu in the confidence interval
for(i in 1:1000){
  x<-rlnorm(n)
  muhat[1]<-mean(log(x))-sqrt(var(log(x)))*qt(1-alpha/2,df=n-1)/sqrt(n)
  muhat[2]<-mean(log(x))+sqrt(var(log(x)))*qt(1-alpha/2,df=n-1)/sqrt(n)
  if(muhat[1]<0 && muhat[2]>0){
    Z[i]=1
  }
}
mean(Z)

## -----------------------------------------------------------------------------
count5test <- function(x,y){
  X <- x-mean(x)
  Y <- y-mean(y)
  outx <- (sum(X>max(Y))+ sum(X<min(Y)))
  outy <- (sum(Y>max(X))+ sum(Y<min(X)))
  return(as.integer(max(c(outx,outy))>5))
}#Count Five test

## -----------------------------------------------------------------------------
alpha<-0.055
Ftest <- function(x,y,n){
  varx<-var(x)
  vary<-var(y)
  if(var(x)/var(y)>=qf(alpha/2,n-1,n-1) && var(x)/var(y)<=qf(1-alpha/2,n-1,n-1)){
    return(0)
  }
  else return(1)
}# F test with alpha=0.055 and n1=n2=n

## -----------------------------------------------------------------------------
set.seed(22091)
m<-10000
n<-20
sigma1<-1
sigma2<-1.5
powercount5<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  count5test(x,y)
}))
# the power of Count Five test
powerF<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  Ftest(x,y,n)
}))
# the power of F test
print(powercount5)
print(powerF)

## -----------------------------------------------------------------------------
set.seed(22091)
m<-10000
n<-100
sigma1<-1
sigma2<-1.5
powercount5<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  count5test(x,y)
}))
# the power of Count Five test
powerF<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  Ftest(x,y,n)
}))
# the power of F test
print(powercount5)
print(powerF)

## -----------------------------------------------------------------------------
set.seed(22091)
m<-10000
n<-1000
sigma1<-1
sigma2<-1.5
powercount5<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  count5test(x,y)
}))
# the power of Count Five test
powerF<-mean(replicate(m,expr={
  x<-rnorm(n,0,sigma1)
  y<-rnorm(n,0,sigma2)
  Ftest(x,y,n)
}))
# the power of F test
print(powercount5)
print(powerF)

## -----------------------------------------------------------------------------
x <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
set.seed(22091)
r <- function(x,i){
  1/mean(x[i])
}
library(boot)
obj <- boot(data=x,statistic=r,R=1e4)
# use bootstrap to estimate MLE
round(c(lambda=obj$t0,bias=mean(obj$t)-obj$t0,se=sd(obj$t)),5)

## -----------------------------------------------------------------------------
x <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
set.seed(22091)
r <- function(x,i){
  mean(x[i])
}
library(boot)
obj <- boot(data=x,statistic=r,R=1e4)
# use bootstrap to estimate MLE
round(c(theta=obj$t0,bias=mean(obj$t)-obj$t0,se=sd(obj$t)),5)

## -----------------------------------------------------------------------------
print(boot.ci(obj,type=c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------
n <- 20;m=1000
sum1 <- numeric(m);sum2 <- numeric(m);sum3 <- numeric(m);sum4 <- numeric(m);
set.seed(22091)
r <- function(x,i){
  mean(x[i])
}
library(boot)
for(i in 1:m){
  x <-rnorm(n)
  obj <- boot(data=x,statistic=r,R=1e4)
  cl <- boot.ci(obj,type=c("norm","basic","perc","bca")) # Compute 95% bootstrap confidence intervals
  if(cl$norm[2]<0 && cl$norm[3]>0){
    sum1[i] <-1
  }
  if(cl$basic[4]<0 && cl$basic[5]>0){
    sum2[i] <-1
  }
  if(cl$perc[4]<0 && cl$perc[5]>0){
    sum3[i] <-1
  }
  if(cl$bca[4]<0 && cl$bca[5]>0){
    sum4[i] <-1
  }
}
round(c(norm=mean(sum1),basic=mean(sum2),perc=mean(sum3),bca=mean(sum4)),3)


## -----------------------------------------------------------------------------
rm=list(c())
library(bootstrap)
n <- nrow(scor)
theta.jack <- numeric(n)
#use a function to calculate thetahat
theta <- function(x,i){
  val=eigen(cov(x[i,]))$values
  return(val[1]/sum(val))
}
theta.hat <- theta(scor,1:n)
for(i in 1:n){
  theta.jack[i] <- theta(scor,(1:n)[-i])
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(thetahat=theta.hat,bias.jack=bias.jack,
se.jack=se.jack),4)

## -----------------------------------------------------------------------------
rm=list(c())
library(DAAG)
attach(ironslag)
n<-length(magnetic)
e1<-e2<-e3<-e4<-c()

for (k in 1:(n-1)) {
  for (l in (k+1):n) {
    y<-magnetic[-c(k,l)]
    x<-chemical[-c(k,l)]
    
    J1<-lm(y~x)
    yhat1a<-J1$coef[1]+J1$coef[2]*chemical[k]
    yhat1b<-J1$coef[1]+J1$coef[2]*chemical[l]
    e1<-c(e1,magnetic[k]-yhat1a,magnetic[l]-yhat1b)
    
    J2<-lm(y~x+I(x^2))
    yhat2a<-J2$coef[1]+J2$coef[2]*chemical[k]+J2$coef[3]*chemical[k]^2
    yhat2b<-J2$coef[1]+J2$coef[2]*chemical[l]+J2$coef[3]*chemical[l]^2
    e2<-c(e2,magnetic[k]-yhat2a,magnetic[l]-yhat2b)
    
    J3<-lm(log(y)~x)
    yhat3a<-exp(J3$coef[1]+J3$coef[2]*chemical[k])
    yhat3b<-exp(J3$coef[1]+J3$coef[2]*chemical[l])
    e3<-c(e3,magnetic[k]-yhat3a,magnetic[l]-yhat3b)
    
    J4<-lm(log(y)~log(x))
    yhat4a<-exp(J4$coef[1]+J4$coef[2]*log(chemical[k]))
    yhat4b<-exp(J4$coef[1]+J4$coef[2]*log(chemical[l]))
    e4<-c(e4,magnetic[k]-yhat4a,magnetic[l]-yhat4b)
  }
}
c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))
detach(ironslag)

## -----------------------------------------------------------------------------
rm=list(c())
set.seed(22091)
x <- rnorm(15)
y <- rnorm(15)
R <- 999 #number of replicate
z <- c(x,y)#pooled sample
K <- 1:30
reps <- numeric(R)# storage for replicate
cor0 <- cor(x,y,method="spearman")
for(i in 1:R){
  k <- sample(K,size=15,replace=FALSE)
  x1 <- z[k]
  y1 <- z[-k] #complement of x1
  reps[i] <- cor(x1,y1,method="spearman")
}
p <- mean(c(cor0,reps)>=cor0)
p
cor.test(x,y)

## -----------------------------------------------------------------------------
set.seed(22091)
# use the function below to calculate the density of  Laplace distribution
f <- function(x){
  return(exp(-abs(x))/2)
}
# use the function below to generate chain with given parameters:sigma,n,X_0,N
mh <- function(sigma,N,x0){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for(i in 2:N){
    xt <- x[i-1]
    y <- rnorm(1,xt,sigma)
    num <- f(y)*dnorm(xt,y,sigma)
    den <- f(xt)*dnorm(y,xt,sigma)
    if(u[i] <= num/den) x[i] <- y else{
      x[i] <- xt 
      k <- k+1 # y is rejected
    }
  }
  return(list(x=x,k=k))
}

# generate chains with four differeint sigma

N <- 15000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
mh1 <- mh(sigma[1],N,x0)
mh2 <- mh(sigma[2],N,x0)
mh3 <- mh(sigma[3],N,x0)
mh4 <- mh(sigma[4],N,x0)
# calculate the acceptance rates
print(c((N-mh1$k)/N,(N-mh2$k)/N,(N-mh3$k)/N,(N-mh4$k)/N))

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi){
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i_th row of X
  pis <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi)
  B <- n*var(psi.means)
  psi.w <- apply(psi,1,"var")
  W <- mean(psi.w)
  v.hat <- W*(n-1)/n + (B/n)
  r.hat <- v.hat / W            # G-R statistic
  return(r.hat)
}

k <- 4
b <- 1000
sigma <- 0.5

# choose overdispersed initial values
x0 <- c(-10,-5,5,10)

# generate the chains
X <- matrix(0,nrow=k,ncol=N)
for(i in 1:k){
  X[i,] <- mh(sigma,N,x0[i])$x
}

# compute diagnostic statistics
psi <- t(apply(X,1,cumsum))
for(i in 1:nrow(psi))
  psi[i,] <- psi[i,]/(1:ncol(psi))
print(Gelman.Rubin(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0,N)
for(j in (b+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N],type="l",xlab="",ylab="R")
abline(h=1.2,lty=2)

## -----------------------------------------------------------------------------
# initialize constants and parameters
N <- 5000          # length of chain
burn <- 1000       # burn-in length
X <- matrix(0,N,2) # the chain, a bivariate sample

rho <- 0.9         # correlation
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

##### generate the chain #####

X[1, ] <- c(mu1,mu2)  # initialize
for(i in 2:N){
  x2 <- X[i-1,2]
  m1 <- mu1 + rho*(x2-mu2)*sigma1/sigma2
  X[i,1] <- rnorm(1,m1,s1)
  x1 <- X[i,1]
  m2 <- mu2 + rho*(x1-mu1)*sigma2/sigma1
  X[i,2] <- rnorm(1,m2,s2)
}

# Remove the first 1000 observations and draw the graph generated by the sample
b <- burn + 1
x <- X[b:N, ]

plot(x,main="",cex=.5,xlab=bquote(X[1]),ylab=bquote(X[2]),ylim=range(X[,2]))

lm <- lm(x[,1]~x[,2])
summary(lm)
qqnorm(x)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi){
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i_th row of X
  pis <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi)
  B <- n*var(psi.means)
  psi.w <- apply(psi,1,"var")
  W <- mean(psi.w)
  v.hat <- W*(n-1)/n + (B/n)
  r.hat <- v.hat / W            # G-R statistic
  return(r.hat)
}

k <- 4
b <- 1000
sigma <- 0.6

# choose overdispersed initial values
x0 <- c(-10,-5,5,10)

# generate the chains
X <- matrix(0,nrow=k,ncol=N)
for(i in 1:k){
  X[i,] <- mh(sigma,N,x0[i])$x
}

# compute diagnostic statistics
psi <- t(apply(X,1,cumsum))
for(i in 1:nrow(psi))
  psi[i,] <- psi[i,]/(1:ncol(psi))
print(Gelman.Rubin(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0,N)
for(j in (b+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N],type="l",xlab="",ylab="R")
abline(h=1.2,lty=2)

## -----------------------------------------------------------------------------
set.seed(22091)
root <- function(N,b1,b2,b3,f0){
  x1 <- rpois(N,1); x2 <- rexp(N,1); x3 <- rbinom(N,1,0.5)
  g <- function(alpha){
    tmp <- exp(alpha+b1*x1+b2*x2+b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
    }
  solution <- uniroot(g,c(0,10))
  solution$root
}
N <- 1e6; b1 <- 0; b2 <- 1; b3  <- -1; f0 <-c(.1,.01,.001,.0001)
alpha <- numeric(4)
for(i in 1:4){
  alpha[i] <- root(N,b1,b2,b3,f0[i])
}
alpha
plot(f0,alpha,xlim=c(0,0.1),ylim=c(0,10))

## -----------------------------------------------------------------------------
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- c(12,9,28,14,17,1,24,11,25,3)
x <- numeric(10)
obj <- function(lambda,u,v){
  for(i in 1: 10){
    x[i] <- (v[i]-u[i])/(exp(lambda*(v[i]-u[i])))-u[i]
  }
  sum(x)
}
u <- uniroot(obj,lower=-10,upper=10,u= u,v= v)
lambda.hat <- u$root
lambda.hat

## -----------------------------------------------------------------------------
data <- data.frame(x=1:5,y=6:10)
scale01 <- function(x) {
     rng <- range(x, na.rm = TRUE)
     (x - rng[1]) / (rng[2] - rng[1])
}
lapply(data,scale01)

## -----------------------------------------------------------------------------
data <- data.frame(x=c(1.1,2.2,3.3),y=c("a","b","c"),z=c(4.4,5.5,6.6))
a=c()
for(i in 1:ncol(data)){
  if(class(data[,i])!='numeric') a=c(a,i)
}
data2=data[,-a]
lapply(data2,scale01)

## -----------------------------------------------------------------------------
data <- data.frame(x=c(1.1,2.2,3.3),y=c(1.1,3.3,5.5))
vapply(data,sd,c(1))

## -----------------------------------------------------------------------------
data <- data.frame(x=c(1.1,2.2,3.3),y=c(1.1,3.3,5.5),z=c('a','b','c'))
a <-vapply(data,is.numeric,c(1))
b <- c()
for(i in 1:length(a)){
  if(a[i]!=1) b=c(b,i)
}
data2=data[,-b]
vapply(data2,sd,c(1))

## -----------------------------------------------------------------------------
N <- 5000
burn <- 1000 # burn-in length
X <- matrix(0,N,2)

rho <-0.9 # correlation
mu1 <- mu2 <-0
sigma1 <- sigma2 <-1


# generate the chain
gibbs_r <- function(N,burn,X,rho,mu1,mu2,sigma1,sigma2){
  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2
  X[1,] <- c(mu1,mu2)
  for(i in 2:N){
    x2 <- X[i-1,2]
    m1 <- mu1+rho*(x2-mu2)*sigma1/sigma2
    X[i,1] <- rnorm(1,m1,s1)
    x1 <- X[i,1]
    m2 <- mu2+rho*(x1-mu1)*sigma2/sigma1
    X[i,2] <- rnorm(1,m2,s2)
  }
  
  b <- burn +1 
  x <- X[b:N,]
}

x1 <- gibbs_r(N,burn,X,rho,mu1,mu2,sigma1,sigma2)

## -----------------------------------------------------------------------------
N <- 5000
burn <- 1000 # burn-in length
X <- matrix(0,N,2)

rho <-0.9 # correlation
mu1 <- mu2 <-0
sigma1 <- sigma2 <-1
library(Rcpp)
sourceCpp('C:/gibbs_cpp.cpp')
x2 <- gibbs_cpp(N,mu1,mu2,sigma1,sigma2,rho)

## -----------------------------------------------------------------------------
qqplot(x1,x2)

## -----------------------------------------------------------------------------
library(microbenchmark)
microbenchmark(gibbs_r(N,burn,X,rho,mu1,mu2,sigma1,sigma2),gibbs_cpp(N,mu1,mu2,sigma1,sigma2,rho))

