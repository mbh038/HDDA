##PH525.2x

## Week 2 exercises

#install.packages("UsingR")

## Falling object

g <- 9.8 ##meters per second
n <- 25
tt <- seq(0,3.4,len=n) ##time in secs, t is a base function
f <- 56.67  - 0.5*g*tt^2
y <-  f + rnorm(n,sd=1)

plot(tt,y,ylab="Distance in meters",xlab="Time in seconds")
lines(tt,f,col=2)

X=cbind(1,tt,tt^2)
Beta=matrix(c(55,0,-5),3,1)

r=y-X%*%Beta

RSS=t(r)%*%r

betahat=solve(crossprod(X))%*%crossprod(X,y)

## MATRIX ALGEBRA EXAMPLES EXERCISES

# Suppose we are analyzing a set of 4 samples. The first two samples are from
# a treatment group A and the second two samples are from a treatment group B.
# This design can be represented with a model matrix like so:
  
X <- matrix(c(1,1,1,1,0,0,1,1),nrow=4)
rownames(X) <- c("a","a","b","b")

# fitted parameters for a linear model give

beta <- c(5, 2)

X%*%beta

# Control a, two treatnents b and c
X <- matrix(c(1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,1,1),nrow=6)
rownames(X) <- c("a","a","b","b","c","c")

beta <- c(10,3,-3) # fitted betas

X%*%beta

## Inference

# Sampling Distribution

g = 9.8 ## meters per second
h0 = 56.67
v0 = 0
n = 25
tt = seq(0,3.4,len=n) ##time in secs, t is a base function
y = h0 + v0 *tt - 0.5* g*tt^2 + rnorm(n,sd=1)
X = cbind(1,tt,tt^2)
A = solve(crossprod(X))%*%t(X)
-2 * (A %*% y) [3]

numTrials<-100000
set.seed(1)
gfind <- function(){
  
  y = h0 + v0 *tt - 0.5* g*tt^2 + rnorm(n,sd=1)
  X = cbind(1,tt,tt^2)
  A = solve(crossprod(X))%*%t(X)
  gval<- -2 * (A %*% y) [3]
}
gg <- replicate(numTrials,gfind())
mean(gg)
sd(gg)

library(UsingR)
x = father.son$fheight
y = father.son$sheight
n = length(y)


numTrials<-10000
N =  50
set.seed(1)
betahats<-replicate(numTrials,{
    index = sample(n,N)
    sampledat = father.son[index,]
    x = sampledat$fheight
    y = sampledat$sheight
    betahat =  lm(y~x)$coef
    betahat[2]
    })
sd(betahats)

mean( (y - mean(y))*(x-mean(x) ) )

## Standard Errors Exercises #1

library(UsingR)

x = father.son$fheight
y = father.son$sheight

n = length(y)
N = 50
set.seed(1)
index = sample(n,N)
sampledat = father.son[index,]
x = sampledat$fheight
y = sampledat$sheight
betahat = lm(y~x)$coef

#se(betahat) = sqrt(var(betahat))
#var(betahat) = sigma^2 * solve((t(X) %*% X))

fit = lm(y ~ x)
fit$fitted.values

r_i = Y_i - Y-hat_i)

# sum of squared residuals

r=y-fit$fitted.values
SSR<-sum(r^2)

sigma2<-SSR/(N-length(fit$coef))

# design matrix
X = cbind(rep(1,N), x) # could also use X<-model.matrix(~x)
solve((t(X) %*% X))

# variance of betahat

var_betahat<-sigma2*diag(solve((t(X) %*% X)))
se_betahat<-sqrt(var_betahat)
se_betahat

summary(fit)
