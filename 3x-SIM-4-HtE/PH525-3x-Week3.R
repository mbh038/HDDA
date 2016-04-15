## PH525-3x

## Statistical Inference and Modelling for High-throughput Experiments

## Week 3 - Statistical Models, Maximum Likelihood, Models for Variance


## Statistical Models Exercises

# The probability of conceiving a girl is 0.49. What is the probability that a family
# with 4 children has 2 girls and 2 boys (you can assume no twins)?

dbinom(2,4,0.49)

# What is the probability that a family with 10 children has 4 girls and 6 boys
# (you can assume no twins)?

dbinom(4,10,0.49)

# The genome has 3 billion bases. About 20% are C, 20% are G, 30% are T and 30% are A.
# Suppose you take a random interval of 20 bases, what is the probability that the 
# GC-content (proportion of Gs or Cs) is strictly above 0.5 in this interval
# (you can assume independence)?

1-pbinom(10,20,0.4)

# The probability of winning the lottery is 1 in 175,223,510. If 189,000,000 randomly
# generated (with replacement) tickets are sold, what is the probability that at least
# one winning tickets is sold? (give your answer as a proportion not percentage)

prob<-1/175223510
N<-189000000
1-pbinom(0,N,prob) # 1 or more winning tickets
1-pbinom(1,N,prob) # 2 or more

# The genome has 3 billion bases. About 20% are C, 20% are G, 30% are T and 30% are A.
# Suppose you take a random interval of 20 bases, what is the exact probability that
# the GC-content (proportion of Gs of Cs) is greater than 0.35 and smaller or equal to
# 0.45 in this interval? HINT: use the binomial distribution.

pbinom(9,20,0.4)-pbinom(7,20,0.4)

# For the question above, what is the normal approximation to the probability?

b <- (9 - 20*.4)/sqrt(20*.4*.6)
a <- (7 - 20*.4)/sqrt(20*.4*.6)
pnorm(b)-pnorm(a)

# Repeat # 3 but with interval of 1000 bases

exact<-pbinom(450,1000,0.4)-pbinom(350,1000,0.4)
b<-(450 - 1000*.4)/sqrt(1000*.4*.6)
a<-(350 - 1000*.4)/sqrt(1000*.4*.6)
normApprox<-pnorm(b)-pnorm(a)
abs(exact-normApprox)



Ns <- c(5,10,30,100)
ps <- c(0.01,0.10,0.5,0.9,0.99)
library(rafalib)
mypar(4,5)
for (N in Ns){
    k<-seq(1,N-1)
    for (p in ps){
        
        exact = dbinom(k,N,p) 
        
        a <- (k+0.5 - N*p)/sqrt(N*p*(1-p))
        b <- (k-0.5 - N*p)/sqrt(N*p*(1-p))
        approx = pnorm(a) - pnorm(b)

        LIM <- range(c(approx,exact))
        plot(exact,approx,main=paste("N =",N," p = ",p),xlim=LIM,ylim=LIM,col=1,pch=16)
        abline(0,1)
    }
}
# lottery - binomial, normal and poisson approximations

# probability of exactly two people winning
N <- 189000000
p <- 1/175223510
dbinom(2,N,p)

# normal approximation to this
a <- (2+0.5 - N*p)/sqrt(N*p*(1-p))
b <- (2-0.5 - N*p)/sqrt(N*p*(1-p))
pnorm(a) - pnorm(b)

# This over estimates the probability - expect this since p is so small

# use poisson approximation, with lambda = N*p = tickets per N sold that win
dpois(2,N*p)

# N is very very large and Np is not 0, so poisson is very good approximation to binomial


1-ppois(1,N*p)


## Maximum Likelihood Exercises

# get latest version of the dagdata library:
library(devtools)
install_github("genomicsclass/dagdata")

# load the palindrome data from the Human cytomegalovirus genome:
library(dagdata)
data(hcmv)

# These are the locations of palindromes on the genome of this virus:
library(rafalib)
mypar()
plot(locations,rep(1,length(locations)),ylab="",yaxt="n")

# These palindromes are quite rare, p is very small. If we break the genome into
# bins of 4000 basepairs, then we have Np not so small and we might be able to
# use Poisson to model the number of palindromes in each bin:
breaks=seq(0,4000*round(max(locations)/4000),4000)
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))

# looks like Poisson would be  reasonable model
hist(counts)

probs <- dpois(counts,4)
likelihood <- prod(probs)
likelihood

logprobs <- dpois(counts,4,log=TRUE)
loglikelihood <- sum(logprobs)
loglikelihood

logll<-function(lambda,x){
    sum(dpois(x,lambda,log=TRUE))
}
lambdas = seq(1,15,len=300)
loglls<-sapply(lambdas,function(lambda) logll(lambda,counts))
plot(lambdas,loglls,type="l",col="blue")
maxll<-lambdas[which.max(loglls)]
maxll
abline(v=maxll)


breaks=seq(0,4000*round(max(locations)/4000),4000)
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))
binLocation=(breaks[-1]+breaks[-length(breaks)])/2
plot(binLocation,counts,type="l",xlab=)

binLocation[which.max(counts)]
max(counts)

lambda = mean(counts[ - which.max(counts) ])

# probability of seeing a count of 14 or more
1-ppois(13,lambda)

# this p value is less than 0.05, but We selected the highest region 
# out of 57 and need to adjust for multiple testing. 

# Use the Bonferroni correction to determine the p-value cut-off that
# guarantees a FWER of 0.05. What is this p-value cutoff?
0.05/57

# Our p-value does satisfy the Bonferroni correction

# Create a qq-plot to see if our Poisson model is a good fit:
ps <- (seq(along=counts) - 0.5)/length(counts)
lambda <- mean( counts[ -which.max(counts)])
poisq <- qpois(ps,lambda)
qqplot(poisq,counts)
abline(0,1)

## MODELS FOR VARIANCE

library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)

data("tissuesGeneExpression")
library(genefilter)
y = e[,which(tissue=="endometrium")]

s2 <- rowVars(y) # between sample variances
datsds<-rowSds(y) # between sample standard devations
library(rafalib)
mypar(1,2)
qqnorm(s2)
qqline(s2)
##To see the square root transformation does not help much:
qqnorm(sqrt(s2))
qqline(sqrt(s2))

#Now fit an F-distribution with 14 df

library(limma)
estimates=fitFDist(datsds^2,14)
estimates$scale

#qq plots
ps <- (seq(along=s2)-0.5)/length(s2)
theoretical<- sqrt(qf(ps, 14, estimates$df2)*estimates$scale)
observed <- datsds

mypar(1,2)
LIM <- range(c(theoretical,observed)) 
qqplot(theoretical,observed,ylim=LIM,xlim=LIM)
abline(0,1)

##close up excluding the upper 5% - is excellent fit
K <- quantile(observed,0.95) 
qqplot(theoretical,observed,ylim=c(0,K),xlim=c(0,K))
abline(0,1)

#histogram
sds=seq(0,2,len=100)
tmp=hist(observed,main=paste("s_0 =", signif(estimates[[1]],2), "d =",
         signif(estimates[[2]],2)), xlab="sd", ylab="density",
         freq=FALSE, nc=100, xlim=c(0,2), ylim=c(0,4))
dd=df(sds^2/estimates$scale,14,estimates$df2)
k=sum(tmp$density)/sum(dd) ##a normalizing constant to assure same area in plot
lines(sds, dd*k, type="l", col=2, lwd=2)
