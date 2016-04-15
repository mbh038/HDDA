## PH525-1x-Week3-Exercises.R

## T-TEST EXERCISES

library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)

library(dplyr)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist

library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)

set.seed(1)
dat.ns<-sample(bwt.nonsmoke,25)
dat.s<-sample(bwt.smoke,25)
tval<-t.test(dat.ns,dat.s)$statistic
tval

# doing this manually (without t.test

N=25
set.seed(1)
dat.ns <- sample(bwt.nonsmoke , N)
dat.s <- sample(bwt.smoke , N)

X.ns <- mean(dat.ns)
sd.ns <- sd(dat.ns) # note use of sd,not popsd since this is a sample.

X.s <- mean(dat.s)
sd.s <- sd(dat.s)

sd.diff <- sqrt(sd.ns^2/N+sd.s^2/N)
tval <- (X.ns - X.s)/sd.diff
tval

pnorm(tval,lower.tail=FALSE)*2

qnorm(.995)*sd.diff

## CONFIDENCE INTERVALS EXERCISES

library (downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)

library(dplyr)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist

library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)

N<-25
set.seed(1)
dat.ns=sample(bwt.nonsmoke,N)
dat.s=sample(bwt.smoke,N)

X.ns <- mean(dat.ns)
sd.ns <- sd(dat.ns) # note use of sd,not popsd since this is a sample.

X.s <- mean(dat.s)
sd.s <- sd(dat.s)

sd.diff <- sqrt(sd.ns^2/N+sd.s^2/N)

tval <- (X.ns - X.s)/sd.diff
tval

tval<-t.test(dat.s,dat.ns)$statistic
tval
2*pt(-abs(tval),48)

qt(.995,48)*sd.diff

N<-5
set.seed(1)
dat.ns=sample(bwt.nonsmoke,N)
dat.s=sample(bwt.smoke,N)
t.test(dat.s,dat.ns)

## POWER CALCULATIONS

library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)

library(dplyr)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist

library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)

numTrials<-10000
N<-5
set.seed(1)
reject <- function(N, alpha=0.01){
  dat.ns=sample(bwt.nonsmoke,N) 
  dat.s=sample(bwt.smoke,N)
  pval<-t.test(dat.s,dat.ns)$p.value
  pval < alpha
}
rejections <- replicate(numTrials,reject(N))
mean(rejections)

numTrials<-10000
Ns<-c(30,60,90,120)
set.seed(1)
power<-sapply(Ns,function(N){
  rejections<-replicate(numTrials,reject(N))
  mean (rejections)
  })

plot(Ns,power,type="b")

## MONTE CARLO EXERCISES

set.seed(1)
n<-500
B<-100
ttestgenerator <- function(n) {
  X<-rnorm(n)
  tstat<-sqrt(n)*mean(X)/sd(X)
  return(tstat)
}
ttests <- replicate(B, ttestgenerator(n))

mean(ttests>2)

ps = seq(1/(B+1), 1-1/(B+1),len=B)
qs = qt(ps,df=n-1)
qsn =qnorm(ps)
mypar(1,2)
qqplot(qs,ttests)
abline(0,1)
qqplot(qsn,ttests)
abline(0,1)

library(rafalib)
mypar(3,2)

Ns<-seq(5,30,5)
B <- 1000
mypar(3,2)
LIM <- c(-4.5,4.5)
for(N in Ns){
  ts <- replicate(B, {
    X <- rnorm(N)
    sqrt(N)*mean(X)/sd(X)
  })
  ps <- seq(1/(B+1),1-1/(B+1),len=B)
  qqplot(qt(ps,df=N-1),ts,main=N,
         xlab="Theoretical",ylab="Observed",
         xlim=LIM, ylim=LIM)
  abline(0,1)
} 

Ns<-seq(5,30,5)
B <- 1000
mypar(3,2)
LIM <- c(-4.5,4.5)
for(N in Ns){
  tstats <- replicate(B, {
    X1 <- rnorm(N)
    X2 <- rnorm(N)
    t.test(X1,X2,var.equal=TRUE)$stat
  })
  ps <- seq(1/(B+1),1-1/(B+1),len=B)
  qqplot(qt(ps,df=2*N-2),tstats,main=N,
         xlab="Theoretical",ylab="Observed",
         xlim=LIM, ylim=LIM)
  abline(0,1)
} 

set.seed(1)
N <- 1000
B <- 10000
tstats <- replicate(B,{
  X <- sample(c(-1,1), N, replace=TRUE)
  sqrt(N)*mean(X)/sd(X)
})
ps=seq(1/(B+1), 1-1/(B+1), len=B) 
qqplot(qt(ps,N-1), tstats, xlim=range(tstats))
abline(0,1)
#The population data is not normal thus the theory does not apply.
#We check with a Monte Carlo simulation. The qqplot shows a large tail. 
#Note that there is a small but positive chance that all the X are the same.
##In this case the denominator is 0 and the t-statistics is not defined

set.seed(1)
N <- 1000
B <- 10000
tstats <- replicate(B,{
  X <-  sample(c(-1,1), N, replace=TRUE)
  sqrt(N)*mean(X)/sd(X)
})
qqnorm(tstats)
abline(0,1)
#With N=1000, CLT kicks in and the t-statistic is approximated with normal 0,1
##Furthermore, t-distribution with df=999 and normal are practically the same.

mypar(3,2)
set.seed(1)
Ns <- seq(5,45,5)
library(rafalib)
mypar(3,3)
for(N in Ns){
  medians <- replicate(10000, median ( rnorm(N) ) )
  title <- paste("N=",N,", avg=",round( mean(medians), 2) , ", sd*sqrt(N)=", round( sd(medians)*sqrt(N),2) )
  qqnorm(medians, main = title )
  qqline(medians)
}
##there is an asymptotic result that says SD is sqrt(N*4*dnorm(0)^2)

## PERMUTATIONS
library(downloader)
library(dplyr)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist

N=10
set.seed(1)
nonsmokers <- sample(bwt.nonsmoke , N)
smokers <- sample(bwt.smoke , N)
obs <- mean(smokers) - mean(nonsmokers)

dat <- c(smokers,nonsmokers)
set.seed(1)
null <- replicate(1000, {
  shuffle <- sample(dat)
  smokersstar <- shuffle[1:N]
  nonsmokersstar <- shuffle[(N+1):(2*N)]
  return(mean(smokersstar)-mean(nonsmokersstar))
})
hist(null)
abline(v=obs, col="red", lwd=2)

#the proportion of permutations with larger difference
(sum(abs(null) > abs(obs)) + 1) / (length(null) + 1)

# Repeat the above exercise, but instead of the differences in mean, consider the differences
# in median obs <- median(smokers) - median(nonsmokers). What is the permutation based p-value?

N=10
set.seed(1)
nonsmokers <- sample(bwt.nonsmoke , N)
smokers <- sample(bwt.smoke , N)
obs <- median(smokers) - median(nonsmokers)


set.seed(1)
null <- replicate(1000, {
  dat <- c(smokers,nonsmokers)
  shuffle <- sample(dat)
  smokersstar <- shuffle[1:N]
  nonsmokersstar <- shuffle[(N+1):(2*N)]
  return(median(smokersstar)-median(nonsmokersstar))
})
hist(null)
abline(v=obs, col="red", lwd=2)

#the proportion of permutations with larger difference
(sum(abs(null) > abs(obs)) + 1) / (length(null) + 1)
(sum(abs(null) > abs(obs))) / (length(null))


## ASSOCIATION TESTS

# lab / textbook example

disease=factor(c(rep(0,180),rep(1,20),rep(0,40),rep(1,10)),
               labels=c("control","cases"))
genotype=factor(c(rep("AA/Aa",200),rep("aa",50)),
                levels=c("AA/Aa","aa"))
dat <- data.frame(disease, genotype)
dat <- dat[sample(nrow(dat)),] #shuffle them up
head(dat)

tab <- table(genotype,disease)
tab

# odds ratio
(tab[2,2]/tab[2,1]) / (tab[1,2]/tab[1,1])

# overall proportion that have the disease
p=mean(disease=="cases")
p

# expected table
expected <- rbind(c(1-p,p)*sum(genotype=="AA/Aa"),
                  c(1-p,p)*sum(genotype=="aa"))
dimnames(expected)<-dimnames(tab)
expected

# chi sq test
chisq.test(tab)$p.value

## ASSOCIATION TESTS EXERCISES

d = read.csv("assoctest.csv")
str(d)
tab <- table(d$allele,d$case)
tab
chisq.test(tab)
fisher.test(tab)