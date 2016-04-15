##PH525.2x

## Week 4 exercises


list.of.packages <- c("contrast", "UsingR","rafalib")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(contrast)
library(rafalib)

species <- factor(c("A","A","B","B"))
condition <- factor(c("control","treated","control","treated"))

SCdata=cbind(species,condition)
X<-model.matrix(~ species + condition)

imagemat(X)

y = rnorm(4)

fit = lm(y ~ species + condition)

contrast(fit, list(species="B",condition="control"), list(species="A",condition="treated"))$X

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)

fitTL <- lm(friction ~ type + leg, data=spider)
summary(fitTL)

(coefs <- coef(fitTL))

Y <- spider$friction
X <- model.matrix(~ type + leg, data=spider)
beta.hat <- solve(t(X) %*% X) %*% t(X) %*% Y
t(beta.hat)

coefs


library(contrast) #Available from CRAN
L3vsL2 <- contrast(fitTL,list(leg="L3",type="pull"),list(leg="L2",type="pull"))
L3vsL2
L4vsL2 <- contrast(fitTL,list(leg="L4",type="pull"),list(leg="L2",type="pull"))
L4vsL2

Sigma.hat <- sum(fitTL$residuals^2)/(nrow(X) - ncol(X)) * solve(t(X) %*% X)
signif(Sigma.hat, 2)

(cT <- L4vsL2$X)
sqrt(cT %*% Sigma.hat %*% t(cT))

L4vsL2$SE


X <- model.matrix(~ type + leg, data=spider)
(Sigma <- sum(fitTL$residuals^2)/(nrow(X) - ncol(X)) * solve(t(X) %*% X))

C <- matrix(c(0,0,-1,0,1),1,5)

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)

spider$log2friction <- log2(spider$friction)
boxplot(log2friction ~ type*leg, data=spider)

fitlogLM<-lm(log2friction~type*leg,data=spider)
summary(fitlogLM)


anova(fitlogLM)

# What is the L2 vs L1 estimate in log2friction for the pull samples?
coef(fitlogLM)["legL2"]


# What is the L2 vs L1 estimate in log2friction for the push samples?
L2push.vs.L1push <- contrast(fitlogLM,
                           list(leg="L2", type = "push"), 
                           list(leg="L1", type = "push"))
L2push.vs.L1push

# '' which is the legL2 coefficient plus the push:L2 interaction:
coef(fitlogLM)["legL2"] + coef(fitlogLM)["typepush:legL2"]


# In the video we briefly mentioned the analysis of variance (or ANOVA,
# performed in R using the anova() function), which allows us to test
# whether a number of coefficients are equal to zero, by comparing a 
# linear model including these terms to a linear model where these terms
# are set to 0.

# The book page for this section has a section,
# "Testing all differences of differences", which explains the ANOVA
# concept and the F-test in some more detail. You can read over that
# section before or after the following question.

# In this last question, we will use Monte Carlo techniques to observe
# the distribution of the ANOVA's "F-value" under the null hypothesis,
# that there are no differences between groups.

# Suppose we have 4 groups, and 10 samples per group, so 40 samples
# overall:
    
N <- 40
p <- 4
group <- factor(rep(1:p,each=N/p))
X <- model.matrix(~ group)

# We will show here how to calculate the "F-value", and then we will use
# random number to observe the distribution of the F-value under the null
# hypothesis.

# The F-value is the mean sum of squares explained by the terms of 
# interest (in our case, the 'group' terms) divided by the mean sum
# of squares of the residuals of a model including the terms of 
# interest. So it is the explanatory power of the terms divided 
# by the leftover variance.

# Intuitively, if this number is large, it means that the group 
# variable explains a lot of the variance in the data, compared to
# the amount of variance left in the data after using group information. We will calculate these values exactly here:

# generate random, null data where the mean is the same for all groups
Y <- rnorm(N,mean=42,7)

# The base model we will compare against is simply Y-hat = mean(Y), 
# which we will call mu0, and the initial sum of squares is the 
# Y values minus mu0:
mu0 <- mean(Y)
initial.ss <- sum((Y - mu0)^2)

# We then need to calculate the fitted values for each group, 
# which is simply the mean of each group, and the residuals from this
# model, which we will call "after.group.ss" for the sum of squares
# after using the group information:
    
s <- split(Y, group)
after.group.ss <- sum(sapply(s, function(x) sum((x - mean(x))^2)))

# Then the explanatory power of the group variable is the initial sum of
# squares minus the residual sum of squares:
(group.ss <- initial.ss - after.group.ss)

# We calculate the mean of these values, but we divide by terms which 
# remove the number of fitted parameters. For the group sum of squares,
# this is number of parameters used to fit the groups (3, because the
# intercept is in the initial model). For the after group sum of 
# squares, this is the number of samples minus the number of parameters
# total (So N - 4, including the intercept).

group.ms <- group.ss / (p - 1)
after.group.ms <- after.group.ss / (N - p)

# The F-value is simply the ratio of these mean sum of squares.

f.value <- group.ms / after.group.ms

# Set the seed to 1, set.seed(1) then calculate the F-value for
# 1000 random versions of Y. What is the mean of these F-values?



N <- 40
p <- 4
group <- factor(rep(1:p,each=N/p))
X <- model.matrix(~ group)

numTrials<-1000
set.seed(1)
Fs=replicate(numTrials,{
    Y <- rnorm(N,mean=42,7)
    mu0 <- mean(Y)
    initial.ss <- sum((Y - mu0)^2)
    s <- split(Y, group)
    after.group.ss <- sum(sapply(s, function(x) sum((x - mean(x))^2)))
    group.ss <- initial.ss - after.group.ss
    group.ms <- group.ss / (p - 1)
    after.group.ms <- after.group.ss / (N - p)
    f.value <- group.ms / after.group.ms
    return(f.value)
})
mean(Fs)

# Plot the distribution of the 1000 F-values:
hist(Fs, col="grey", border="white", breaks=50, freq=FALSE)

# Overlay the theoretical F-distribution, with parameters 
# df1=p - 1, df2=N - p.
xs <- seq(from=0,to=6,length=100)
lines(xs, df(xs, df1 = p - 1, df2 = N - p), col="red")

# This is the distribution which is used to calculate the p-values
# for the ANOVA table produced by anova(). 

## COLLINEARITY

# The rank of a matrix is the number of columns that are
#independent of all the others.

# qr() is the rank function

# matrix for which rank = ncol called a full rank matrix +> no collinearity

# if rank < ncol then one of more columns are collinear

# In this example, C and D are collinear with Sex
Sex <- c(0,0,0,0,1,1,1,1)
A <-   c(1,1,0,0,0,0,0,0)
B <-   c(0,0,1,1,0,0,0,0)
C <-   c(0,0,0,0,1,1,0,0)
D <-   c(0,0,0,0,0,0,1,1)
X <- model.matrix(~Sex+A+B+C+D-1)
cat("ncol=",ncol(X),"rank=", qr(X)$rank,"\n")

# We can remove the collinearity here by balancing the sexes across the treatments:

Sex <- c(0,1,0,1,0,1,0,1)
A <-   c(1,1,0,0,0,0,0,0)
B <-   c(0,0,1,1,0,0,0,0)
C <-   c(0,0,0,0,1,1,0,0)
D <-   c(0,0,0,0,0,0,1,1)
X <- model.matrix(~Sex+A+B+C+D-1)
cat("ncol=",ncol(X),"rank=", qr(X)$rank,"\n")

# Let's use the example from the lecture to visualize how there is not
# a single best beta-hat, when the design matrix has collinearity 
# of columns.

sex <- factor(rep(c("female","male"),each=4))
trt <- factor(c("A","A","B","B","C","C","D","D"))

X <- model.matrix( ~ sex + trt)
qr(X)$rank
Y <- 1:8

a<-1
b<-2

makeYstar <- function(a,b) Y - X[,2] * a - X[,5] * b

fitTheRest <- function(a,b) {
    Ystar <- makeYstar(a,b)
    Xrest <- X[,-c(2,5)]
    betarest <- solve(t(Xrest) %*% Xrest) %*% t(Xrest) %*% Ystar
    residuals <- Ystar - Xrest %*% betarest
    sum(residuals^2)
}

# What is the sum of squared residuals when the male coefficient is
# 1 and the D coefficient is 2, and the other coefficients are fit
# using the linear model solution?

fitTheRest(1,2)

betas = expand.grid(-2:8,-2:8)
rss = apply(betas,1,function(x) fitTheRest(x[1],x[2]))

## Note that all pairs add to 6
themin= min(rss)
betas[which(rss==themin),]

library(rafalib)
## plot the pairs what are minimum
themin=min(rss)
plot(betas[which(rss==themin),])

# There is clearly not a single beta which optimizes the least squares
# equation, due to collinearity, but an infinite line of solutions
# which produce an identical sum of squares values.

imagemat(rss)
