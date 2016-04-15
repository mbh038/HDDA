## Random Variables Exercises

library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- paste("../data/",basename(url),sep="")
download(url, destfile=filename)
x <- unlist( read.csv(filename) )

mean(x)

set.seed(1)
s1<-sample(x,5)
abs(mean(x)-mean(s1))

set.seed(5)
s2<-sample(x,5)
abs(mean(x)-mean(s2))

set.seed(1)
n <- 1000
sAve <- vector("numeric",n)
for (i in 1:n) {
  s5 <- sample(x,5)
  sAve[i] <- mean(s5)
}
hist(sAve)
mean(abs(sAve-mean(x))>1)

set.seed(1)
n <- 10000
sAve <- vector("numeric",n)
for (i in 1:n) {
  s5 <- sample(x,5)
  sAve[i] <- mean(s5)
}
hist(sAve)
mean(abs(sAve-mean(x))>1)

set.seed(1)
n <- 1000
sAve <- vector("numeric",n)
for (i in 1:n) {
  s50 <- sample(x,50)
  sAve[i] <- mean(s50)
}
hist(sAve)
mean(abs(sAve-mean(x))>1)

## Probability Distributions

#install.packages("gapminder")

library(gapminder)
data(gapminder)
head(gapminder)

x<-gapminder$lifeExp


dat1952 = gapminder[ gapminder$year == 1952, ]
x = dat1952$lifeExp
mean(x <= 40) # finds proportion of x that is less than or equal to 40.

dat1952 = gapminder[ gapminder$year == 1952, ]
x = dat1952$lifeExp
mean(x <= 60) - mean(x <= 40) # finds proportion of x that is between 40 and 60.

# proportion function
prop = function(q) {
  mean(x <= q)
}

props = sapply(qs, prop)
plot(qs, props)
props = sapply(qs, function(q) mean(x <= q))

## Central Limit Theorem

## NORMAL DISTRIBUTION EXERCISES

library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
x <- unlist( read.csv(filename) )

set.seed(1)
n <- 1000
averages5 <- vector("numeric",n)
for(i in 1:n){
  X <- sample(x,5)
  averages5[i] <- mean(X)
}
hist(averages5) ##take a look
mean( abs( averages5 - mean(x) ) > 1)

set.seed(1)
n <- 1000
averages50 <- vector("numeric",n)
for(i in 1:n){
  X <- sample(x,50)
  averages50[i] <- mean(X)
}
hist(averages50) ##take a look
mean( abs( averages50 - mean(x) ) > 1)

mean(averages50 <= 25) - mean(averages50 <= 23)
mean( averages50 < 25 & averages50 > 23) # same thing

qnorm(.975)
pnorm((25-23.9)/0.43)-pnorm((23-23.9)/0.43)

## POPULATION, SAMPLES, AND ESTIMATES EXERCISES

library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- read.csv(filename) 
#We will remove the lines that contain missing values:
dat <- na.omit( dat )
str(dat)

## Males from populations

library(dplyr)
# Use dplyr to create a vector x with the body weight of all males on the control (chow) diet. What is this population's average?
x<-dat %>% filter(Diet=="chow" & Sex=="M") %>% select(Bodyweight) %>% unlist()
mean(x)

# find x population standard deviation
library(rafalib)
popsd(x)

#Set the seed at 1. Take a random sample X of size 25 from x. What is the sample average?
set.seed(1)
X<-sample(x,25)
mean(X)

library(dplyr)
# Use dplyr to create a vector y with the body weight of all males on the high fat (hf) diet. What is this population's average?
y<-dat %>% filter(Diet=="hf" & Sex=="M") %>% select(Bodyweight) %>% unlist()
mean(y)

# find y population standard deviation
library(rafalib)
popsd(y)

# Set the seed at 1. Take a random sample Y of size 25 from y. What is the sample average?
set.seed(1)
Y<-sample(y,25)
mean(Y)

# What is the difference in absolute value between y¯???x¯ and X¯???Y¯?
abs(abs(mean(x)-mean(y))-abs(mean(X)-mean(Y)))

## Repeat for Females from populations

library(dplyr)
# Use dplyr to create a vector x with the body weight of all males on the control (chow) diet. What is this population's average?
x<-dat %>% filter(Diet=="chow" & Sex=="F") %>% select(Bodyweight) %>% unlist()
mean(x)

# find x population standard deviation
library(rafalib)
popsd(x)

#Set the seed at 1. Take a random sample X of size 25 from x. What is the sample average?
set.seed(1)
X<-sample(x,25)
mean(X)

library(dplyr)
# Use dplyr to create a vector y with the body weight of all males on the high fat (hf) diet. What is this population's average?
y<-dat %>% filter(Diet=="hf" & Sex=="F") %>% select(Bodyweight) %>% unlist()
mean(y)

# find y population standard deviation
library(rafalib)
popsd(y)

# Set the seed at 1. Take a random sample Y of size 25 from y. What is the sample average?
set.seed(1)
Y<-sample(y,25)
mean(Y)

# What is the difference in absolute value between y¯???x¯ and X¯???Y¯?
abs(abs(mean(x)-mean(y))-abs(mean(X)-mean(Y)))

## CENTRAL LIMIT THEOREM EXERCISES

library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- na.omit( read.csv(filename) )
str(dat)


# If a list of numbers has a distribution that is well approximated by the 
# normal distribution, what proportion of these numbers are within one 
# standard deviation away from the list's average?
pnorm(1)-pnorm(-1)

# What proportion of these numbers are within two standard deviations away
# from the list's average?
pnorm(2)-pnorm(-2)
# .. and 3 sd?
pnorm(3)-pnorm(-3)

# Define y to be the weights of males on the control diet. What proportion
# of the mice are within one standard deviation away from the average weight
# (remember to use popsd for the population sd)?
library(dplyr)
library(rafalib)
y<-dat %>% filter(Diet=="chow" & Sex=="M") %>% select(Bodyweight) %>% unlist()
mean(y<=(mean(y)+popsd(y)))-mean(y<=(mean(y)-popsd(y)))
# better...
z <- ( y - mean(y) ) / popsd(y)
mean( abs(z) <=1 )
mean( abs(z) <=2 )
mean( abs(z) <=3 )

## qq plots
qqnorm(z)
abline(0,1) # intercept 0, slope 1


# m,f, chow,hf => 4 sub populations
mypar(2,2)
y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="M" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)

y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
set.seed(1)
z<-replicate(10000,mean(sample(y,25)))
mypar(2,2)
hist(z)
z<-( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)

#..or
set.seed(1)
avgs <- replicate(10000, mean( sample(y, 25)))
#mypar(2,2)
hist(avgs)
qqnorm(avgs)
qqline(avgs)

mean(avgs)
popsd(avgs)

## CLT AND T-DISTRIBUTION IN PRACTICE EXERCISES #1
library(downloader)
library(rafalib)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv"
if(!file.exists("femaleMiceWeights.csv")) download(url,destfile=filename)
dat <- read.csv(filename)

n=100
set.seed(1)
dieRolls <- replicate(10000, mean(sample(1:6, n, replace=TRUE)==6))
p<-1/6
var<-p*(1-p)/n
z<-(dieRolls - p) / sqrt(p*(1-p)/n)
mypar(2,2)
hist(z)
qqnorm(z)
qqline(z)
mean(abs(z)>2)

# model answer
set.seed(1)
n <- 30
sides <- 2
p <- 1/sides
zs <- replicate(10000,{
    x <- sample(1:sides,n,replace=TRUE)
    (mean(x==6) - p) / sqrt(p*(1-p)/n)
}) 
qqnorm(zs)
abline(0,1)#confirm it's well approximated with normal distribution
mean(abs(zs) > 2)

tryDie<-function(n,sides){
    set.seed(1)
    p <- 1/sides
    zs <- replicate(10000,{
        x <- sample(1:sides,n,replace=TRUE)
        (mean(x==1) - p) / sqrt(p*(1-p)/n)
    })
    hist(zs)
    qqnorm(zs)
    abline(0,1)#confirm it's well approximated with normal distribution
    mean(abs(zs) > 2)
}
mypar(4,2)
tryDie(5,2) # p=0.5
tryDie(30,2) # p=0.5
tryDie(30,100) # p=0.01
tryDie(100,100) # p=0.01

# model answer
ps <- c(0.5,0.5,0.01,0.01)
ns <- c(5,30,30,100)
library(rafalib)
mypar(4,2)
for(i in 1:4){
    p <- ps[i]
    sides <- 1/p
    n <- ns[i]
    zs <- replicate(10000,{
        x <- sample(1:sides,n,replace=TRUE)
        (mean(x==1) - p) / sqrt(p*(1-p)/n)
    }) 
    hist(zs,nclass=7)
    qqnorm(zs)
    abline(0,1)
}

library(dplyr)
X <- filter(dat, Diet=="chow") %>% select(Bodyweight) %>% unlist
Y <- filter(dat, Diet=="hf") %>% select(Bodyweight) %>% unlist

mean(X)
popsd(X) # not this one, since we dont know the population
sd(X) # this one, since we have a sample.

# Use the CLT to approximate the probability that our estimate X¯ is off 
# by more than 2 grams from ??X.

se<-sd(X)/sqrt(12)
z<-2/se
pnorm(z,lower.tail=FALSE)*2

# sd of mean(X) - mean(Y)
se<-sqrt(var(X)/12 + var(Y)/12)
se

# t stat for difference of means of X and Y
diff<-mean(Y)-mean(X)
diff
tstat<-diff/se
tstat

1 - pt(3,df=3)
1 - pt(3,df=15)
1 - pt(3,df=30)
1 - pnorm(3)

# p value if t statistic were normally distributed (mean 0, sd 1)
pnorm(2.055174,lower.tail=FALSE)*2

# actual p value
t.test(X,Y)
