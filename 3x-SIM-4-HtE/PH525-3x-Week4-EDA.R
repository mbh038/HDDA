## PH525-3x

## Statistical Inference and Modelling for High-throughput Experiments

## Week 4 - EDA


source("http://www.bioconductor.org/biocLite.R")
biocLite("SpikeInSubset")
library(SpikeInSubset)
data(mas133)

library(Biobase)
e <- exprs(mas133)
plot(e[,1],e[,2],main=paste0("corr=",signif(cor(e[,1],e[,2]),3)),cex=0.5)
k <- 3000
b <- 1000 #a buffer
polygon(c(-b,k,k,-b),c(-b,-b,k,k),col="red",density=0,border="red")

# What proportion of the points are inside the box?

edf<-as.data.frame(e)
all<-length(edf[,1])
library(dplyr)
filter(edf,edf[,1]<k,edf[,2]<k) %>%
    nrow()/all

# better way to find this!
mean(e[,1]<k & e[,2]<k)

# now make sample plot with log
plot(log2(e[,1]),log2(e[,2]),main=paste0("corr=",signif(cor(log2(e[,1]),log2(e[,2])),2)),cex=0.5)
k <- log2(3000)
b <- log2(0.5)
polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")
# ensures that 95% of data no longer in tiny section  in corner of plot

# Make an MA-plot
e <- log2(exprs(mas133))
plot((e[,1]+e[,2])/2,e[,2]-e[,1],cex=0.5)

# What is the standard deviation of the log ratios for this comparison?
sd((e[,1]+e[,2])/(e[,2]-e[,1]))
sd(e[,1]/e[,2])

logratio <- (e[,1]-e[,2])
sum(abs(logratio)>1)

