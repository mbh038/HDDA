##PH525.2x

## Week 4 exercises

## QR Decomposition


list.of.packages <- c("contrast","rafalib")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(contrast)
library(rafalib)


url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)

fit <- lm(friction ~ type + leg, data=spider)
betahat <- coef(fit)

Y <- matrix(spider$friction, ncol=1)
X <- model.matrix(~ type + leg, data=spider)

# Solution for beta is

# R beta = Q^T Y so beta = R^-1 Q^T Y

QR <- qr(X)
Q <- qr.Q( QR )
R <- qr.R( QR )
(betahat <- backsolve(R, crossprod(Q,Y) ) )

solve(R)%*%crossprod(Q,Y)
