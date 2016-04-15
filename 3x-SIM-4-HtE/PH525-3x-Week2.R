## PH525-3x

## Statistical Inference and Modeling for High-throughput Experiments

## Week 2 - Vectorizing Code

## ERROR RATES AND PROCEDURES EXERCISES

library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
set.seed(1)
population = unlist( read.csv("femaleControlsPopulation.csv") )

alpha <- 0.05
N <- 12
m <- 10000
p0 <- 0.90 ##10% of diets work, 90% don't
m0 <- m*p0
m1 <- m-m0
nullHypothesis <- c( rep(TRUE,m0), rep(FALSE,m1))
delta <- 3

B <- 10 ##number of simulations 
system.time(
    VandS <- replicate(B,{
        calls <- sapply(1:m, function(i){
            control <- sample(population,N)
            treatment <- sample(population,N)
            if(!nullHypothesis[i]) treatment <- treatment + delta
            t.test(treatment,control)$p.val < alpha
        })
        c(sum(nullHypothesis & calls),sum(!nullHypothesis & calls))
    })
)


# now vectorize this - runs much faster
library(rafalib)
install_bioc("genefilter")
library(genefilter) ##rowttests is here
set.seed(1)
##Define groups to be used with rowttests
g <- factor( c(rep(0,N),rep(1,N)) )
B <- 10 ##number of simulations
system.time(
    VandS <- replicate(B,{
        ##matrix with control data (rows are tests, columns are mice)
        controls <- matrix(sample(population, N*m, replace=TRUE),nrow=m)
        
        ##matrix with control data (rows are tests, columns are mice)
        treatments <-  matrix(sample(population, N*m, replace=TRUE),nrow=m)
        
        ##add effect to 10% of them
        treatments[which(!nullHypothesis),]<-treatments[which(!nullHypothesis),]+delta
        
        ##combine to form one matrix
        dat <- cbind(controls,treatments)
        
        calls <- rowttests(dat,g)$p.value < alpha
        
        c(sum(nullHypothesis & calls),sum(!nullHypothesis & calls))
    })
)

library(genefilter) ##rowttests is here
set.seed(1)
##Define groups to be used with rowttests
g <- factor( c(rep(0,N),rep(1,N)) )
B <- 1000 ##number of simulations
Qs <- replicate(B,{
  ##matrix with control data (rows are tests, columns are mice)
  controls <- matrix(sample(population, N*m, replace=TRUE),nrow=m)
  
  ##matrix with control data (rows are tests, columns are mice)
  treatments <-  matrix(sample(population, N*m, replace=TRUE),nrow=m)
  
  ##add effect to 10% of them
  treatments[which(!nullHypothesis),]<-treatments[which(!nullHypothesis),]+delta
  
  ##combine to form one matrix
  dat <- cbind(controls,treatments)
  
  calls <- rowttests(dat,g)$p.value < alpha
  R=sum(calls)
  Q=ifelse(R>0,sum(nullHypothesis & calls)/R,0)
  return(Q)
})


# BONFERRONI CORRECTION EXERCISES #1 (BONFERONNI VERSUS SIDAK)

# model answer
alphas <- seq(0,0.25,0.01)
par(mfrow=c(2,2))
for(m in c(2,10,100,1000)){
    plot(alphas,alphas/m - (1-(1-alphas)^(1/m)),type="l")
    abline(h=0,col=2,lty=2)
}
par(mfrow=c(1,1))
pvals <- runif(8793,0,1)
hist(pvals)


# Bonferroni FWER
set.seed(1)
B <- 10000
N <- 8793
bfc <- .05 / N

sigs <- replicate(B, {
    pvals <- runif(N, 0, 1)
    sum(pvals < bfc)
})
mean(sigs > 0)

# matrix way of doing this
set.seed(1)
B <- 10000
m <- 8793
alpha <- 0.05
pvals <- matrix(runif(B*m,0,1),B,m)
k <- alpha/m
mistakes <- rowSums(pvals<k) 
mean(mistakes>0)

# Sidak FWER
set.seed(1)
B <- 10000
N <- 8793
sdkc <- 1-(1-alpha)^(1/m)

sigs <- replicate(B, {
    pvals <- runif(N, 0, 1)
    sum(pvals < sdkc)
})
mean(sigs > 0)

## FDR EXERCISES

library(devtools)
library(rafalib)
# install_github("genomicsclass/GSE5859Subset")
load("./data/GSE5859Subset.rda") # cloned from github
data(GSE5859Subset)
install_bioc("genefilter")
install_bioc("qvalue")

# Load gene expression data

#library(GSE5859Subset)


# We are interested in comparing gene expression between the two groups defined 
# in the sampleInfo table.

# Compute a p-value for each gene using the function rowttests from the genefilter
# package in Bioconductor.
library(genefilter)
?rowttests

alpha<-0.05
g <- factor(sampleInfo$group)
pvals<-rowttests(geneExpression,g)$p.value
sum(pvals<alpha)

# Apply the Bonferroni correction to the p-values obtained in question #1 to achieve
# a FWER of 0.05. How many genes are called significant under this procedure?

m<-nrow(pvals)
k <- alpha / m
sum(pvals<k)
# gives us list of 10 genes that we are quite certain has no false positives

# q values and the FDR (False Discovery Rate)
?p.adjust
g = factor(sampleInfo$group)
pvals = rowttests(geneExpression,g)$p.value
fdr = p.adjust(pvals,method="fdr")
sum(fdr<alpha)
# gives us list of 13 genes of which we think about 5% are false positives

# Now use the qvalue function, in the Bioconductor qvalue package,
# to estimate q-values using the procedure described by Storey.

library(qvalue)
?qvalue
g = factor(sampleInfo$group)
pvals = rowttests(geneExpression,g)$p.value
res<-qvalue(pvals)
sum(res$qvalues<alpha)
# gives us a list of 22 genes

# Read the help file for qvalue and report the estimated proportion of genes for which
# the null hypothesis is true ??0=m0/m
res$pi0

# Note that we have the number of genes passing the q-value <0.05 threshold is larger
# with the qvalue function than the p.adjust difference.

# This is because the qvalue function estimates the proportion of genes for which the
# null hypothesis is true and provides a less conservative estimate

plot(qvalue(pvals)$qvalue/p.adjust(pvals,method="fdr"))
abline(h=qvalue(pvals)$pi0,col=2)

hist(pvals,breaks=seq(0,1,len=21))
expectedfreq <- length(pvals)/20 #per bin
abline(h=expectedfreq*qvalue(pvals)$pi0,col=2,lty=2)

## FDR EXERCISES #7
# Create a Monte Carlo Simulation in which you simulate measurements from 8,793 genes
# for 24 samples: 12 cases and 12 controls.
n <- 24
m <- 8793
mat <- matrix(rnorm(n*m),m,n)

# Now for 500 genes, there is a difference of 2 between cases and controls:

delta <- 2
positives <- 500
mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta

# So the null hypothesis is true for 8793-500 genes.
# Using the notation from the videos m=8793, m0=8293 and m1=500

m0<-m-positives
m1<-positives
nullHypothesis <- c( rep(FALSE,m1), rep(TRUE,m0))


# Set the seed at 1, set.seed(1) and run this experiment 1,000 times with a
# Monte Carlo simulation. For each instance compute p-values using a t-test
# (using rowttests in the genefilter package) and create three lists of genes
# using:
    
# Bonferroni correction to achieve an FWER of 0.05,
# p-adjust estimates of FDR to achieve an FDR of 0.05, and
# qvalue estimates of FDR to to achieve an FDR of 0.05.

library(genefilter)
alpha<-0.05
g <- factor( c(rep(0,n/2),rep(1,n/2)) )



B <- 1000 ##number of simulations
n <- 24
m <- 8793
delta <- 2
positives <- 500
m0<-m-positives
m1<-positives
nullHypothesis <- c( rep(FALSE,m1), rep(TRUE,m0))
library(genefilter)
alpha<-0.05
g <- factor( c(rep(0,n/2),rep(1,n/2)) )

set.seed(1)
# bonferroni
ptab <- replicate(B,{
    mat <- matrix(rnorm(n*m),m,n)
    mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
    # Bonferroni
    calls <- rowttests(mat,g)$p.value  < alpha/m # Bonferroni
    V1<-sum(calls & nullHypothesis)
    S1<-sum(calls & !nullHypothesis)
    #  qvalues from p.adjust
    pvals = rowttests(mat,g)$p.value
    calls <- p.adjust(pvals,method="fdr") < alpha # fdr
    V2<-sum(calls & nullHypothesis)
    S2<-sum(calls & !nullHypothesis)
    # qvalues from qvalue
    pvals = rowttests(mat,g)$p.value
    res<-qvalue(pvals)
    calls<-res$qvalues<alpha
    V3<-sum(calls & nullHypothesis)
    S3<-sum(calls & !nullHypothesis)    

    c(V1,S1,V2,S2,V3,S3)
})
FPR1<-sum(ptab[1,])/(B*m0)
FPR1
FNR1<-((B*m1)-sum(ptab[2,]))/(B*m1)
FNR1

FPR2<-sum(ptab[3,])/(B*m0)
FPR2
FNR2<-((B*m1)-sum(ptab[4,]))/(B*m1)
FNR2

FPR3<-sum(ptab[5,])/(B*m0)
FPR3
FNR3<-((B*m1)-sum(ptab[6,]))/(B*m1)
FNR3

