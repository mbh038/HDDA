## PH525-3x

## Statistical Inference and Modeling for High-throughput Experiments

## Week 1 - Refersher in R

library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset) ##this loads the three tables

# How many samples where processed on 2005-06-27?

head(sampleInfo)
theDate="2005-06-27"
thisDate<-sampleInfo[sampleInfo$date==theDate,]
nrow(thisDate)

# model answer
sum(sampleInfo$date=="2005-06-27")

# How many of the genes represented in this particular technology are on chromosome Y?
sum(geneAnnotation$CHR=="chrY",na.rm=TRUE)

# What is the log expression value of the for gene ARPC1A on the one subject
# that we measured on 2005-06-10?

i = which(geneAnnotation$SYMBOL=="ARPC1A")
j = which(sampleInfo$date=="2005-06-10")
geneExpression[i,j]


# Use the function apply to find the median value of each column.
# What is the median value of these values?

colMedians=apply(geneExpression,2,median) # "2" because we act on each column
median(colMedians)

# creating functions to use with apply

# Write a function that takes a vector of values e and a binary vector group
# coding two groups, and returns the p-value from a t-test: 
# t.test( e[group==1], e[group==0])$p.value.

g <- factor(sampleInfo$group)
myttest <- function(e){
    x <- e[g==1]
    y <- e[g==0]
    return( t.test(x,y)$p.value )
}
pVals=apply(geneExpression,1,myttest) # "1" because we act on each row
min(pVals)
hist(pVals)

# INFERENCE IN PRACTICE
set.seed(1)
library(downloader)
url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population = read.csv(filename)
pvals <- replicate(1000,{
    control = sample(population[,1],12)
    treatment = sample(population[,1],12)
    t.test(treatment,control)$p.val
})
head(pvals)
hist(pvals)

# proportion of p-values < 0.05
mean(pvals<0.05)
# proportion of p-values < 0.01
mean(pvals<0.01)

cases = rnorm(10,30,2)
controls = rnorm(10,30,2)
t.test(cases,controls)

set.seed(100)
pvals <- replicate(2000,{
    cases = rnorm(10,30,2)
    controls = rnorm(10,30,2)
    t.test(cases,controls)$p.val
})
head(pvals)
hist(pvals) # get a uniform distribution of p-values
x<-seq(0,1,0.1)
y<-sapply(x,function(x){mean(pvals<x)})
plot(x,y)
abline(0,1,col="red")



set.seed(100)
B=1000
sigp<-replicate(B,{
    pvals <- replicate(20,{
        cases = rnorm(10,30,2)
        controls = rnorm(10,30,2)
        t.test(cases,controls)$p.val
    })
    sum(pvals<=0.05)
})
head(sigp)
table(sigp) ##just for illustration
mean(plessthan)
hist(sigp)
mean(sigp)

#reject null hypothesis at least once
mean(sigp>0)

