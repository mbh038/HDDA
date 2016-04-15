## PH525-3x

## Statistical Inference and Modelling for High-throughput Experiments

## Week 4 - Hierarchical Modelling


## HIERARCHICAL MODELS IN PRACTICE

## to install:
library(rafalib)
install_bioc("SpikeInSubset")
library(Biobase)
library(SpikeInSubset)
data(rma95)
y <- exprs(rma95)

# What proportion of genes with a p-value < 0.01 (no multiple comparison correction)
# are not part of the artificially added (false positive)?
# sum(!spike & smallp)/sum(smallp) # works!

library(genefilter)
g <- factor(rep(0:1,each=3))

# create index of which rows are associated with the artificially added genes:
spike <- rownames(y) %in% colnames(pData(rma95))

rtt = rowttests(y,g)
index = rtt$p.value < 0.01 
print (mean( !spike[index] ))# is the answer

## We can make a volcano plot to visualize this:
mask <- with(rtt, abs(dm) < .2 & p.value < .01)
cols <- ifelse(mask,"red",ifelse(spike,"dodgerblue","black"))
with(rtt,plot(-dm, -log10(p.value), cex=.8, pch=16,
              xlim=c(-1,1), ylim=c(0,5),
              xlab="difference in means",
              col=cols))
abline(h=2,v=c(-.2,.2), lty=2)

library(dplyr)
g <- factor(rep(0:1,each=3))
spike <- rownames(y) %in% colnames(pData(rma95))
rtt = rowttests(y,g)
index = rtt$p.value < 0.01

wg1sds<-data.frame("sds"=rowSds(y[,1:3]),spike,index)

mypar(1,4)
boxplot(sds~spike,data=wg1sds)
boxplot(sds~index,data=wg1sds)

wg1sds$cat<-ifelse(wg1sds$spike & wg1sds$index,"TP",ifelse(!wg1sds$spike & wg1sds$index,"FP",ifelse(wg1sds$spike & !wg1sds$index,"FN","TN")))

table(wg1sds$cat)
mypar(1,1)
boxplot(sds~cat,data=wg1sds)

# model answer
library(genefilter)
sds <- rowSds(y[,g==0])
index <- paste0( as.numeric(spike), as.numeric(rtt$p.value<0.01))
index <- factor(index,levels=c("11","01","00","10"),labels=c("TP","FP","TN","FN"))
boxplot(split(sds,index))

library(limma)
fit <- lmFit(y, design=model.matrix(~ g))
colnames(coef(fit))
fit <- eBayes(fit)

sampleSD = fit$sigma
posteriorSD = sqrt(fit$s2.post)

LIM = range( c(posteriorSD,sampleSD))
plot(sampleSD, posteriorSD,ylim=LIM,xlim=LIM)
abline(0,1)
abline(v=sqrt(fit$s2.prior))

library(limma)
fit = lmFit(y, design=model.matrix(~ g))
fit = eBayes(fit)
##second coefficient relates to diffences between group
pvals = fit$p.value[,2] 
newindex = pvals < 0.01
print (mean( !spike[newindex] ))# is the answer

## We can make a volcano plot to visualize this:
mask <- abs(fit$coef[,2]) < .2 & fit$p.value[,2] < .01
cols <- ifelse(mask,"red",ifelse(spike,"dodgerblue","black"))
plot(fit$coef[,2], -log10(fit$p.value[,2]), cex=.8, pch=16,
     xlim=c(-1,1), ylim=c(0,5),
     xlab="difference in means",
     col=cols)
abline(h=2,v=c(-.2,.2), lty=2)