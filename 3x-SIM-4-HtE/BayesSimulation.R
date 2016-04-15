# Bayes Simulation
# from Irizarry and Love

set.seed(3)
prev <- 1/20
##Later, we are arranging 1000 people in 80 rows and 20 columns
M <- 50 ; N <- 30
##do they have the disease?
d<-rbinom(N*M,1,p=prev)

accuracy <- 0.9
test <- rep(NA,N*M)
##do controls test positive?
test[d==1]  <- rbinom(sum(d==1), 1, p=accuracy)
##do cases test positive?
test[d==0]  <- rbinom(sum(d==0), 1, p=1-accuracy)

cols <- c("grey","red")
people <- expand.grid(1:M,N:1)
allcols <- cols[d+1] ##Cases will be red

positivecols <- allcols
positivecols[test==0] <- NA ##remove non-positives

negativecols <- allcols
negativecols[test==1] <- NA ##remove non-positives

library(rafalib)
mypar()
layout(matrix(c(1,2,1,3),2,2),width=c(0.35,0.65))
###plot of all people
plot(people,col=allcols,pch=16,xaxt="n",yaxt="n",xlab="",ylab="",
     main=paste0("Population: ",round(mean(d)*100),"% are red"))

plot(people,col=positivecols,pch=16,xaxt="n",yaxt="n",xlab="",ylab="",
     main=paste("Tested Positive:",round(mean(d[test==1])*100),"% are red"))

plot(people,col=negativecols,pch=16,xaxt="n",yaxt="n",xlab="",ylab="",
     main=paste("Tested Negative:",round(mean(d[test==0])*100,1),"% are red"))