##PH525.2x

## Week 3 exercises

#install.packages("UsingR")

nx<-5
ny<-7
X = cbind(rep(1,nx + ny),rep(c(0,1),c(nx, ny)))

crossprod(X) # same as t(x) %*% x



