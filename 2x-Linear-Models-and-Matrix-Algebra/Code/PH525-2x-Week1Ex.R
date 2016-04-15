##PH525.2x

## Week 1 exercises

#install.packages("UsingR")

library(UsingR)
data("father.son",package="UsingR")

mean(father.son$sheight)

# mean height of sons whose fathers are 71 inches to nearest inch.
library(dplyr)
father.son %>%
    filter(round(fheight,0)==71) %>%
    summarize(mean(sheight))



## MATRIX NOTATION EXERCISES

x1=1:10
X=cbind(x1,2*x1,3*x1,4*x1,5*x1)
sum(X[7,])

X=matrix(seq(1,100),20,5)

X %*% matrix(1,ncol(X) )

# Solve set of linear equationa
X <- matrix(c(1,3,2,1,-2,1,1,1,-1),3,3)
y <- matrix(c(6,2,1),3,1)
solve(X)%*%y #equivalent to solve(X,y)

X<-matrix(c(3,2,1,5,4,2,-1,0,-5,2,5,0,1,-1,-5,1),4,4)
Y <- matrix(c(10,5,7,4),4,1)
solve(X)%*%Y

a <- matrix(1:12, nrow=4)
b <- matrix(1:15, nrow=3)

# matrix multiplication
c<-a%*%b
c[3,2]

# element wise multiplication
sum(a[3,]*b[,2])
