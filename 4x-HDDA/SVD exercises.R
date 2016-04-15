## PH525 4x HDDA

## SVD

m <- 100
n <- 2
x <- rnorm(m)
e <- rnorm(n*m,0,0.01)
Y <- cbind(x,x)+e
cor(Y)

str(Y)
str(e)


# SVD exercises

library(tissuesGeneExpression)
data(tissuesGeneExpression)

s = svd(e)
signflips = sample(c(-1,1),ncol(e),replace=TRUE)
signflips

#Now we switch the sign of each column and check that we get the same answer.
# We do this using the function sweep. If x is a matrix and a is a vector, then
# sweep(x,1,y,FUN="*") applies the fun FUN to each row i FUN(x[i,],a[i]), in this
# case x[i,]*a[i]. If instead of 1 we use 2 sweep applies this to columns.
# To learn about sweep read ?sweep. 

newu= sweep(s$u,2,signflips,FUN="*")
newv= sweep(s$v,2,signflips,FUN="*" )
all.equal( s$u %*% diag(s$d) %*% t(s$v), newu %*% diag(s$d) %*% t(newv))

# SVD Exercises #1
s = svd(e)
m = rowMeans(e)

# What is the correlation between the first column of U and m?
cor(m,s$u[,1])


# SVD Exercises #2

newmeans = rnorm(nrow(e)) ##random values we will add to create new means
newe = e+newmeans ##we change the means
# changing the means does not change the distances:
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(newe[,3]-newe[,45])) 

# So we might as well make the mean of each row 0 since it does not help us
# approximate the column distances. We will define y as the detrended e and
# recompute the SVD:

y = e - rowMeans(e)
s = svd(y)

# We showed that UDV??? is equal to y up to numerical error
resid = y - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))

# The above can be made more efficient in two ways. First, using the crossprod
# and second not creating a diagonal matrix. Note that in R we can multiply a 
# matrix x by vector a. The result is a matrix with row i equal to x[i,]*a[i].
# Here is an example to illustrate this.

x=matrix(rep(c(1,2),each=5),5,2)
x

x*c(1:5)

# this code is equivalent to

sweep(x,1,1:5,"*")

# This means that we don't have to convert s$d into a matrix to obtain DV???.

a<- diag(s$d)%*%t(s$v)
b <- s$d * t(s$v)
identical(a,b)

# SVD Exercises #3

# If we define vd = t(s$d * t(s$v)) then which of the following is not the same UDV??? :
udvt= s$u %*% (s$d * t(s$v))
vd = t(s$d * t(s$v))
b<-tcrossprod(s$u,vd)
identical(udvt,b)


# SVD Exercises #4

z = s$d * t(s$v)

sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(y[,3]-y[,45]))
sqrt(crossprod(z[,3]-z[,45]))

approxdist<-sqrt(crossprod(z[1:2,3]-z[1:2,45]))
truedist<-sqrt(crossprod(e[,3]-e[,45]))
abs(truedist-approxdist)

# SVD Exercises #5

ks = 1:189
realdistance = sqrt(crossprod(e[,3]-e[,45]))
approxdistances = sapply(ks,function(k){
    sqrt(crossprod(z[1:k,3,drop=FALSE]-z[1:k,45,drop=FALSE] )) 
})
percentdiff = 100*abs(approxdistances - realdistance)/realdistance
plot(ks,percentdiff) ##take a look
min(ks[which(percentdiff < 10)])

# SVD Exercises #6

distances = sqrt(apply(e[,-3]-e[,3],2,crossprod))
approxdistances= sqrt(apply(z[1:2,-3]-z[1:2,3],2,crossprod))
plot(distances,approxdistances) ##take a look
cor(distances,approxdistances,method="spearman")

# this shows how just two dimensions can be useful to get an idea about the actual distances.
