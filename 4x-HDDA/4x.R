## Distance Exercises

library(devtools)
install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e)
head(tissue)

# How many biological replicates for hippocampus?

str(tissue)
str(e)

table(tissue)

# What is the distance between samples 3 and 45?

x<-e[,3]
y<-e[,45]

# a
sqrt(sum((x-y)^2))

# b
sqrt(crossprod(x-y))

# c
d<-dist(t(e))
as.matrix(d)[3,45]
as.matrix(d)[45,3]

# What is the distance between gene 210486_at and 200805_at

x<-e["210486_at",]
y<-e["200805_at",]
sqrt(crossprod(x-y))

# If I run the command (don't run it!)
# d = as.matrix(dist( e))
# how many cells does d have?

nrow(e)^2

# How many distances are stored in d? (Hint: What is the length of d)?
d = dist(t(e))
num<-(nrow(as.matrix(d))^2-nrow(as.matrix(d)))/2
num
length(d) # does the same!

