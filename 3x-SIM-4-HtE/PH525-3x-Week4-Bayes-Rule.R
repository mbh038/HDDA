## PH525-3x

## Statistical Inference and Modelling for High-throughput Experiments

## Week 4 - Bayes Rule, Hierarchical Modelling


## BAYES RULE IN PRACTICE

# download baseball data
tmpfile <- tempfile()
tmpdir <- tempdir()
download.file("http://seanlahman.com/files/database/lahman-csv_2014-02-14.zip",tmpfile)
##this shows us files
filenames <- unzip(tmpfile,list=TRUE)
players <- read.csv(unzip(tmpfile,files="Batting.csv",exdir=tmpdir),as.is=TRUE)
unlink(tmpdir)
file.remove(tmpfile)

str(players)

# Find the batting averages (AVG) for players with more than 500 at bats (AB)
# in 2012:

library(dplyr)
filter(players,yearID==2012) %>%
    mutate(AVG=H/AB) %>%
    filter(AB>=500) %>%
    select(AVG)

# Edit the command above to obtain all the batting averages from 2010, 2011,
# 2012 and removing rows with AB < 500.

dat<-filter(players,yearID>=2010,yearID<=2012) %>%
    mutate(AVG=H/AB) %>%
    filter(AB>=500) %>%
    select(AVG)
mean(dat$AVG)

# What is the standard deviation of these batting averages?

sd(dat$AVG)

# Use EDA to identify best approximation ot odistribution of batting averages

hist(dat$AVG,nc=100,freq=FALSE)
# looks normal
mypar(1,1)
#check
qqnorm(dat$AVG)
qqline(dat$AVG)
