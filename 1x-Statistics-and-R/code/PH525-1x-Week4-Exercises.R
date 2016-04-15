## PH525-1x-Week4-Exercises.R

## Exploratory Data Analysis

## QQ PLOTS AND HISTOGRAMS - identifying skew

load("skew.RData")
dim(dat)

# qq plot on each columns, to compare with normal distribution
par(mfrow = c(3,3))

for (i in 1:9) {
  qqnorm(dat[,i])
}

for (i in 1:9) {
  hist(dat[,i])
}

par(mfrow=c(1,1))

## BOX PLOT EXERCISES

head(InsectSprays)
boxplot(split(InsectSprays$count, InsectSprays$spray))
boxplot(InsectSprays$count~InsectSprays$spray)

library(dplyr)
data(nym.2002, package="UsingR")

str(nym.2002)

summary(nym.2002)
nym.2002 %>% group_by(gender) %>% summarise(mean(time))
boxplot(nym.2002$time~nym.2002$gender)
hist(nym.2002$time~nym.2002$gender)

# model answer code
library(rafalib)
mypar(1,3)
males <- filter(nym.2002, gender=="Male") %>% select(time) %>% unlist
females <- filter(nym.2002, gender=="Female") %>% select(time) %>% unlist
boxplot(females, males)
hist(females,xlim=c(range( nym.2002$time)))
hist(males,xlim=c(range( nym.2002$time)))

## SCATTERPLOT EXERCISES 

library(dplyr)
data(nym.2002, package="UsingR")
males <- filter(nym.2002, gender=="Male")
females <- filter(nym.2002, gender=="Female")

cor(males$age,males$time)
cor(females$age,females$time)

library(rafalib)
mypar(2,2)
plot(females$age, females$time)
plot(males$age, males$time)
group <- floor(females$age/5) * 5
boxplot(females$time~group)
group <- floor(males$age/5) * 5
boxplot(males$time~group)


## SYMMETRY OF LOG RATIOS

time = sort(nym.2002$time)

min(time)/median(time)
max(time)/median(time)

plot(time/median(time), ylim=c(1/4,4))
abline(h=c(1/2,1,2))

plot(log2(time/median(time)),ylim=c(-2,2))
abline(h=-1:1)

## ROBUST SUMMARIES

data(ChickWeight)
head(ChickWeight)
plot( ChickWeight$Time, ChickWeight$weight, col=ChickWeight$Diet)

chick = reshape(ChickWeight, idvar=c("Chick","Diet"), timevar="Time",
                direction="wide")
head(chick)
chick = na.omit(chick)

# add an outlier to day 4 and see its effect:
summary(chick$weight.4)
summary(c(chick$weight.4,3000))

mean(c(chick$weight.4,3000))/mean(chick$weight.4) # is NOT robust to outliers
median(c(chick$weight.4,3000))/median(chick$weight.4) # IS robust to outliers
sd(c(chick$weight.4,3000))/sd(chick$weight.4) # is NOT robust to outliers
mad(c(chick$weight.4,3000))/mad(chick$weight.4) # IS robust to outliers

mypar(1,1)
plot(chick$weight.4,chick$weight.21)

# calculate Pearson correlation
pearson_No_Outlier<-cor(chick$weight.4,chick$weight.21,method="pearson")
# now with outlier
pearson_With_Outlier<-cor(c(chick$weight.4,3000),c(chick$weight.21,3000),method="pearson")
pearson_With_Outlier/pearson_No_Outlier

# calculate Spearman correlation
spearman_No_Outlier<-cor(chick$weight.4,chick$weight.21,method="spearman")
# now with outlier
spearman_With_Outlier<-cor(c(chick$weight.4,3000),c(chick$weight.21,3000),method="spearman")
spearman_With_Outlier/spearman_No_Outlier

## MANN-WHITNEY-WILCOXON TEST

data(ChickWeight)
head(ChickWeight)
plot( ChickWeight$Time, ChickWeight$weight, col=ChickWeight$Diet)
chick = reshape(ChickWeight, idvar=c("Chick","Diet"), timevar="Time",direction="wide")

head(chick)
chick = na.omit(chick)

# Mann-Whitney robust to outliers, while t-test is not.
library(dplyr)
x<-chick %>% filter(Diet==1) %>% select(weight.4) %>%unlist
y<-chick %>% filter(Diet==4) %>% select(weight.4) %>%unlist
t.test(x,y)$p.value
wilcox.test(x,y)$p.value
t.test(c(x,200),y)$p.value
wilcox.test(c(x,200),y)$p.value

# possible downside to Mann-Whitney
library(rafalib)
mypar(1,3)
boxplot(x,y)
boxplot(x,y+10)
boxplot(x,y+100)

t.test(x,y+10)$statistic-t.test(x,y+100)$statistic
wilcox.test(x,y+10)$statistic-wilcox.test(x,y+100)$statistic

wilcox.test(c(1,2,3),c(4,5,6))$p.value
wilcox.test(c(1,2,3),c(400,500,600))$p.value
