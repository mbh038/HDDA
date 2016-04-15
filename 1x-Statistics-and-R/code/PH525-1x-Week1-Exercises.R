library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv" 
download(url, destfile=filename)

mice<-read.csv(filename,stringsAsFactors=FALSE)
head(mice)
names(mice)[2]

mice[12,2]

mice$Bodyweight[11]

nrow(mice)

mean(mice[13:24,2])
sapply(split(mice$Bodyweight,mice$Diet),mean)

set.seed(1)

sample(mice$Bodyweight[13:24],1)

# or
set.seed(1)
i<-sample(13:24,1)
mice$Bodyweight[i]

## Dplyr exercises
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
filename <- "msleep_ggplot2.csv"
if (!file.exists(filename)) download(url,filename)
msleep <- read.csv("msleep_ggplot2.csv")
head(msleep)
class(msleep)

library(dplyr)
primates<-filter(msleep,order=="Primates")
nrow(primates)

msleep %>%
    filter(order=="Primates") %>%
    select(sleep_total) %>%
    class

msleep %>%
    filter(order=="Primates") %>%
    select(sleep_total) %>%
    unlist %>%
    class

msleep %>%
    filter(order=="Primates") %>%
    select(sleep_total) %>%
    unlist %>%
    mean

msleep %>%
    filter(order=="Primates") %>%
    summarise(mean(sleep_total))
