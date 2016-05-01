# In the blood composition assessment we will be downloading a dataset from GEO directly from R. 
# This can take several minutes depending on the download speed. If you want to start downloading
# this file ahead of time here are the lines of code:
  
library(minfi)
grset=getGenomicRatioSetFromGEO("GSE32148")
# If you have a slow internet connection you might want to save this file to disk for future upload.

##after download  you can save the object like this:
save(grset,file="grset.rda")
##then to load later
load("grset.rda")