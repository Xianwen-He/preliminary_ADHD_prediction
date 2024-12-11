### This script aims to process the data required for the modeling based on connectome data.
### @Xianwen He, Yu Chen, Dec 7, 2024

library(dplyr)

### load full data set and processed traits
alldata <- readRDS("./Data/alldata.rds")
load('./Data/demographic.data.RData')

### Variable selection
# load PCA data
connectome.PCA.dat <- alldata %>% select("id", "fc_matrix", "sc_matrix",
                                     paste0("s", 1:60), paste0("f", 1:60))
# merge with traits
connectome.dat <- merge(connectome.PCA.dat, demographic.data, by='id', all=T)

# attach fcp and scp all
connectome.dat$fcpcall <- alldata$fcpcall
connectome.dat$scpcall <- alldata$scpcall

# save the data
saveRDS(connectome.dat, file = "./Data/connectome.dat.rds")
