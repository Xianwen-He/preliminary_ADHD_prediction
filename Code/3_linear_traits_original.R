### This script is compiled to build up a linear regression model
### for the prediction on DSM_Adh_Raw
### based on the traits in interest
### @Xianwen He, Nov 25, 2025

library(dplyr)
library(ggplot2)

### preparation

# function for data split
prepare.dataframe.func <- function(dat, proportion, seed=0){
  # split the data set
  set.seed(seed)  # randomness control
  train.idx <- sample(1:nrow(dat), size = proportion * nrow(dat))
  train.dat <- dat[train.idx, ]
  test.dat <- dat[-train.idx, ]
  
  return (list('train.dat'=train.dat, 
               'test.dat'=test.dat))
}



### process data

# load data
load('./Data/traits.dat.RData')

# remove ids
mod.dat <- traits.dat[, -1]

# split the data set
dat.lst <- prepare.dataframe.func(mod.dat, 0.8)


### model training
model.linear.original <- lm(DSM_Adh_Raw~., data = dat.lst$train.dat)
summary(model.linear.original)


### prediction
preds <- predict.lm(model.linear.original, newdata = dat.lst$test.dat)
train.mse <- mean((model.linear.original$residuals)^2)
mse <- mean((preds - dat.lst$test.dat$DSM_Adh_Raw)^2)
mst.test <- mean((dat.lst$test.dat$DSM_Adh_Raw - mean(dat.lst$test.dat$DSM_Adh_Raw))^2)
cat("Training RMSE", sqrt(train.mse), '\n')
cat("RMSE:", sqrt(mse), "\n")
cat('Variance explained:', 1.0-mse/mst.test, "\n")


### save model checkpoints
save(model.linear.original, file='./Data/Checkpoints/model.linear.original.RData')
