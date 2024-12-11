### This script is compiled to build up an XGBoost model
### for the prediction on DSM_Adh_Raw_log
### based on the traits in interest
### @Xianwen He, Nov 25, 2025


### preparation

# import libraries
library(xgboost)
library(ggplot2)

# customized theme
mytheme <- theme(
  panel.background = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  axis.title = element_text(size = 14, hjust=0.5, family='serif'),
  axis.text = element_text(size = 12, family='serif'),  
  plot.title = element_text( # Center the title
    hjust = 0.5, size = 16, face = "bold", family='serif')
)
myblue <- 'cornflowerblue'

# function for data split
prepare.dat.func <- function(dat, proportion, seed=0){
  # Split the data set and convert the subsets into matrices
  # The label should be recorded in the last column
  
  # split the data set
  set.seed(seed)
  train.idx <- sample(1:nrow(dat), size = proportion * nrow(dat))
  train.dat <- dat[train.idx, ]
  test.dat <- dat[-train.idx, ]
  
  # convert into matrices
  train.mat <- as.matrix(train.dat[, -ncol(train.dat)])
  train.label <- train.dat[, ncol(train.dat)]
  test.mat <- as.matrix(test.dat[, -ncol(test.dat)])
  test.label <- test.dat[, ncol(test.dat)]
  
  return (list('train.mat'=train.mat, 
               'train.label'=train.label,
               'test.mat'=test.mat,
               'test.label'=test.label))
}



### process data

# load data
load('./Data/traits.log.dat.RData')

# convert to numeric
mod.dat <- traits.log.dat[, -1]
mod.dat.numeric <- data.frame(lapply(mod.dat, function(x){
  if (is.factor(x)){
    return (as.numeric(x))
  }else{
    return (x)
  }
}))

# prepare the data set
dat.lst <- prepare.dat.func(mod.dat.numeric, 0.8)

# convert to DMatrix
dtrain <- xgb.DMatrix(data = dat.lst$train.mat, label = dat.lst$train.label)
dtest <- xgb.DMatrix(data = dat.lst$test.mat, label = dat.lst$test.label)



### model training
# set seed to control randomness
set.seed(0)
model.depth6 <- xgboost(data = dtrain, nrounds = 50, objective = "reg:squarederror",
                        eta = 0.1, max_depth = 6)  # set the depth as 6


### prediction
preds <- predict(model.depth6, dtest)
mse <- mean((preds - dat.lst$test.label)^2)
cat("RMSE:", sqrt(mse), "\n")
cat('Variance explained:', 1.0-mse/mean((dat.lst$test.label-mean(dat.lst$test.label))^2), "\n")

# feature importance
importance <- xgb.importance(feature_names = colnames(dat.lst$train.mat), model = model.depth6)
importance.plot.depth6 <- xgb.ggplot.importance(importance) + mytheme + theme(legend.position = 'none')
ggsave('./Data/img/xgb_depth6.png', importance.plot.depth6, width=10, height=12, dpi=200)

### tune the parameters
set.seed(0)
model.depth10 <- xgboost(data = dtrain, nrounds = 50, objective = "reg:squarederror",
                         eta = 0.1, max_depth = 10)  # set depth as 10
# prediction
preds <- predict(model.depth10, dtest)
mse <- mean((preds - dat.lst$test.label)^2)
cat("RMSE:", sqrt(mse), "\n")
cat('Variance explained:', 1.0-mse/mean((dat.lst$test.label-mean(dat.lst$test.label))^2), "\n")
# model summary
importance <- xgb.importance(feature_names = colnames(dat.lst$train.mat), model = model.depth10)
importance.plot.depth10 <- xgb.ggplot.importance(importance) + mytheme + theme(legend.position = 'none')
ggsave('./Data/img/xgb_depth10.png', importance.plot.depth10, width=10, height=12, dpi=200)

### save model checkpoints
save(model.depth6, file='./Data/Checkpoints/model.depth6.RData')
save(model.depth10, file='./Data/Checkpoints/model.depth10.RData')
