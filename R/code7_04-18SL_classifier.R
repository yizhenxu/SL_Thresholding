#load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/Risk Score Wrap Up/04-01withCVCfromISRES10000POP_goodSeed.RData")
memory.limit(size=21264)

library("SuperLearner")
library("randomForest")
library("gam")
library("rpart")
library("dplyr")
library("plyr")
library("ggplot2")
library("nloptr")
library("lpSolve")

setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up")
setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up")
source("functions.r")

########################################################################################
########################################################################################
########################################################################################

# Define SL templates for classifier

# GLM
SL.glm.classifier <- function(Y, X, newX, family, obsWeights, ...) {
  fit.glm <- glm(Y ~ ., data = X, family = family, weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type="response")
  fit <- list(object = fit.glm)
  
  # get training prediction
  train.p <- predict(fit.glm, newdata = X, type="response")
  # get optimal cutoff based on the training prediction
  opt <- Opt.nonpar.rule(Y,train.p,phi=0,lambda)
  cutoff <- as.numeric(opt)[1]
  # prediction is compared with training cutoff
  pred.c <- as.numeric(pred > cutoff)
  
  class(fit) <- "SL.glm"
  out <- list(pred = pred.c, fit = list(fit, cutoff))
  return(out)
}

predict.SL.glm.classifier <- function(object, newdata, ...){
  pred <- predict(object = object[[1]]$object, newdata = newdata, type = "response")
  pred.c <- as.numeric(pred > object[[2]])
  return(pred.c)
}

########################################################################################

# GAM
SL.gam.classifier <- function(Y, X, newX, family, obsWeights, deg.gam = 2, cts.num = 4, ...) {
  require('gam')
  if("mgcv" %in% loadedNamespaces()) warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
  # create the formula for gam with a spline for each continuous variable
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  if (sum(!cts.x) > 0) { 
    gam.model <- as.formula(paste("Y~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,")", sep=""), collapse = "+"), "+", paste(colnames(X[, !cts.x, drop=FALSE]), collapse = "+")))
  } else {
    gam.model <- as.formula(paste("Y~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, ")", sep=""), collapse = "+")))
  }
  # fix for when all variables are binomial
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(paste("Y~", paste(colnames(X), collapse = "+"), sep = ""))
  }
  fit.gam <- gam::gam(gam.model, data = X, family = family, control = gam::gam.control(maxit = 50, bf.maxit = 50), weights = obsWeights)
  pred <- gam::predict.gam(fit.gam, newdata = newX, type = "response")
  fit <- list(object = fit.gam)
  
  # get training prediction
  train.p <- gam::predict.gam(fit.gam, newdata = X, type = "response")
  # get optimal cutoff based on the training prediction
  opt <- Opt.nonpar.rule(Y,train.p,phi=0,lambda)
  cutoff <- as.numeric(opt)[1]
  # prediction is compared with training cutoff
  pred.c <- as.numeric(pred > cutoff)
  
  class(fit) <- c("SL.gam")
  out <- list(pred = pred.c, fit = list(fit, cutoff))
  return(out)
}

predict.SL.gam.classifier <- function(object, newdata, ...){
  .SL.require('gam')
  pred <- gam::predict.gam(object = object[[1]]$object, newdata = newdata, type = "response")
  pred.c <- as.numeric(pred > object[[2]])
  return(pred.c)
}

########################################################################################

# Random Forest

SL.randomForest.classifier <- function(Y, X, newX, family, mtry = ifelse(family$family=="gaussian", max(floor(ncol(X)/3), 1), floor(sqrt(ncol(X)))), ntree=1000, nodesize = ifelse(family$family=="gaussian", 5, 1), ...) {
  require('randomForest')
  if(family$family=="gaussian"){
    fit.rf <- randomForest::randomForest(Y ~ ., data = X, ntree = ntree, xtest = newX, keep.forest = TRUE, mtry = mtry, nodesize = nodesize)
    pred <- fit.rf$test$predicted
    fit <- list(object = fit.rf)
    train.p <- predict(fit.rf, newdata = X, type = 'response')
  }
  if(family$family=="binomial"){
    fit.rf <- randomForest::randomForest(y = as.factor(Y), x = X, ntree = ntree, xtest = newX, keep.forest = TRUE, mtry = mtry, nodesize = nodesize)
    pred <- fit.rf$test$votes[, 2]
    fit <- list(object = fit.rf)
    train.p <- predict(fit.rf, newdata = X, type = 'vote')[,2]
  }
  
  # get optimal cutoff based on the training prediction
  opt <- Opt.nonpar.rule(Y,train.p,phi=0,lambda)
  cutoff <- as.numeric(opt)[1]
  # prediction is compared with training cutoff
  pred.c <- as.numeric(pred > cutoff)
  
  class(fit) <- c("SL.randomForest")
  out <- list(pred = pred.c, fit = list(fit, cutoff))
  return(out)
  
}

predict.SL.randomForest.classifier <- function(object, newdata, family, ...){
  .SL.require('randomForest')
  if(family$family=="gaussian"){
    pred <- predict(object[[1]]$object, newdata = newdata, type = 'response')
  }
  if(family$family=="binomial"){
    pred <- predict(object[[1]]$object, newdata = newdata, type = 'vote')[,2]
  }
  pred.c <- as.numeric(pred > object[[2]])
  return(pred.c)
}
########################################################################################

# Rpart

SL.rpart.classifier <- function(Y, X, newX, family, obsWeights, cp = 0.01, minsplit = 20, xval = 10, maxdepth = 30, minbucket = round(minsplit/3), ...) {
  require('rpart')
  if(family$family == "gaussian"){
    fit.rpart <- rpart::rpart(Y~., data = data.frame(Y, X), control = rpart::rpart.control(cp = cp, minsplit = minsplit, xval = xval, maxdepth = maxdepth, minbucket = minbucket), method = "anova", weights = obsWeights)
    pred <- predict(fit.rpart, newdata = newX)
    train.p <- predict(fit.rpart, newdata = X)
  }
  if(family$family == "binomial") {
    fit.rpart <- rpart::rpart(Y ~ ., data = data.frame(Y, X), control = rpart::rpart.control(cp = cp, minsplit = minsplit, xval = xval, maxdepth = maxdepth, minbucket = minbucket), method = "class", weights = obsWeights)
    pred <- predict(fit.rpart, newdata = newX)[, 2]
    train.p <- predict(fit.rpart, newdata = X)[, 2]
  }
  
  # get optimal cutoff based on the training prediction
  opt <- Opt.nonpar.rule(Y,train.p,phi=0,lambda)
  cutoff <- as.numeric(opt)[1]
  # prediction is compared with training cutoff
  pred.c <- as.numeric(pred > cutoff)
  
  fit <- list(object = fit.rpart)
  class(fit) <- c("SL.rpart")
  out <- list(pred = pred.c, fit = list(fit, cutoff))
  return(out)
  
}

# 
predict.SL.rpart.classifier <- function(object, newdata, family, ...) {
  .SL.require('rpart')
  if(family$family=="gaussian") { 
    pred <- predict(object[[1]]$object, newdata = newdata)
  }
  if(family$family=="binomial") {
    pred <- predict(object[[1]]$object, newdata = newdata)[, 2]
  }
  pred.c <- as.numeric(pred > object[[2]])
  return(pred.c)
}

########################################################################################
SL.library = c("SL.randomForest.classifier","SL.glm.classifier","SL.gam.classifier", "SL.rpart.classifier")

lambda=0.5
test = SuperLearner(Y=Y, X=W,newX=W, SL.library = SL.library, family = binomial(), method = "method.NNLS",verbose = FALSE)
aaa = test$Z
lambda=0.05
test = SuperLearner(Y=Y, X=W,newX=W, SL.library = SL.library, family = binomial(), method = "method.NNLS",verbose = FALSE)
bbb = test$Z
all(aaa == bbb) # FALSE -> GREAT!

########################################################################################
########################################################################################
########################################################################################
### cross validation

### IMPORTANT NOTICE
### since lambda is coorporated in the built in SL candidate classifiers
### we need to save a CV Yhat for each lambda and calculate its corresponding FNR FPR Risk
N = length(Y)
#K=3
K = 10

N0 = table(Y)[1]
l0 = floor(N0/K)
t0 = rep(1:K,rep(l0,K)+c(rep(1,N0%%l0),rep(0,K-N0%%l0)))

N1 = table(Y)[2]
l1 = floor(N1/K)
t1 = rep(1:K,rep(l1,K)+c(rep(1,N1%%l1),rep(0,K-N1%%l1)))

data = cbind(Y,W)

set.seed(100)
ind0 = sample(1:N0,replace=FALSE)
ind1 = sample(1:N1,replace=FALSE)
t0 = t0[ind0]
t1 = t1[ind1]

cv.fold = rep(0,length(Y))
cv.fold[data$Y==0] = t0
cv.fold[data$Y==1] = t1

for(i in 1:length(lambdas)){
  lambda = lambdas[i]
  for(k in 1:K){
    train.ind <- (cv.fold!=k)
    val.ind <- (cv.fold==k)
    fit.data.SLL <- SuperLearner(Y=Y[train.ind], X=W[train.ind,],newX=W, SL.library = SL.library, family = binomial(), method = "method.NNLS",verbose = FALSE)
    sl.pred <- fit.data.SLL$SL.predict #prediction from super learner 
    lib.pred <- fit.data.SLL$library.predict #prediction from library algorithms
    pred <- cbind(sl.pred,lib.pred) #all predictions
    colnames(pred) <- c("SuperLearner","Random Forest","GLM","GAM","If-then Trees") # needs to change when SL.library varies
    val.S <- pred[val.ind,]
    val.Z <- Z[val.ind]

    assign(paste("CLSL",k,"_", i, sep=""),fit.data.SLL)
    assign(paste("CLval.S",k,"_", i, sep=""),val.S)
    assign(paste("CLval.Z", k,"_", i, sep=""),val.Z)
    print(c(lambda, k))
  }  
}

save.image("04-18Classification_temp.RData")


CVZ1 = CLval.Z1_1
for(k in 2:K){
  CVZ1 = c(CVZ1,get(paste("CLval.Z", k,"_1", sep="")))
}

all(CVZ == CVZ1) #same?




CLFPR = matrix(NA,ncol = 5, nrow = length(lambdas))
CLFNR = matrix(NA,ncol = 5, nrow = length(lambdas))
CLTPR = matrix(NA,ncol = 5, nrow = length(lambdas))
CLrisk = matrix(NA,ncol = 5, nrow = length(lambdas))

for(i in 1:length(lambdas)){
  lambda = lambdas[i]
  for(k in 1:K){
    val.Ski = get(paste("CLval.S",k,"_", i, sep=""))
    if(k == 1){
      CLCVSi = val.Ski
    } else {
      CLCVSi = rbind(CLCVSi, val.Ski)
    }
  }# for k fold
  CLCVYhati = (CLCVSi > 0.5)*1
  CLFPR[i,] = apply(CLCVYhati ,2,function(x) mean((x==1)*(1-CVZ))/mean(1-CVZ))
  CLFNR[i,] = apply(CLCVYhati ,2,function(x) mean((x==0)*(CVZ))/mean(CVZ))
  CLTPR[i,] = 1-CLFNR[i,]
  CLrisk[i,] = lambda*mean(CVZ)*CLFNR[i,] + (1-lambda)*(mean(1-CVZ))*CLFPR[i,]
}#for i lambda

CLrisk[,1] # classification risk of Super Learner
CombineRisk = cbind(CLrisk[,1], CombineRisk)
colnames(CombineRisk)[1] = "SLClassifier"
colnames(CombineRisk)[10] = "EnsembleCSuperLearner"
look = c(1,5,10,2)
dat = data.frame(cbind(lambdas,c(CombineRisk[,look])))
colnames(dat) <- c("lambda","risk")
dat$rules = rep(colnames(CombineRisk)[look],rep(length(lambdas),length(look)))
p <- ggplot(dat, aes(x=lambda, y=risk)) + geom_line(data=dat,aes(x=lambda, y=risk,group=rules,col = rules),lwd=1)
print(p)

temp = CombineRisk[,1]


look = c(1,5:9,11:14,10,2)
ind = c(1,5,9,13,17,21,25,29,39,49,59)
t(round(cbind(lambdas[ind],CombineRisk[ind,look]),3))


