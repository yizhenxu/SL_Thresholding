# 11-29-2016
# clean code for larger base learner library
# same code as CleanCode1, basically changing workspace names and specify new library

# 8 base learners:
# random forest, generalized linear regression, quadratic splines regression, CART, 
# 10 nearest neighbors, generalized boosting, support vector machine, Bagging Classification

# note that for Kenyan data, the svm need to redefine tuning parameter nu,
# this is done in the data preparation section, new wrapper with small
# nu is called SL.svm_nu

# Controlled Random Search use initiation from SL NNLS and constrain to [0, 5]^K

# CRS setting: maxeval = 10000, pop.size = 10000*(length(x0)+1)

# 1. prepare data (and save Kenyan data to a workspace)
# 2. clean code
###############################################################
# Some IMPORTANT notes for applications:

# Sample size is very important in this method application.
# Even the theoretical result is derived under asymptotic assumptions.
# For PIMA Indians Diabetes data, the sample size is less than 400,
# different seeds on the cross validation process has quite some influence on
# the cross-validated risk curves, in addition, it is important in this case to 
# re-tune the CRS parameters (e.g greatly increase maxeval) in order to 
# have satisfying performance.

########################################################################################
# 1. Prepare Data
###############################################################
###############################################################
# Kenyan data

setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\data")
load("analysis_data_10_7_2015.Rdata")
load("Long Crisis Study CD4.Rdata")
load("Long Second Line Study CD4.Rdata")
load("Long TDF Study CD4.Rdata")

colnames(crisis.cd4)[2] = "Study.ID"
second.cd4[,1] = as.factor(second.cd4[,1])
prodata = rbind(second.cd4[,c(1,3)], crisis.cd4[,c(2,5)], tdf.cd4[,c(1,3)])

# define min CD4 and max CD4 for each patient in every study
library(plyr)
minCD4 = ddply(prodata,"Study.ID",function(x) min(x$CD4Count))
maxCD4 = ddply(prodata,"Study.ID",function(x) max(x$CD4Count))

# define the full data with all complete cases
FullData = kenyadata 
slope.count <- (FullData[,5]-FullData[,17])/FullData[,16]
slope.perc <- (FullData[,6]-FullData[,19])/FullData[,18]
sndline = (FullData[,2]=="line2")
tdf = (FullData[,2]=="TDF")
pftadh = (FullData[,10]=="None")
female = (FullData[,11]=="F")
mydata <- cbind(FullData[,c(1,7,2,5,6,12,14,22)], female,pftadh,tdf,sndline,slope.count, slope.perc)

mydata <- mydata[complete.cases(mydata),]
dim(mydata) # 899 14

# merge minCD4 and maxCD4 into mydata
colnames(minCD4)[2] = "minCD4"
colnames(maxCD4)[2] = "maxCD4"
mydata = merge(minCD4, mydata, by = "Study.ID")
mydata = merge(maxCD4, mydata, by = "Study.ID")
dim(mydata)

# select covariates for model fitting 
Z = (mydata$ViralLoad_E>1000)*1 #define True Status
Y = Z #response used in finding risk score
W = as.data.frame(mydata[,c(8,11,3,6,7,12,9,16)]*1)
W[,c(1,3,4,5,7,8)]=apply(W[,c(1,3,4,5,7,8)],2,function(x) (x-mean(x))/sd(x))
W$female[W$female==0] = -1
W$pftadh[W$pftadh==0] = -1

setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\CRS constrained")
save(Y, W, file = "KenyanData.RData")

SL.svm_nu <- function(Y, X, newX, family, type.reg = "nu-regression", type.class = "nu-classification", kernel =
                     "radial", nu = 0.01, degree = 3, cost = 1, coef0 = 0, ...) {
  require('e1071')
  if(family$family == "gaussian") {
    fit.svm <- e1071::svm(y = Y, x = X, nu = nu, type = type.reg, fitted = FALSE, kernel = kernel, degree = degree, cost = cost, coef0 = coef0)
    pred <- predict(fit.svm, newdata = newX)
    fit <- list(object = fit.svm)
  }
  if(family$family == "binomial") {
    fit.svm <- e1071::svm(y = as.factor(Y), x = X, nu = nu, type = type.class, fitted = FALSE, probability = TRUE, kernel = kernel, degree = degree, cost = cost, coef0 = coef0)
    pred <- attr(predict(fit.svm, newdata = newX, probability = TRUE), "prob")[, "1"] # assumes Y is 0/1 numeric
    fit <- list(object = fit.svm)
  }
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.svm")
  return(out)
}

predict.SL.svm_nu <- function(object, newdata, family,...){
  require('e1071')
  if(family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata)
  }
  if(family$family == "binomial") {
    pred <- attr(predict(object$object, newdata = newdata, probability = TRUE), "prob")[, "1"]
  }
  return(pred)
}
###############################################################
###############################################################
# Wisconsin Breast Cancer

#setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\data")

#data = read.table("wdbc.data",sep = ",")
#colnames(data) = c("ID","diag",
#                   "MRadius","MTexture","MPerimeter","MArea","MSmooth",
#                   "MCompact","MConcavity","MConcaveP","MSymmetry","MFracDim",
#                   "SERadius","SETexture","SEPerimeter","SEArea","SESmooth",
#                   "SECompact","SEConcavity","SEConcaveP","SESymmetry","SEFracDim",
#                   "WRadius","WTexture","WPerimeter","WArea","WSmooth",
#                   "WCompact","WConcavity","WConcaveP","WSymmetry","WFracDim")

#W = matrix(unlist(data[,3:32]),ncol=dim(data[,3:32])[2])
#colnames(W) = colnames(data)[3:32]
#W = apply(W,2,function(x) (x-mean(x))/sd(x))
#W = as.data.frame(W)

#Y = rep(0,dim(data)[1])
#Y[data[,2]=="B"] = 1
#Z = Y
#table(Y)

###############################################################
###############################################################
# PIMA


setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\data")
setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\data")

data = read.table("pima-indians-diabetes.data",sep = ",")
dim(data) #768   9
head(data)
summary(data)

Y = data[,ncol(data)]
W = data[,-ncol(data)]
summary(W)

# we need to make all variables mean 0 sd 1
# for continuous variables, simply standardize 
# for categorical varialbes, (0,1)->(-1,1), (0,1,2)->(-1,0,1)

# continuous columns: all W


W  = apply(W,2,function(x) (x-mean(x))/sd(x))

W = as.data.frame(W)


table(Y)
#0   1 
#500 268 
mean(Y)
#  prevelence 0.3489583


########################################################################################

# 2.clean code
# load in functions_CRSnes.R -- changes made to crs.fit: 
#   (1)initiation using NNLS
#   (2)region [0, 5]^K
#   (3)normalize to sum up to 1 before output and making predictions
# simply change the load data step for other runs
# data format: Y -- outcome; W -- covariates in dataframe
# other parameters for tuning:
#     number of folds K, SL library, number of minimization alg

setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\CRS constrained")
setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\CRS constrained")

#load("KenyanData.RData")# load data

# load packages
library("SuperLearner")
library("randomForest")
library("gam")
library("rpart")
library("dplyr")
library("plyr")
library("ggplot2")
library("nloptr")
library("lpSolve")
library("nnls")

source("functions_CRSnew.r")# load functions
# svm_nu for Kenyan Daya
SL.library = c("SL.randomForest","SL.glm","SL.gam", "SL.rpart","SL.knn","SL.gbm","SL.svm_nu","SL.ipredbagg")
# regular svm for PIMA
SL.library = c("SL.randomForest","SL.glm","SL.gam", "SL.rpart","SL.knn","SL.gbm","SL.svm","SL.ipredbagg")
# random forest, generalized linear regression, quadratic splines regression, CART, 
# 10 nearest neighbors, generalized boosting, support vector machine, Bagging Classification


###### cross validation
###### balance the ratio of Y=0 and Y=1's in each fold

N = length(Y)
K = 10 # number of folds for calculating cross validated risk

N0 = table(Y)[1] # number of {Y=0}
l0 = floor(N0/K) # number of {Y=0} obs in each fold
# evenly distribute the leftovers, fold label vector:
t0 = rep(1:K,rep(l0,K)+c(rep(1,N0%%l0),rep(0,K-N0%%l0))) 

N1 = table(Y)[2] # number of {Y=1}
l1 = floor(N1/K) # number of {Y=0} obs in each fold
t1 = rep(1:K,rep(l1,K)+c(rep(1,N1%%l1),rep(0,K-N1%%l1)))

data = cbind(Y,W)

#as.integer(runif(1)*2e9)
#set.seed(1818715699)#Kenyan
#set.seed(100) #WDBC
set.seed(852632362)#PIMA
# also ran seed = 49836901 (riskPIMA_bigLib1.RData) and seed = 1566626083 (riskPIMA_bigLib2.RData)
ind0 = sample(1:N0,replace=FALSE) 
ind1 = sample(1:N1,replace=FALSE)
# permute the fold label vector
t0 = t0[ind0] 
t1 = t1[ind1]

# permute the Y=0 and Y=1's separately to make balance in each fold
# cv.fold is the validation fold index for the whole dataset
cv.fold = rep(0,length(Y))
cv.fold[Y==0] = t0
cv.fold[Y==1] = t1


for(i in 1:K){
  print(i)
  train.ind <- (cv.fold!=i)
  val.ind <- (cv.fold==i)
  fit.data.SLL <- SuperLearner(Y=Y[train.ind], X=W[train.ind,],newX=W, SL.library = SL.library, family = binomial(), method = "method.NNLS",verbose = FALSE)
  sl.pred <- fit.data.SLL$SL.predict #prediction from super learner 
  lib.pred <- fit.data.SLL$library.predict #prediction from library algorithms
  pred <- cbind(sl.pred,lib.pred) #all predictions
  colnames(pred) <- c("SuperLearner",SL.library) 
  train.S <- pred[train.ind,] # trained predictions
  val.S <- pred[val.ind,] # validation predictions
  train.CVS <- fit.data.SLL$Z # cross-validated library predictions 
  trainCVS_SL = fit.data.SLL$Z %*% fit.data.SLL$coef # SL CV predictions from within SL
  train.CVS.wSL = cbind(trainCVS_SL, train.CVS) # add the first col to be the SL CV predictions
  
  train.Z <- Y[train.ind] # trained outcome
  val.Z <- Y[val.ind] # validation outcome
  
  assign(paste("train.CVS.wSL", i, sep=""),train.CVS.wSL)
  assign(paste("train.CVS", i, sep=""),train.CVS)
  assign(paste("train.S", i, sep=""),train.S)
  assign(paste("val.S", i, sep=""),val.S)
  assign(paste("train.Z", i, sep=""),train.Z)
  assign(paste("val.Z", i, sep=""),val.Z)
}  

###### Solving for (alpha, c) using alg number of methods
alg=3
#as.integer(runif(1)*2e9)
# randomly generated number, for replicating results
#set.seed(1988420473)#Kenyan data
#set.seed(1920413227)#WDBC
set.seed(17187750)#PIMA

lambdas = unique(seq(0.1,0.9,0.01))
#lambdas = unique(c(seq(0,0.7,0.025),seq(0.7,1,0.01)))
# alg algorithms, length(lambdas) lambdas, 10 folds
cutoff = array(NA, dim = c(alg,length(lambdas),K))

# true status vector that follows the order of the stacked cross validated predictions
CVZ = val.Z1
for(k in 2:K){
  CVZ = c(CVZ,get(paste("val.Z", k, sep="")))
}

FPR = matrix(NA,ncol = alg, nrow = length(lambdas))
FNR = matrix(NA,ncol = alg, nrow = length(lambdas))
TPR = matrix(NA,ncol = alg, nrow = length(lambdas))
risk = matrix(NA,ncol = alg, nrow = length(lambdas))

deci =  vector("list", length(lambdas)) 

for(k in 1:K){
  
  train.CVS = get(paste("train.CVS", k, sep=""))
  train.S = get(paste("train.S", k, sep=""))
  val.S = get(paste("val.S", k, sep=""))
  train.Z = get(paste("train.Z", k, sep=""))
  val.Z = get(paste("val.Z", k, sep=""))
  train.CVS.wSL = get(paste("train.CVS.wSL", k, sep=""))
  XX = train.S[,-1] # training fold library predictions
  
  for(i in 1:length(lambdas)){
    print(c(i,k))
    lambda = lambdas[i]
    crs.seed = as.integer(runif(1)*2e9)
    cbs = crs.fit(seed = crs.seed, lambda,train.CVS,train.Z, val.S[,-1], val.Z)
    cutoff[3, i, k] = cbs$c #CBS cutoff for lambda i fold k
    
    # SL Common
    opt = Opt.nonpar.rule(train.Z,train.S[,1],phi=0,lambda)
    cutoff[1, i, k] = as.numeric(opt)[1]
    # SL Proposal Iterative
    opt = Opt.nonpar.rule(train.Z,train.CVS.wSL[,1],phi=0,lambda)
    cutoff[2, i, k] = as.numeric(opt)[1]
    
    cut = matrix(rep(cutoff[, i, k],nrow(val.S)),nrow=nrow(val.S),byrow=T)
    val = cbind(val.S[,1], val.S[,1], cbs$score)
    dki = (val > cut)*1  
    # one decision matrix for each lambda value, matrix size n x alg
    if(k == 1){
      deci[[i]] = dki
    }
    if(k>1){
      deci[[i]] = rbind(deci[[i]],dki) # each row is validation decision arranged by folds, col is algorithm
    }
  }#i
}#k


for(i in 1:length(lambdas)){
  lambda = lambdas[i]
  dec = deci[[i]]
  FPR[i,] = apply(dec,2,function(x) mean((x==1)*(1-CVZ))/mean(1-CVZ))
  FNR[i,] = apply(dec,2,function(x) mean((x==0)*(CVZ))/mean(CVZ))
  risk[i,] = lambda*mean(CVZ)*FNR[i,] + (1-lambda)*(mean(1-CVZ))*FPR[i,]
  
}

toplot = risk
colnames(toplot) = c("Conditional Thresholding", "Two-Step Minimization", "CRS Minimization")
dat1 = data.frame(cbind(lambdas,c(toplot)))
colnames(dat1) <- c("lambda","risk")
dat1$Approach = rep(colnames(toplot),rep(length(lambdas),ncol(toplot)))
p1 <- ggplot(dat1, aes(x=lambda, y=risk)) + geom_line(data=dat1,aes(x=lambda, y=risk,group=Approach,col = Approach),lwd=1)

print(p1)


ind = which(lambdas %in% seq(0.1,0.9,0.1))
round(t(cbind(lambdas,risk)[ind,]),3)


#############################################################################
# save workspace


riskPIMA = toplot
lambdasPIMA = lambdas
save(lambdasPIMA, riskPIMA, file = "riskPIMA_bigLib.RData")
save.image("CRSnew_PIMAData_bigLib.RData")


#setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\CRS constrained")
#load("CRSnew_KenyanData.RData")
riskKenya = toplot
lambdaKenya = lambdas
save(lambdaKenya, riskKenya, file = "riskKenya_bigLib.RData")
save.image("CRSnew_KenyanData_bigLib.RData")


## the WDBC is already done in wdbc_largerLibrary.R (11-28-2016)
# riskWDBC = toplot
# lambdasWDBC = lambdas
# save(lambdasWDBC, riskWDBC, file = "riskWDBC1.RData")
# save.image("CRSnew_WDBCData_bigLib.RData")

###############################################################

###############################################################

###############################################################

###############################################################

###############################################################

###############################################################

###############################################################

rm(list = ls())
###############################################################



setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\CRS constrained")
setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\CRS constrained")
# load packages
library("SuperLearner")
library("randomForest")
library("gam")
library("rpart")
library("dplyr")
library("plyr")
library("ggplot2")
library("nloptr")
library("lpSolve")
library("nnls")

load("riskKenya_bigLib.RData")
load("riskPIMA_bigLib.RData")
load("riskWDBC_bigLib.RData")


toplot = riskKenya
colnames(toplot) = c("Conditional Thresholding", "Two-Step Minimization", "CRS Minimization")
lambdas = lambdaKenya
dat3 = data.frame(cbind(lambdas,c(toplot)))
colnames(dat3) <- c("lambda","risk")
dat3$Approach = rep(colnames(toplot),rep(length(lambdas),ncol(toplot)))
p3 <- ggplot(dat3, aes(x=lambda, y=risk)) + geom_line(data=dat3,aes(x=lambda, y=risk,group=Approach,col = Approach),lwd=1) + ggtitle("Kenyan Study")
p3 <- p3 + xlab(expression(paste(lambda)))
print(p3)

colnames(riskWDBC) = c("Conditional Thresholding", "Two-Step Minimization", "CRS Minimization")
toplot = riskWDBC
lambdas = lambdasWDBC
dat4 = data.frame(cbind(lambdas,c(toplot)))
colnames(dat4) <- c("lambda","risk")
dat4$Approach = rep(colnames(toplot),rep(length(lambdas),ncol(toplot)))
p4 <- ggplot(dat4, aes(x=lambda, y=risk)) + geom_line(data=dat4,aes(x=lambda, y=risk,group=Approach,col = Approach),lwd=1) + ggtitle("Breast Cancer Study")
p4 <- p4 + xlab(expression(paste(lambda)))
print(p4)

toplot = riskPIMA
colnames(toplot) = c("Conditional Thresholding", "Two-Step Minimization", "CRS Minimization")
lambdas = lambdasPIMA
dat5 = data.frame(cbind(lambdas,c(toplot)))
colnames(dat5) <- c("lambda","risk")
dat5$Approach = rep(colnames(toplot),rep(length(lambdas),ncol(toplot)))
p5 <- ggplot(dat5, aes(x=lambda, y=risk)) + geom_line(data=dat5,aes(x=lambda, y=risk,group=Approach,col = Approach),lwd=1) + ggtitle("PIMA Diabetes Study")
p5 <- p5 + xlab(expression(paste(lambda)))
print(p5)

library(ggplot2)
library(gridExtra)
library(grid)
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


png(filename = "RiskPlot_Application_bigLib.png", width = 550, height = 850)

grid_arrange_shared_legend(p3, p4,p5)

dev.off()

which(lambdaKenya %in% seq(0.1,0.9, 0.1))
ind = c(11, 21, 31,41, 51, 61, 71, 81, 91)
round(t(cbind(lambdaKenya,riskKenya)[ind,]),3)

ind = c(5,9,13,17,21,25,29,33,37)
round(t(cbind(lambdasWDBC,riskWDBC)[ind,]),3)

ind = c(11, 21, 31,41, 51, 61, 71, 81, 91)
round(t(cbind(lambdasPIMA,riskPIMA)[ind,]),3)

ind = which(lambdaKenya %in% c(0.2,0.5,0.8))
round(t(cbind(lambdaKenya,riskKenya)[ind,]),3)
ind = which(lambdasWDBC %in% c(0.2,0.5,0.8))
round(t(cbind(lambdasWDBC,riskWDBC)[ind,]),3)
