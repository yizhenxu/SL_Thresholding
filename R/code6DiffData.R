
library("SuperLearner")
library("randomForest")
library("gam")
library("rpart")
library("dplyr")
library("plyr")
library("ggplot2")
library("nloptr")
library("lpSolve")

# Data Preparation
setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\data")
setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up\\data")

###############################################################
###############################################################
# MUSHROOM
data = read.table("agaricus-lepiota.data",sep = ",")
library(ade4)
mushroom = acm.disjonctif(data)[,-c(1,3,9,13,23,25,34,36,38,40,52,54,59,63,67,76,86,90,93,98,107,113)]
for(i in 2:ncol(mushroom)){
  mushroom[mushroom[,i]==0,i] = -1
}
head(mushroom)
Z = mushroom[,1]
Y = Z #response used in finding risk score
W = mushroom[,2:ncol(mushroom)]
###############################################################
###############################################################

# Wisconsin Breast Cancer
data = read.table("wdbc.data",sep = ",")
colnames(data) = c("ID","diag",
                       "MRadius","MTexture","MPerimeter","MArea","MSmooth",
                       "MCompact","MConcavity","MConcaveP","MSymmetry","MFracDim",
                       "SERadius","SETexture","SEPerimeter","SEArea","SESmooth",
                       "SECompact","SEConcavity","SEConcaveP","SESymmetry","SEFracDim",
                       "WRadius","WTexture","WPerimeter","WArea","WSmooth",
                       "WCompact","WConcavity","WConcaveP","WSymmetry","WFracDim")

W = matrix(unlist(data[,3:32]),ncol=dim(data[,3:32])[2])
colnames(W) = colnames(data)[3:32]
W = apply(W,2,function(x) (x-mean(x))/sd(x))
W = as.data.frame(W)

Y = rep(0,dim(data)[1])
Y[data[,2]=="B"] = 1
Z = Y
table(Y)


###############################################################
###############################################################

setwd("C:\\Users\\yxu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up")
setwd("C:\\Users\\Yizhen Xu\\Google Drive\\Desktop\\2015 summer\\Risk Score Wrap Up")
source("functions.r")
SL.library = c("SL.randomForest","SL.glm","SL.gam", "SL.rpart")

###############################################################
###############################################################
###############################################################
###############################################################
### cross validation

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
cv.fold[Y==0] = t0
cv.fold[Y==1] = t1

for(i in 1:K){
  
  train.ind <- (cv.fold!=i)
  val.ind <- (cv.fold==i)
  fit.data.SLL <- SuperLearner(Y=Y[train.ind], X=W[train.ind,],newX=W, SL.library = SL.library, family = binomial(), method = "method.NNLS",verbose = FALSE)
  sl.pred <- fit.data.SLL$SL.predict #prediction from super learner 
  lib.pred <- fit.data.SLL$library.predict #prediction from library algorithms
  pred <- cbind(sl.pred,lib.pred) #all predictions
  colnames(pred) <- c("SuperLearner","Random Forest","GLM","GAM","If-then Trees") # needs to change when SL.library varies
  train.S <- pred[train.ind,]
  val.S <- pred[val.ind,]
  train.CVS <- fit.data.SLL$Z
  trainCVS_SL = fit.data.SLL$Z %*% fit.data.SLL$coef
  train.CVS.wSL = cbind(trainCVS_SL, train.CVS) # add the first col to be the SL CV predictions
  
  train.Z <- Z[train.ind]
  val.Z <- Z[val.ind]
  
  assign(paste("train.CVS.wSL", i, sep=""),train.CVS.wSL)
  assign(paste("train.CVS", i, sep=""),train.CVS)
  assign(paste("train.S", i, sep=""),train.S)
  assign(paste("val.S", i, sep=""),val.S)
  assign(paste("train.Z", i, sep=""),train.Z)
  assign(paste("val.Z", i, sep=""),val.Z)
}  


###########################################################

as.integer(runif(1)*2e9)

#set.seed(530575739)
set.seed(1920413227)

# 8 algorithms, length(lambdas) lambdas, 10 folds
lambdas = unique(c(seq(0,0.7,0.025),seq(0.7,1,0.01)))
alg=8
cutoff = array(NA, dim = c(alg,length(lambdas),K))

# calculate the training and validation scores

for(k in 1:K){
  train.CVS = get(paste("train.CVS", k, sep=""))
  train.S = get(paste("train.S", k, sep=""))
  val.S = get(paste("val.S", k, sep=""))
  train.Z = get(paste("train.Z", k, sep=""))
  val.Z = get(paste("val.Z", k, sep=""))
  
  for(i in 1:length(lambdas)){
    lambda = lambdas[i]
    crs.seed = as.integer(runif(1)*2e9)
    cbs = crs.fit(seed = crs.seed, lambda,train.CVS,train.Z, val.S[,-1], val.Z)
    isr = isres.fit(lambda,train.CVS,train.Z, val.S[,-1], val.Z)
    srg = srgFun(lambda,train.CVS, train.Z, val.S[,-1], val.Z) # l1 svm
    XX = train.S[,-1]
    X.SL = as.matrix(XX[,cbs$ord])
    trainS_cbs = (X.SL%*%cbs$b) #training score
    trainS_isr = (X.SL%*%isr$b)
    trainS_srg = (as.matrix(XX)%*%srg$b)
    train.Si = cbind(trainS_cbs, trainS_isr, trainS_srg, train.S)
    val.Si = cbind(cbs$score, isr$score, srg$score, val.S)
    colnames(train.Si) = colnames(val.Si) = c("CRS", "ISRES", "L1 SVM", colnames(train.S))
    cutoff[1,i,k] = cbs$c #CBS cutoff for lambda i fold k
    cutoff[2,i,k] = isr$c #ISRES cutoff for lambda i fold k
    cutoff[3,i,k] = srg$c # L1 cutoff for lambda i fold k
    assign(paste("train.S", k,"_", i, sep=""),train.Si)
    assign(paste("val.S", k,"_", i, sep=""),val.Si)
  }#i
}#k

# calculate the cutoffs from training score for each lambda i, fold k, alrogithm j
# using Opt.nonpar.rule
# CBS and L1 cutoffs are given with the score calculations
# apply the trained cutoffs to the corresponding validation score
# obtain the validation DECISIONS, stack up among the 10 folds
# Then we have CV decisions for all the 7 algorithms at the length(lambdas) lambdas
# Calculate FPR, TPR, Risk for the 7xlength(lambdas) CV decisions (length N),
# add (FPR, TPR) as points to ROC plot
# Risk would be as tables

# true status vector that follows the order of the stacked cross validated predictions
CVZ = val.Z1
for(k in 2:K){
  CVZ = c(CVZ,get(paste("val.Z", k, sep="")))
}

FPR = matrix(NA,ncol = alg, nrow = length(lambdas))
FNR = matrix(NA,ncol = alg, nrow = length(lambdas))
TPR = matrix(NA,ncol = alg, nrow = length(lambdas))
risk = matrix(NA,ncol = alg, nrow = length(lambdas))
for(i in 1:length(lambdas)){
  lambda = lambdas[i]
  for(k in 1:K){
    train.Ski = get(paste("train.S", k,"_", i, sep=""))
    val.Ski = get(paste("val.S", k,"_", i, sep=""))
    train.Z = get(paste("train.Z", k, sep=""))
    val.Z = get(paste("val.Z", k, sep=""))
    for(j in 1:ncol(train.Ski)){
      if(j>3){
        opt = Opt.nonpar.rule(train.Z,train.Ski[,j],phi=0,lambda)
        cutoff[j,i,k] = as.numeric(opt)[1]
      }
    } # j algorithm
    cut = matrix(rep(cutoff[,i,k],nrow(val.Ski)),nrow=nrow(val.Ski),byrow=T)
    dki = (val.Ski > cut)*1  
    if(k == 1){
      deci = dki
    }
    if(k>1){
      deci = rbind(deci,dki) # each row is a fold, col is algorithm
    }
  }# for k fold
  FPR[i,] = apply(deci,2,function(x) mean((x==1)*(1-CVZ))/mean(1-CVZ))
  FNR[i,] = apply(deci,2,function(x) mean((x==0)*(CVZ))/mean(CVZ))
  TPR[i,] = 1-FNR[i,]
  risk[i,] = lambda*mean(CVZ)*FNR[i,] + (1-lambda)*(mean(1-CVZ))*FPR[i,]
  assign(paste("decision_", i, sep=""),deci)
}#for i lambda

colnames(risk) = colnames(FNR) = colnames(FPR) = colnames(TPR) = colnames(train.Ski)

CValg=ncol(train.CVS.wSL)
CVcutoff = array(NA, dim = c(CValg,length(lambdas),10))

CVCFPR = matrix(NA,ncol = CValg, nrow = length(lambdas))
CVCFNR = matrix(NA,ncol = CValg, nrow = length(lambdas))
CVCTPR = matrix(NA,ncol = CValg, nrow = length(lambdas))
CVCrisk = matrix(NA,ncol = CValg, nrow = length(lambdas))

for(i in 1:length(lambdas)){
  lambda = lambdas[i]
  for(k in 1:K){
    train.CVS.wSLk = get(paste("train.CVS.wSL", k, sep=""))
    #train.Ski = get(paste("train.S", k,"_", i, sep=""))
    val.Ski = get(paste("val.S", k,"_", i, sep=""))
    train.Z = get(paste("train.Z", k, sep=""))
    val.Z = get(paste("val.Z", k, sep=""))
    for(j in 1:ncol(train.CVS.wSLk)){
      opt = Opt.nonpar.rule(train.Z,train.CVS.wSLk[,j],phi=0,lambda)
      CVcutoff[j,i,k] = as.numeric(opt)[1]
    } # j algorithm in CVpred with SL
    CVCcut = matrix(rep(CVcutoff[,i,k],nrow(val.Ski)),nrow=nrow(val.Ski),byrow=T)
    CVCdki = (val.Ski[,4:8] > CVCcut)*1  
    if(k == 1){
      CVCdeci = CVCdki
    }
    if(k>1){
      CVCdeci = rbind(CVCdeci, CVCdki) # each row is a fold, col is algorithm
    }
  }# for k fold
  CVCFPR[i,] = apply(CVCdeci,2,function(x) mean((x==1)*(1-CVZ))/mean(1-CVZ))
  CVCFNR[i,] = apply(CVCdeci,2,function(x) mean((x==0)*(CVZ))/mean(CVZ))
  CVCTPR[i,] = 1-CVCFNR[i,]
  CVCrisk[i,] = lambda*mean(CVZ)*CVCFNR[i,] + (1-lambda)*(mean(1-CVZ))*CVCFPR[i,]
  assign(paste("CVCdecision_", i, sep=""),CVCdeci)
}#for i lambda

colnames(CVCrisk) = colnames(CVCFNR) = colnames(CVCFPR) = colnames(CVCTPR) = colnames(pred)

round(CVCrisk,3)
round(risk,3)
#save.image("mushroom_all.RData")
save.image("wdbc_all.RData")


CombineRisk = cbind(risk,CVCrisk)
colnames(CombineRisk) = c("CRS", "ISRES", "L1 SVM","SuperLearner","Random Forest","GLM","GAM","If-then Trees","CVCSuperLearner","CVCRandom Forest","CVCGLM","CVCGAM","CVCIf-then Trees")

look = c(4,9)
look = c(5,10)
look = c(6,11)
look = c(7,12)
look = c(8,13)
look = c(1,2,3,4,5,6,7,8,9)
look = c(1:3,9)

look = c(1,4,5,9,10)
look = c(1,4,9)

look = c(4,9,1)
temp = CombineRisk[,look]
colnames(temp) = c("SL Common Approach", "SL Proposal Iterative", "SL Proposal Controlled Search")
dat = data.frame(cbind(lambdas,c(temp)))
colnames(dat) <- c("lambda","risk")
dat$rules = rep(colnames(temp),rep(length(lambdas),length(look)))
p <- ggplot(dat, aes(x=lambda, y=risk)) + geom_line(data=dat,aes(x=lambda, y=risk,group=rules,col = rules),lwd=1)
print(p)

ind = c(5,9,13,17,21,25,29,39,49)
round(t(temp[ind,]),3)

