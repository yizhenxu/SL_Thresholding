# 8 base learners:
# random forest, generalized linear regression, quadratic splines regression, CART, 
# 10 nearest neighbors, generalized boosting, support vector machine, Bagging Classification

# load in functions_CRSnes.R -- changes made to crs.fit: 
#   (1)initiation using NNLS
#   (2)region [0, 5]^K
#   (3)normalize to sum up to 1 before output and making predictions
# simply change the load data step for other runs
# data format: Y -- outcome; W -- covariates in dataframe
# other parameters for tuning:
#     number of folds K, SL library, number of minimization alg
################################################################

### Data Processing

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
#table(Y)

################################################################

### Analysis

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
SL.library = c("SL.randomForest","SL.glm","SL.gam", "SL.rpart","SL.knn","SL.gbm","SL.svm","SL.ipredbagg")

# cross validation
# balance the ratio of Y=0 and Y=1's in each fold
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

set.seed(100) #for reproducibility, set seed = 100

# permute the fold label vector
ind0 = sample(1:N0,replace=FALSE) 
ind1 = sample(1:N1,replace=FALSE)

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

# Solving for (alpha, c) using alg number of methods
alg=3 # CRS, Two steps, conditional thresholding
#as.integer(runif(1)*2e9)
set.seed(1920413227)#WDBC


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
    
    # SL CRS
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





