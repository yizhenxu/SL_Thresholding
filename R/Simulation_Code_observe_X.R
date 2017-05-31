###############################################################
# 05-30-2017
# clean code for simulations case 2

# 4 base learners:
# random forest, generalized linear regression, quadratic splines regression, CART

# source in functions.R 
# CRS: crs.fit function in functions.R for parameter changes: 
#   (1)initiation using NNLS
#   (2)region [0, 5]^K
#   (3)normalize to sum up to 1 before output and making predictions

###############################################################

rm(list = ls())

### Case 2

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

source("functions.r")# load functions
SL.library = c("SL.randomForest","SL.glm","SL.gam", "SL.rpart")

set.seed(1000)

NN = (10^4)
nsize = 2*NN # generate  D1 and D2, each of size NN
kz = 4
set.seed(1)
z = rnorm(kz*nsize, mean=0, sd=1)
z = matrix(z, ncol=kz)

dimnames(z)[[2]] <- c("Z1", "Z2", "Z3", "Z4")

W = matrix( rep(0,nsize*4), ncol=4 )

W[,1] = exp(z[,1] / 2)
W[,2] = z[,2] / (1 + exp(z[,1])) + 10
W[,3] = ( z[,1]*z[,3]/25 + .6)^3
W[,4] = ( z[,2] + z[,4] + 20 )^2
W[,1:4]   = apply(W[,1:4],2,function(x) (x-mean(x))/sd(x))
dimnames(W)[[2]] <- c("X1", "X2", "X3", "X4")

W = as.data.frame(W)

beta      = matrix( c(210, 27.4, 13.7, 13.7, 13.7), ncol=1 )
epsilon   = matrix( rnorm(nsize, mean=0, sd=100), ncol=1 )
YY        = cbind(1,z) %*% beta + epsilon # underlying score

hist(YY)


temp1 = unique(sort(YY))
cut1 = floor(temp1[floor(nsize*0.7)])
cut1
b = c(27.4, 13.7, 13.7, 13.7)
sig = sqrt(sum(b^2)+100^2)
#cut1 = qnorm(0.7, 210, sig)
#cut1
Y = (YY >= cut1)*1 # outcome
table(Y)
hist(YY)
abline(v = cut1)



# Fit the algorithm on the first half of the data
# Fit on D1 = (W1, Y1), get predictions on the whole simulated data D = D1 U D2 = (W, Y)
fit.data.SLL <- SuperLearner(Y=Y[1:NN], X=W[1:NN,], newX=W, SL.library = SL.library, family = binomial(), method = "method.NNLS",verbose = FALSE)
sl.pred <- fit.data.SLL$SL.predict #prediction from super learner 
lib.pred <- fit.data.SLL$library.predict #prediction from library algorithms
pred <- cbind(sl.pred,lib.pred) #all predictions
colnames(pred) <- c("SuperLearner","Random Forest","GLM","GAM","If-then Trees") # needs to change when SL.library varies
train.S <- pred[1:NN,] # predictions of SL and candidate learners on D1
val.S <- pred[(NN+1):(2*NN),] # predictions of SL and candidate learners on D2
train.CVS <- fit.data.SLL$Z # Z of D1
trainCVS_SL = fit.data.SLL$Z %*% fit.data.SLL$coef # Z\alpha of D1
train.CVS.wSL = cbind(trainCVS_SL, train.CVS) # (Z\alpha, Z) of D1

train.Z <- Y[1:NN] # outcome in D1
val.Z <- Y[(NN+1):(2*NN)] # outcome in D2


set.seed(1920413227)
lambdas = unique(seq(0,1,0.025))
alg=3 # crs, SLProposal, SLCommon
cutoff = matrix(NA,nrow = alg, ncol = length(lambdas))

# Apply to the second half of the data to get the true risk

CVZ = val.Z

FPR = matrix(NA,ncol = alg, nrow = length(lambdas))
FNR = matrix(NA,ncol = alg, nrow = length(lambdas))
TPR = matrix(NA,ncol = alg, nrow = length(lambdas))
risk = matrix(NA,ncol = alg, nrow = length(lambdas))

XX = train.S[,-1]
for(i in 1:length(lambdas)){
  print(i)
  lambda = lambdas[i]
  
  # SL CRS
  crs.seed = as.integer(runif(1)*2e9)
  cbs = crs.fit(seed = crs.seed, lambda,train.CVS,train.Z, val.S[,-1], val.Z) # alpha and c estimated based on Z and Y of D1 (train.CVS and train.z)
  cutoff[3,i] = cbs$c #CBS cutoff for lambda i 
  
  # SL Common
  opt = Opt.nonpar.rule(train.Z,train.S[,1],phi=0,lambda) # c estimated by line search based on Y and \hat{Y}_{SL} of D1
  cutoff[1,i] = as.numeric(opt)[1]
  
  # SL Proposal Iterative
  opt = Opt.nonpar.rule(train.Z,train.CVS.wSL[,1],phi=0,lambda) # c estimated by line search based on Y and Z\hat{\alpha} of D1 
  cutoff[2,i] = as.numeric(opt)[1]
  
  cut = matrix(rep(cutoff[,i],nrow(val.S)),nrow=nrow(val.S),byrow=T)
  val = cbind(val.S[,1], val.S[,1], cbs$score)
  deci = (val > cut)*1  
  
  FPR[i,] = apply(deci,2,function(x) mean((x==1)*(1-CVZ))/mean(1-CVZ))
  FNR[i,] = apply(deci,2,function(x) mean((x==0)*(CVZ))/mean(CVZ))
  TPR[i,] = 1-FNR[i,]
  risk[i,] = lambda*mean(CVZ)*FNR[i,] + (1-lambda)*(mean(1-CVZ))*FPR[i,]
  assign(paste("decision_", i, sep=""),deci)
}#for i lambda


# calculate true probability score \Psi_0(U) 
epsig = 100
MY = cbind(1,z) %*% beta
S = pnorm(cut1-MY, mean = 0, sd = epsig,lower.tail=FALSE)
S2 = S[(NN+1):(2*NN)]

# true risk calculation (risk for the optimal classification rule)
TR2 = c()
for(i in 1:length(lambdas)){
  TR2[i] = mean(lambdas[i]*Y[(NN+1):(2*NN)]*(S2<(1-lambdas[i]))+(1-lambdas[i])*(1-Y[(NN+1):(2*NN)])*(S2>=(1-lambdas[i])))
}

draw = cbind(risk,TR2)
colnames(draw) = c("Conditional Thresholding", "Two-Step Minimization", "CRS Minimization","True Classification Rule")
round(draw,4)
toplot = draw
dat = data.frame(cbind(lambdas,c(toplot)))
colnames(dat) <- c("lambda","risk")
dat$rules = rep(colnames(toplot),rep(length(lambdas),ncol(toplot)))
p <- ggplot(dat, aes(x=lambda, y=risk)) + geom_line(data=dat,aes(x=lambda, y=risk,group=rules,col = rules),lwd=1)
print(p)

ind = which(lambdas %in% seq(0.1,0.9,0.1))
round(t(cbind(lambdas,draw)[ind,]),3)


