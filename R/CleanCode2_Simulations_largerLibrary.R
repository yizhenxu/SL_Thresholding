###############################################################
# 11-30-2016
# clean code for simulations
# same code as CleanCode1, basically changing workspace names and specify new library

# 8 base learners:
# random forest, generalized linear regression, quadratic splines regression, CART, 
# 10 nearest neighbors, generalized boosting, support vector machine, Bagging Classification

# Controlled Random Search use initiation from SL NNLS and constrain to [0, 5]^K

# CRS setting: maxeval = 10000, pop.size = 10000*(length(x0)+1)

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

source("functions_CRSnew.r")# load functions
SL.library = c("SL.randomForest","SL.glm","SL.gam", "SL.rpart","SL.knn","SL.gbm","SL.svm","SL.ipredbagg")


###############################################################
###############################################################
set.seed(1000)
# results r1
NN = (10^4)
nsize = 2*NN
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
YY        = cbind(1,z) %*% beta + epsilon

hist(YY)


###############################################################

temp1 = unique(sort(YY))
cut1 = floor(temp1[floor(nsize*0.7)])
cut1
b = c(27.4, 13.7, 13.7, 13.7)
sig = sqrt(sum(b^2)+100^2)
cut1 = qnorm(0.7, 210, sig)
cut1
Y = (YY >= cut1)*1
table(Y)
hist(YY)
abline(v = cut1)



###############################################################

# Fit the algorithm on the first half of the data
# Fit on (W1, Y1), get predictions on the whole simulated data (W, Y)
fit.data.SLL <- SuperLearner(Y=Y[1:NN], X=W[1:NN,], newX=W, SL.library = SL.library, family = binomial(), method = "method.NNLS",verbose = FALSE)
sl.pred <- fit.data.SLL$SL.predict #prediction from super learner 
lib.pred <- fit.data.SLL$library.predict #prediction from library algorithms
pred <- cbind(sl.pred,lib.pred) #all predictions
colnames(pred) <- c("SuperLearner","Random Forest","GLM","GAM","If-then Trees") # needs to change when SL.library varies
train.S <- pred[1:NN,]
val.S <- pred[(NN+1):(2*NN),]
train.CVS <- fit.data.SLL$Z
trainCVS_SL = fit.data.SLL$Z %*% fit.data.SLL$coef
train.CVS.wSL = cbind(trainCVS_SL, train.CVS) # add the first col to be the SL CV predictions

train.Z <- Y[1:NN]
val.Z <- Y[(NN+1):(2*NN)]


set.seed(1920413227)
# 3 algorithms, length(lambdas) lambdas
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
  crs.seed = as.integer(runif(1)*2e9)
  cbs = crs.fit(seed = crs.seed, lambda,train.CVS,train.Z, val.S[,-1], val.Z)
  cutoff[3,i] = cbs$c #CBS cutoff for lambda i 
  
  # SL Common
  opt = Opt.nonpar.rule(train.Z,train.S[,1],phi=0,lambda)
  cutoff[1,i] = as.numeric(opt)[1]
  # SL Proposal Iterative
  opt = Opt.nonpar.rule(train.Z,train.CVS.wSL[,1],phi=0,lambda)
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


epsig = 100

MY = cbind(1,z) %*% beta
#S = pnorm(cut1, mean = MY, sd = epsig,lower.tail=FALSE)
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

ind = c(5,9,13,17,21,25,29,33,37)
round(t(cbind(lambdas,draw)[ind,]),3)

trisk_X = draw
tlambda_X = lambdas
save(tlambda_X, trisk_X, file = "TCR_X1_bigLib.RData")
save.image("CRSnew_SimX_bigLib.RData")
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

source("functions_CRSnew.r")# load functions
SL.library = c("SL.randomForest","SL.glm","SL.gam", "SL.rpart","SL.knn","SL.gbm","SL.svm","SL.ipredbagg")

###############################################################
###############################################################
set.seed(2000)
# results r2
NN = (10^4)
nsize = 2*NN
kz = 4
set.seed(1)
z = rnorm(kz*nsize, mean=0, sd=1)
z = matrix(z, ncol=kz)

dimnames(z)[[2]] <- c("Z1", "Z2", "Z3", "Z4")


W = as.data.frame(z)

beta      = matrix( c(210, 27.4, 13.7, 13.7, 13.7), ncol=1 )
epsilon   = matrix( rnorm(nsize, mean=0, sd=100), ncol=1 )
YY        = cbind(1,z) %*% beta + epsilon

hist(YY)


###############################################################

temp1 = unique(sort(YY))
cut1 = floor(temp1[floor(nsize*0.7)])
cut1
b = c(27.4, 13.7, 13.7, 13.7)
sig = sqrt(sum(b^2)+100^2)
cut1 = qnorm(0.7, 210, sig)
cut1
Y = (YY >= cut1)*1
table(Y)
hist(YY)
abline(v = cut1)


###############################################################

# Fit the algorithm on the first half of the data
# Fit on (W1, Y1), get predictions on the whole simulated data (W, Y)
fit.data.SLL <- SuperLearner(Y=Y[1:NN], X=W[1:NN,], newX=W, SL.library = SL.library, family = binomial(), method = "method.NNLS",verbose = FALSE)
sl.pred <- fit.data.SLL$SL.predict #prediction from super learner 
lib.pred <- fit.data.SLL$library.predict #prediction from library algorithms
pred <- cbind(sl.pred,lib.pred) #all predictions
colnames(pred) <- c("SuperLearner","Random Forest","GLM","GAM","If-then Trees") # needs to change when SL.library varies
train.S <- pred[1:NN,]
val.S <- pred[(NN+1):(2*NN),]
train.CVS <- fit.data.SLL$Z
trainCVS_SL = fit.data.SLL$Z %*% fit.data.SLL$coef
train.CVS.wSL = cbind(trainCVS_SL, train.CVS) # add the first col to be the SL CV predictions

train.Z <- Y[1:NN]
val.Z <- Y[(NN+1):(2*NN)]


set.seed(1920413227)
# 3 algorithms, length(lambdas) lambdas
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
  crs.seed = as.integer(runif(1)*2e9)
  cbs = crs.fit(seed = crs.seed, lambda,train.CVS,train.Z, val.S[,-1], val.Z)
  cutoff[3,i] = cbs$c #CBS cutoff for lambda i 
  
  # SL Common
  opt = Opt.nonpar.rule(train.Z,train.S[,1],phi=0,lambda)
  cutoff[1,i] = as.numeric(opt)[1]
  # SL Proposal Iterative
  opt = Opt.nonpar.rule(train.Z,train.CVS.wSL[,1],phi=0,lambda)
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


epsig = 100

MY = cbind(1,z) %*% beta
S = pnorm(cut1, mean = MY, sd = epsig,lower.tail=FALSE)
S = pnorm(cut1-MY, mean = 0, sd = epsig,lower.tail=FALSE)
S2 = S[(NN+1):(2*NN)]

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

ind = c(5,9,13,17,21,25,29,33,37)
round(t(cbind(lambdas,draw)[ind,]),3)

trisk_U = draw
tlambda_U = lambdas
save(tlambda_U, trisk_U, file = "TCR_U1_bigLib.RData")
save.image("CRSnew_SimU_bigLib.RData")


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
library("ggplot2")


load("TCR_U1_bigLib.RData")
load("TCR_X1_bigLib.RData")

toplot = trisk_U
colnames(toplot) = c("Conditional Thresholding", "Two-Step Minimization", "CRS Minimization","Optimal Classification Rule")
lambdas = tlambda_U
dat1 = data.frame(cbind(lambdas,c(toplot)))
colnames(dat1) <- c("lambda","risk")
dat1$Approach = rep(colnames(toplot),rep(length(lambdas),ncol(toplot)))
p1 <- ggplot(dat1, aes(x=lambda, y=risk)) + geom_line(data=dat1,aes(x=lambda, y=risk,group=Approach,col = Approach),lwd=1) + ggtitle("Simulation 1")
p1 <- p1 + xlab(expression(paste(lambda)))
print(p1)

toplot = trisk_X
colnames(toplot) = c("Conditional Thresholding", "Two-Step Minimization", "CRS Minimization","Optimal Classification Rule")
lambdas = tlambda_X
dat2 = data.frame(cbind(lambdas,c(toplot)))
colnames(dat2) <- c("lambda","risk")
dat2$Approach = rep(colnames(toplot),rep(length(lambdas),ncol(toplot)))
p2 <- ggplot(dat2, aes(x=lambda, y=risk)) + geom_line(data=dat2,aes(x=lambda, y=risk,group=Approach,col = Approach),lwd=1) + ggtitle("Simulation 2")
p2 <- p2 + xlab(expression(paste(lambda)))
print(p2)


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


png(filename = "RiskPlot_Simulation_bigLib.png", width = 550, height = 650)

grid_arrange_shared_legend(p1, p2)

dev.off()

# make tables
ind = c(5,9,13,17,21,25,29,33,37)
round(t(cbind(lambdas,trisk_U)[ind,]),3)

ind = c(5,9,13,17,21,25,29,33,37)
round(t(cbind(lambdas,trisk_X)[ind,]),3)

ind = c(9,21,33)
round(t(cbind(lambdas,trisk_U)[ind,]),3)
ind = c(9,21,33)
round(t(cbind(lambdas,trisk_X)[ind,]),3)