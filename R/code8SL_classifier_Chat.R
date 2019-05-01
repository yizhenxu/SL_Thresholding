
library("SuperLearner")
library("randomForest")
library("gam")
library("rpart")
library("dplyr")
library("plyr")
library("ggplot2")
library("nloptr")
library("lpSolve")
memory.limit(size=22264)
load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/Risk Score Wrap Up/04-18Classification_temp.RData")
load("C:/Users/Yizhen Xu/Google Drive/Desktop/2015 summer/Risk Score Wrap Up/04-22Classification_Chat.RData")

#training C
CTCLFPR = c()
CTCLFNR = c()
CTCLTPR = c()
CTCLrisk = c()
#ensemble C
CECLFPR = c()
CECLFNR = c()
CECLTPR = c()
CECLrisk = c()

for(i in 1:length(lambdas)){
  lambda = lambdas[i]
  for(k in 1:K){
    train.ind <- (cv.fold!=k)
    val.ind <- (cv.fold==k)
    CLSLki = get(paste("CLSL",k,"_", i, sep=""))
    sl.pred <-  CLSLki$SL.predict #prediction from super learner 
    train.CVSki = CLSLki$Z %*% CLSLki$coef
    val.Ski <- sl.pred[val.ind]
    train.Ski <- sl.pred[train.ind]
    val.Zki <- Z[val.ind]
    train.Zki <- Z[train.ind]
    
    CToptki = Opt.nonpar.rule(train.Zki,train.Ski,phi=0,lambda)
    CT = as.numeric(CToptki)[1]
    CEoptki = Opt.nonpar.rule(train.Zki,train.CVSki,phi=0,lambda)
    CE = as.numeric(CEoptki)[1]
    print(c(CT,CE))
    
    if(k == 1){
      CTCLCVSi = 1*(val.Ski > CT)
      CECLCVSi = 1*(val.Ski > CE)
      #CLCVSi = val.Ski
    } else {
      CTCLCVSi = c(CTCLCVSi, 1*(val.Ski > CT))
      CECLCVSi = c(CECLCVSi, 1*(val.Ski > CE))
      #CLCVSi = rbind(CLCVSi, val.Ski)
    }
  }# for k fold
  #CLCVYhati = (CLCVSi > 0.5)*1
  CTCLFPR[i] = mean((CTCLCVSi==1)*(1-CVZ))/mean(1-CVZ)
  CTCLFNR[i] = mean((CTCLCVSi==0)*(CVZ))/mean(CVZ)
  CTCLTPR[i] = 1-CTCLFNR[i]
  CTCLrisk[i] = lambda*mean(CVZ)*CTCLFNR[i] + (1-lambda)*(mean(1-CVZ))*CTCLFPR[i]
  CECLFPR[i] = mean((CECLCVSi==1)*(1-CVZ))/mean(1-CVZ)
  CECLFNR[i] = mean((CECLCVSi==0)*(CVZ))/mean(CVZ)
  CECLTPR[i] = 1-CECLFNR[i]
  CECLrisk[i] = lambda*mean(CVZ)*CECLFNR[i] + (1-lambda)*(mean(1-CVZ))*CECLFPR[i]
}#for i lambda

#CLrisk[,1] # classification risk of Super Learner
CTCLrisk
CECLrisk

colnames(CombineRisk)
dim(CombineRisk)
CombineRisk = cbind(CTCLrisk, CECLrisk,CombineRisk)
colnames(CombineRisk)[1:2] = c("CTSLClassifier","CESLClassifier")
colnames(CombineRisk)[12] = "EnsembleCSL"
look = c(1,2,3,7,12,4)
dat = data.frame(cbind(lambdas,c(CombineRisk[,look])))
colnames(dat) <- c("lambda","risk")
dat$rules = rep(colnames(CombineRisk)[look],rep(length(lambdas),length(look)))
p <- ggplot(dat, aes(x=lambda, y=risk)) + geom_line(data=dat,aes(x=lambda, y=risk,group=rules,col = rules),lwd=1)
print(p)


look = c(3,1,2,7:11,13:16,12,4)
ind = c(1,5,9,13,17,21,25,29,39,49,59)
toplot = CombineRisk[ind,look]
t(round(cbind(lambdas[ind],toplot),3))
look = c(1,2,3,7,12,4)
toplot = CombineRisk[,look]
colnames(toplot) = c("Decision Fusion Common", "Decision Fusion Proposal", "Decision Fusion c=0.5",
                     "Score Fusion Common", "Score Fusion Proposal Opt2", "Score Fusion Proposal Opt1")
dat = data.frame(cbind(lambdas,c(toplot)))
colnames(dat) <- c("lambda","risk")
dat$rules = rep(colnames(toplot),rep(length(lambdas),length(look)))
p <- ggplot(dat, aes(x=lambda, y=risk)) + geom_line(data=dat,aes(x=lambda, y=risk,group=rules,col = rules),lwd=1)
print(p)


plotRisk = cbind(lambdas,CombineRisk)
save(plotRisk,file = "04-26CombineRisk.RData")
#load("...")
lambdas = plotRisk[,1]
CombineRisk = plotRisk[,-1]

ind = c(1,5,9,13,17,21,25,29,39,49,59)
look = c(1,2,3,7,12,4)
toplot = CombineRisk[,look]
colnames(toplot) = c("Decision Fusion Common", "Decision Fusion Proposal", "Decision Fusion c=0.5",
                     "Score Fusion Common", "Score Fusion Proposal Opt2", "Score Fusion Proposal Opt1")
toplot = toplot[,4:6]
colnames(toplot) = c("SL Common Approach", "SL Proposal Iterative", "SL Proposal Controlled Search")
#toplot = toplot[,c(2,5)]
dat = data.frame(cbind(lambdas,c(toplot)))
colnames(dat) <- c("lambda","risk")
dat$rules = rep(colnames(toplot),rep(length(lambdas),ncol(toplot)))
p <- ggplot(dat, aes(x=lambda, y=risk)) + geom_line(data=dat,aes(x=lambda, y=risk,group=rules,col = rules),lwd=1)
print(p)






