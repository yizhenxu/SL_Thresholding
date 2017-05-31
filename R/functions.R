# X is cross validated predictions from single learners
# Xnew is predictions from single learners
isres.fit = function(lambda,X,Z,Xnew,Znew){
  # initial values
  mod1 = lm(Z ~ X)
  ord = order(abs(coef(mod1)[-1]),decreasing=TRUE) # order the columns by abs(coef)
  nna = sum(!is.na(coef(mod1)[-1]))
  ord = ord[1:nna]
  X.SL = X[,ord]
  mod1 = lm(Z ~ X.SL)
  b1 = coef(mod1)[3:length(coef(mod1))]/abs(coef(mod1)[2]) #initial b
  fc = as.numeric(sign(coef(mod1)[2])) #first  b coef
  G1 = fc*X.SL[,1] + X.SL[,-1]%*%b1
  cR1 = as.numeric(Opt.nonpar.rule(Z,G1,phi=0,lambda)[1]) #initial c
  r = range((predict(mod1) - coef(mod1)[1])/abs(coef(mod1)[2]))
  
  Rtest <- function(t, lambda, X, Y, fc) { 
    t = matrix(unlist(t),ncol=1) 
    G =  fc*X[,1]+X[,-1]%*%t[-1] # note that beta_gam is positive
    result = sum(lambda*(G<=t[1])*Y+(1-lambda)*(G>t[1])*(1-Y))
    return(result)
  }
  
  fn = function(t) Rtest(t, lambda, X.SL, Z, fc) # prepare for CRS (b,c)
  x0=c(cR1,b1)
  low = c(r[1]-0.5,rep(-5,length(b1)))
  upp = c(r[2]+0.5,rep(5,length(b1)))
  #crssol = crs2lm(x0 , fn , lower=low, upper=upp,
  #                maxeval = 10000, pop.size = 10000*(length(x0)+1), ranseed = seed,
  #                xtol_rel = 1e-6, nl.info = FALSE)
  isressol = isres(x0, fn, lower=low, upper=upp, hin = NULL, heq = NULL,
        maxeval = 10000, pop.size = 10000*(length(x0)+1),
        xtol_rel = 1e-6, nl.info = FALSE)
  bcrs = isressol$par[-1]
  ccrs = isressol$par[1]
  
  Xnew = as.matrix(Xnew[,ord])
  bcrs = matrix(bcrs,ncol=1)
  Gt = fc*Xnew[,1] + Xnew[,-1]%*%bcrs #crs predictions
  risk = sum(lambda*(Gt<=ccrs)*Znew+(1-lambda)*(Gt>ccrs)*(1-Znew)) 
  fnr = mean((Gt<=ccrs)*Znew)/mean(Znew)
  fpr = mean((Gt>ccrs)*(1-Znew))/mean(1-Znew)
  return(list(fc = fc, order = ord, c = ccrs, b = c(fc,bcrs), score = Gt, risk = risk, FPR = fpr, TPR = 1-fnr, sol=isressol))
  
}



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



hybridFun = function(lambda, phi, W, Y, maxrun=20){
  
  A = W[Y==1,]
  B = W[Y==0,]
  m = dim(A)[1]
  k = dim(B)[1]
  obj.fun = c(rep(lambda,m), rep(1-lambda,k),rep(0,2*ncol(W)))
  constr1 = cbind(diag(m), matrix(0,nrow=m,ncol=k), A, -A)
  constr2 = cbind(matrix(0,nrow=k,ncol=m), - diag(k), B, -B)
  constr = rbind(constr1, constr2)
  constr.dir = c(rep(">=",m), rep("<=",k))
  
  # initial cutoff c0 from lm full model
  mod3a = lm(Y ~ W)
  bh = bh0 = coef(mod3a)[-1] # b_0
  G = W%*%bh  
  OptR = Opt.nonpar.rule(Y,G,phi,lambda)
  g = g0 = OptR[5] # initial risk
  cut = cut0 = OptR[1:2] # c_0
  print(round(c(bh, cut, g),5))
  
  for(i in 1:maxrun){
    gprev = g
    bprev = round(bh,4)
    rhs = c(rep(cut[1]+1,m),rep(cut[2]-1,k)) # input c_(i-1)
    prod.sol = lp("min", obj.fun , constr , constr.dir , rhs ,compute.sens = TRUE )
    bh = prod.sol$solution[(m+k+1):(m+k+ncol(W))]-prod.sol$solution[(m+k+ncol(W)+1):(m+k+2*ncol(W))] #b_i
    
    if(all(bh == rep(0,length(bh)))){
      print("Null solution !")
      return(list(b = bh, cutoffs = c(NA,NA), risk = NA, bLS = bh0, cutoffLS = cut0, riskLS = g0))
    }
    
    G = W%*%bh
    OptR = Opt.nonpar.rule(Y,G,phi,lambda) 
    cut = OptR[1:2] # new cutoffs c_i
    g = OptR[5] # risk(b_i, c_i)
    print(i)
    print(round(c(bh, cut, g),5)) # b_i, c_i
    print(gprev)
    if(gprev >= g) {
      print(i)
      if(g < g0){
        cat("Solution is b = ",bh,", risk = ", g,"\n")
        return(list(b = bh, cutoffs = cut, risk = g, bLS = bh0, cutoffLS = cut0, riskLS = g0))
      } else {
        cat("No improvement, return the initial value","\n")
        return(list(b = bh0, cutoffs = cut0, risk = g0, bLS = bh0, cutoffLS = cut0, riskLS = g0))
      }
    }#if
  }#for
}#function



srgFun = function(lambda,  W, Y, Wnew,Ynew){
  
  A = W[Y==1,]
  B = W[Y==0,]
  m = dim(A)[1]
  k = dim(B)[1]
  obj.fun = c(rep(lambda,m), rep(1-lambda,k),rep(0,ncol(W)),0,0)
  constr1 = cbind(diag(m), matrix(0,nrow=m,ncol=k), A, matrix(-1,nrow=m,ncol=1),matrix(1,nrow=m,ncol=1) )
  constr2 = cbind(matrix(0,nrow=k,ncol=m), - diag(k), B, matrix(-1,nrow=k,ncol=1),matrix(1,nrow=k,ncol=1))
  constr3 = c(rep(0,m+k),rep(1,ncol(W)),0,0)
    constr = rbind(constr1, constr2, constr3)
  constr.dir = c(rep(">=",m), rep("<=",k), "==")
  rhs = c(rep(0,m),rep(0,k), 1)
  #rhs = c(rep(1,m),rep(-1,k), 1)
  prod.sol = lp("min", obj.fun , constr , constr.dir , rhs ,compute.sens = TRUE )
  bh = prod.sol$solution[(m+k+1):(m+k+ncol(W))]
  ch = prod.sol$solution[m+k+ncol(W)+1]-prod.sol$solution[m+k+ncol(W)+2]
  
  Wnew = as.matrix(Wnew)
  bh = matrix(bh,ncol=1)
  Gt = Wnew%*%bh #convex surrogate predictions
  risk = sum(lambda*(Gt<=ch)*Ynew+(1-lambda)*(Gt>ch)*(1-Ynew)) 
  fnr = mean((Gt<=ch)*Ynew)/mean(Ynew)
  fpr = mean((Gt>ch)*(1-Ynew))/mean(1-Ynew)
  return(list(c = ch, b = bh, score = Gt,  risk = risk, FPR = fpr, TPR = 1-fnr))
  
}#function


# X is cross validated predictions from single learners
# Xnew is predictions from single learners
crs.fit = function(seed,lambda,X,Z,Xnew,Znew){
  # initial values by non-negative least squares
  mod1 = nnls(as.matrix(X), Z) 
  initCoef = coef(mod1)
  initCoef[is.na(initCoef)] = 0
  #mod1 = lm(Z ~ X)
  ord = order(initCoef,decreasing=TRUE) # order the columns by initCoef
  X.SL = X[,ord]
  initCoef = initCoef[ord]
  b1 = initCoef/initCoef[1] #initial b
  G1 = X.SL%*%b1
  cR1 = as.numeric(Opt.nonpar.rule(Z,G1,phi=0,lambda)[1]) #initial c
  r = range(G1)
  
  Rtest <- function(t, lambda, X, Y) { #t[1] is cutoff, t[-1] is alpha
    t = matrix(unlist(t),ncol=1) 
    G =  X%*%t[-1] 
    result = sum(lambda*(G<=t[1])*Y+(1-lambda)*(G>t[1])*(1-Y))
    return(result)
  }
  
  fn = function(t) Rtest(t, lambda, X.SL, Z) # prepare for CRS (c,alpha)
  x0=c(cR1,b1)
  low = c(r[1]-0.5,rep(0,length(b1)))
  upp = c(r[2]+0.5,rep(5,length(b1)))
  crssol = crs2lm(x0 , fn , lower=low, upper=upp,
                  maxeval =  10000, pop.size = 100000*(length(x0)+1), ranseed = seed,
                  xtol_rel = 1e-6, nl.info = FALSE)
  bcrs = crssol$par[-1]
  norm = sum(bcrs)
  ccrs = crssol$par[1]/norm #normalize to make coefficients sum up to one
  bcrs = bcrs/norm
  
  Xnew = as.matrix(Xnew[,ord])
  bcrs = matrix(bcrs,ncol=1)
  Gt = Xnew%*%bcrs #crs predictions
  risk = mean(lambda*(Gt<=ccrs)*Znew+(1-lambda)*(Gt>ccrs)*(1-Znew)) 
  fnr = mean((Gt<=ccrs)*Znew)/mean(Znew)
  fpr = mean((Gt>ccrs)*(1-Znew))/mean(1-Znew)
  return(list(order = ord, c = ccrs, b = bcrs, score = Gt, risk = risk, FPR = fpr, TPR = 1-fnr, sol=crssol))
  
}




# Run these functions 

SL.bart <- function(Y, X, newX, family, ntree = 300, sigdf = 3, sigquant = 0.90, k = 2, power = 2, base = 0.95, binaryOffset = 0, ndpost = 1000, nskip = 100, ...) { 
  require('BayesTree')
  if(family$family == "gaussian") {
    fitBart <- bart(x.train = X, y.train = Y, x.test = newX, ntree = ntree, sigdf = sigdf, sigquant = sigquant, k = k, power = power, base = base, binaryOffset = binaryOffset, ndpost = ndpost, nskip = nskip, verbose = FALSE)
    pred  <- fitBart$yhat.test.mean
  }
  if(family$family == "binomial") {
    fitBart <- bart(x.train = X, y.train = as.factor(Y), x.test = newX, ntree = ntree, sigdf = sigdf, sigquant = sigquant, k = k, power = power, base = base, binaryOffset = binaryOffset, ndpost = ndpost, nskip = nskip, verbose = FALSE)
    pred <- pnorm(apply(fitBart$yhat.test, 2, mean))
  }
  fit <- list(object = fitBart)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.bart")
  return(out)
}

# 
predict.SL.bart <- function(object, newdata, ...) {
  stop("no predict method currently available for bart")
}

###########function for each round of cross validation############################

cal.each <- function(tZ,tS,vZ,vS,phi,lambda){
  
  n.opt.rule <- Opt.nonpar.rule(tZ,tS,phi,lambda)
  n.fnr.fpr <- nonpar.fnr.fpr(vZ,vS,n.opt.rule[1],n.opt.rule[2])
  
  s.opt.rule <- Opt.semipar.rule(tZ,tS,phi,lambda)
  s.fnr.fpr <- semipar.fnr.fpr(vZ,vS,s.opt.rule[1],s.opt.rule[2])
  
  p <- mean(vZ) 
  
  n.risk <- n.fnr.fpr[1]*p*lambda + n.fnr.fpr[2]*(1-p)*(1-lambda)
  n.TMR <- n.fnr.fpr[1]*p + n.fnr.fpr[2]*(1-p)
  
  s.risk <- s.fnr.fpr[1]*p*lambda + s.fnr.fpr[2]*(1-p)*(1-lambda)
  s.TMR <- s.fnr.fpr[1]*p + s.fnr.fpr[2]*(1-p)
  
  rules <- nonpar.rules(tZ,tS,phi)
  
  n.auc <- s.auc <- cal.AUC(vZ,vS,rules[,1],rules[,2])
  
  result <- c(n.fnr.fpr, n.risk, n.TMR, n.auc, s.fnr.fpr, s.risk, s.TMR, s.auc)
  names(result) <- c("nFNR","nFPR","nrisk","nTMR","nAUC","sFNR","sFPR","srisk","sTMR","sAUC")
  return(result)
}



###########function for making table 2 and 3############################
#NNLS

make.tables = function(phi, lambda, Z, yhat){ 
  
  ###calculate general l and u###
  temp <- ncol(yhat)  #number of algorithms
  r = matrix(NA,nrow=temp, ncol=4)
  for (i in 1:temp){
    r[i,1:2] = tryCatch(Opt.nonpar.rule(Z,yhat[,i],phi,lambda)[1:2], error=function(e) c(NA,NA))
    r[i,3:4] = tryCatch(Opt.semipar.rule(Z,yhat[,i],phi,lambda)[1:2], error=function(e) c(NA,NA))
  }
  rownames(r) = colnames(yhat)
  colnames(r) = c("nl","nu","sl","su")
  cutoffs <- round(r,3)
  ######################################################################################
  
  nFNR <- nFPR <- nAUC <- nTMR <- nrisk <- matrix(NA, ncol = temp, nrow = K)
  sFNR <- sFPR <- sAUC <- sTMR <- srisk <- matrix(NA, ncol = temp, nrow = K)
  
  K = 10
  
  set.seed(1)
  ind = sample(1:length(Z))
  fold = rep(1:10,table(cv.fold))
  fld.flg = fold[ind]
  
  for (i in 1:K){#loop through folds
    tZ <- Z[fld.flg!=i]
    tS <- yhat[fld.flg!=i,]
    vZ <- Z[fld.flg==i]
    vS <- yhat[fld.flg==i,]
    
    for (ia in 1:temp){#loop through algorithms
      result <- tryCatch(cal.each(tZ, tS[,ia], vZ, vS[,ia], phi, lambda), error=function(e) NULL)
      if(!is.null(result)){
        nFNR[i,ia] <- result[1]
        nFPR[i,ia] <- result[2]
        nrisk[i,ia] <- result[3]
        nTMR[i,ia] <- result[4]
        nAUC[i,ia] <- result[5]
        
        sFNR[i,ia] <- result[6]
        sFPR[i,ia] <- result[7]
        srisk[i,ia] <- result[8]
        sTMR[i,ia] <- result[9]
        sAUC[i,ia] <- result[10]
      }#if
    }#for ia
  }#for i
  
  outpt <- apply(cbind(nFNR,nFPR,nrisk,nTMR,nAUC,sFNR,sFPR,srisk,sTMR,sAUC),2,function(t) mean(t,na.rm=TRUE))
  o <- matrix(outpt,nrow=temp)
  colnames(o) <- c("nFNR","nFPR","nrisk","nTMR","nAUC","sFNR","sFPR","srisk","sTMR","sAUC")
  rownames(o) <- colnames(yhat)
  o <- round(o,3)
  
  ######################################################################################
  #make tables
  table2 <- as.data.frame(matrix(NA, ncol=7, nrow=temp))
  rownames(table2) <- colnames(yhat)
  colnames(table2) <- c("l",	"u",	"CV-FNR",	"CV-FPR",	"CV-risk",	"CV-TMR",	"CV-AUC")
  table2[,1:2] <- cutoffs[,1:2]
  table2[,3:7] <- o[,1:5]
  
  table3 <- as.data.frame(matrix(NA, ncol=7, nrow=temp))
  rownames(table3) <- colnames(yhat)
  colnames(table3) <- c("l",	"u",	"CV-FNR",	"CV-FPR",	"CV-risk",	"CV-TMR",	"CV-AUC")
  table3[,1:2] <- cutoffs[,3:4]
  table3[,3:7] <- o[,6:10]
  
  tab <- cbind(table2,table3)
  return(tab)
}

 

##########################TVLT functions###############################################

nonpar.rules <- function(Z,S,phi){
  if(length(Z)!=length(S))
    cat("***** Warning: Disease status and risk score vector lengths do not match. \n")
  data <- cbind(Z,S)
  Z <- data[complete.cases(data),1]
  S <- data[complete.cases(data),2]
  
  Z <- 1*Z #make logical values into {0,1}
  if(phi>1 || phi<0)
    cat("***** Warning: Invalid phi. \n")
  if(cor(Z,S)<0)
    cat("***** Warning: Disease status is negatively associated with risk score.
        Suggest using (-) value for risk score. \n")
  
  # "total.sam.sz" stands for total sample size
  total.sam.sz <- length(Z)
  # total unique sort risk.score
  S.srt <- sort(unique(S))
  # length(S)=645, length(risk.sc.srt)=457
  n.unique.S <- length(S.srt)
  # empirical cdf of S, cum.F=c(0,....,1)
  cum.F <- 0
  for(i in 1:n.unique.S){
    cum.F <- c(cum.F, mean(S<=S.srt[i],na.rm=TRUE))
  }
  # this can be also achieved by 
  # cdf.S=ecdf(S)
  # cum.F = c(0,cdf.S(S.srt))
  
  # identify all cut-off points that allow the percent in the middle
  # to be no more than phi
  
  # by cum.F[i+1]=P(S<=S.srt[i])=G(S.srt[i]) and cum.F[1]=0
  # the following code returns bounds=c(i,j), where
  # cum.F[j+1]-cum.F[i]=G(S.srt[j])-G(S.srt[i-1])
  # =P(S \in ( S.srt[i-1] , S.srt[j] ])=P(S \in [ S.srt[i] , S.srt[j] ])
  # <=\phi 
  flg <- T; i <- 1;j1 <- 0; bounds <- NULL
  while(flg){
    j <- sum((cum.F-cum.F[i]) <=phi)-1
    if(j>=i){
      if(i == 1){
        bounds <- c(i, j) 
        j1 <- j
      }
      if(i>1){
        if(j>j1){
          bounds <- rbind(bounds, c(i,j))
        }
        j1 <- j
      }
    }
    #c(i,j)
    i <- i+1
    if(i>n.unique.S || j1==n.unique.S) flg=F
  }
  if(phi==0){
    bounds <- cbind(1:n.unique.S, 1:n.unique.S)
    #cat("***** Warning: 0 patient taking viral load test. \n")
  }
  l <- S.srt[bounds[,1]]
  u <- S.srt[bounds[,2]]
  return(cbind(l,u))
}


nonpar.fnr.fpr <- function(Z,S,l,u){
  if(length(l)!=length(u))   #l is lower cutoff, u is upper cutoff
    cat("***** Warning: Wrong rules set. \n")
  #l is lower cutoff, u is upper cutoff
  n.bounds <-  length(l) #number of all possible rules
  mean.Z <- mean(Z,na.rm=TRUE)
  fnr.fpr <- NULL
  for(i in 1:n.bounds){
    fnr.fpr <- rbind(fnr.fpr, c(mean((S<l[i])*Z,na.rm=TRUE)/mean(Z,na.rm=TRUE),
                                mean((S>u[i])*(1-Z),na.rm=TRUE)/mean(1-Z,na.rm=TRUE)))
  }
  return(fnr.fpr)
}





semipar.fnr.fpr <- function(Z,S,l,u){
  if(length(l)!=length(u))   #l is lower cutoff, u is upper cutoff
    cat("***** Warning: Wrong rules set. \n")
  p <- mean(Z)
  temp <- density(S) #marginal density (S)
  fit <- glm(Z~ S, family=binomial)
  beta0star <- fit$coef[1]-log(p/(1-p))
  t <- exp(beta0star+temp$x*fit$coef[2]) #g1=t*g0 under exp tilt assumption
  g1 <- temp$y/(p+(1-p)/t)
  g0 <- temp$y/(p*t+1-p)
  x <- temp$x
  
  len <- length(x)
  dif <- x[2:len]-x[1:(len-1)]
  
  cal.fnr <- function(dens,a){  
    if( a>max(x) ){
      area <- 1
    } else if( a<min(x) ){
      area <- 0
    } else {
      diff <- a-x
      diff1 <- diff[diff<=0][1] 
      indx <- which(diff==diff1)#return index of nearest right endpoint
      area <- sum(dens[1:(indx-2)]*dif[1:(indx-2)])+dens[indx-1]*(a-x[indx-1])
    }
    return(area)
  } 
  
  fnr.fpr <- NULL
  K <- length(l)
  for( i in 1:K){
    fnr <- cal.fnr(g1,l[i])
    fpr <- 1-cal.fnr(g0,u[i])
    fnr.fpr <- rbind(fnr.fpr,c(fnr,fpr))        
  }
  
  return(fnr.fpr)
  
}





cal.AUC <- function(Z,S,l,u){
  ## AUC
  #Write the kth rule in Rule.set as (i_k,j_k), let j_0=0
  #Hphi(u)=argmin_w {G(u)-G(w)<=phi}
  #For Sj in (j_{k-1},j_k], Hphi(Sj)=i_k
  n = length(Z)
  p <- mean(Z)
  Hphi <- function(Sj,bounds=cbind(l,u)){
    diff <- Sj-bounds[,2]
    diff1 <- diff[diff<=0][1] #Sj-j_k, where Sj in (j_{k-1},j_k]
    indx <- which(diff==diff1)
    return(bounds[indx,1])
  }
  #calculate AUC from eqn (10) pg1177
  auc <- 0
  for(j in 1:n){
    auc <- auc+sum(Z*(1-Z[j])*((S>Hphi(S[j]))+(S==Hphi(S[j]))/2))
  }
  auc <- auc/(n^2*p*(1-p))  
}



Opt.nonpar.rule <- function(Z,S,phi,lambda){
  rules <- nonpar.rules(Z,S,phi)
  fnr.fpr <- nonpar.fnr.fpr(Z,S,rules[,1],rules[,2])
  fnr <- fnr.fpr[,1]
  fpr <- fnr.fpr[,2]
  p <- mean(Z,na.rm=TRUE)
  risk <- fnr*p*lambda + fpr*(1-p)*(1-lambda)
  index <- which.min(risk)
  opt.risk <- risk[index]
  opt.rule <- rules[index,]
  opt.fnr.fpr <- fnr.fpr[index,]
  TMR <- fnr.fpr[index,1]*p + fnr.fpr[index,2]*(1-p)
  z <- c(opt.rule,opt.fnr.fpr,opt.risk,TMR)
  names(z) <- c("lower.cutoff","upper.cutoff","FNR","FPR","opt.risk","TMR")
  return(z)
}




Opt.semipar.rule <- function(Z,S,phi,lambda){
  rules <- nonpar.rules(Z,S,phi)
  fnr.fpr <- semipar.fnr.fpr(Z,S,rules[,1],rules[,2])
  fnr <- fnr.fpr[,1]
  fpr <- fnr.fpr[,2]
  p <- mean(Z,na.rm=TRUE)
  risk <- fnr*p*lambda + fpr*(1-p)*(1-lambda)
  index <- which.min(risk)
  opt.risk <- risk[index]
  opt.rule <- rules[index,]
  opt.fnr.fpr <- fnr.fpr[index,]
  TMR <- fnr.fpr[index,1]*p + fnr.fpr[index,2]*(1-p)
  z <- c(opt.rule,opt.fnr.fpr,opt.risk,TMR)
  names(z) <- c("lower.cutoff","upper.cutoff","FNR","FPR","opt.risk","TMR")
  return(z)
}
