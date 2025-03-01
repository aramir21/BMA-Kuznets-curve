#### Kuznets curve ####

################# Polynomial of degree 5 #########################

rm(list=ls())
dataI <- read.csv("regionalKuznets.csv", sep = ",", header = TRUE)
attach(dataI)
dataI$lnGPDpcXfederal <- federal*lnGDPpc
summary(dataI)
data <- na.omit(dataI)
summary(data)
attach(data)

set.seed(010101)
FixedEff <- as.matrix(fastDummies::dummy_cols(country)[,-1])
TimeEff <- cbind(t1, t2, t3, t4)
X <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3, lnGDPpc4, lnGDPpc5, trade, fdi, school, rents, gasoline,
           land, aid, areaXgasoline, ethnic_gini, lnGPDpcXfederal, polity2)
FixedReg <- lm(gini~FixedEff+TimeEff+X-1)
summary(FixedReg)
ResFixedReg <- summary(FixedReg)

FixedRegNew <- lm(gini~FixedEff+TimeEff+X[,-c(4:5)]-1)
summary(FixedRegNew)

# ############# LASSO ####################
# library(glmnet)
# Xlasso <- scale(cbind(FixedEff, TimeEff, X[,-c(4,5)]))
# cv_lasso <- cv.glmnet(Xlasso, gini, alpha = 1)
# # Plot cross-validation results
# plot(cv_lasso)
# # Optimal lambda
# best_lambda <- cv_lasso$lambda.min
# # Fit the model with the best lambda
# final_model <- glmnet(Xlasso, gini, alpha = 1, lambda = best_lambda)
# # Coefficients of the final model
# coef(final_model)

########################################
M <- expand.grid(c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0))
colnames(M) <- colnames(X)
n <- length(gini)

####################################### Equal model probability #########################################
############################### Fixed effects OLS Reg: Individual + Time ################################
logMargLikeFunctNormal <- function(Xr){
  Reg <- lm(gini~FixedEff+TimeEff+Xr)
  ns <- length(Reg$residuals)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}

logMargLikOLS <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~FixedEff+TimeEff)
    ns <- length(Reg$residuals)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS[m] <- logMargLike
  }else{
    logMargLikOLS[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+TimeEff+Xm-1)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS[m]-logMargLikOLS[which.max(logMargLikOLS)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
plot(ProbModelOrd)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLSfixedIndTim <- list(Top10 = Top10, PIP = PIP)
save(OLSfixedIndTim, file = "OLSfixedIndTimNewV2.RData")
# load(file = "OLSfixedIndTimNewV2.RData")

############################### Fixed effects OLS Reg: Individual ################################
logMargLikeFunctNormal <- function(Xr){
  Reg <- lm(gini~FixedEff+Xr)
  ns <- length(Reg$res)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}

logMargLikOLS1 <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~FixedEff)
    ns <- length(Reg$res)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS1[m] <- logMargLike
  }else{
    logMargLikOLS1[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS1),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+Xm-1)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS1[m]-logMargLikOLS1[which.max(logMargLikOLS1)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLSfixedInd <- list(Top10 = Top10, PIP = PIP)
save(OLSfixedInd, file = "OLSfixedIndNewV2.RData")
# load(file = "OLSfixedIndNewV2.RData")

############################### Fixed effects OLS Reg: Time ################################
logMargLikeFunctNormal <- function(Xr){
  Reg <- lm(gini~TimeEff+Xr)
  ns <- length(Reg$res)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}

logMargLikOLS2 <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~TimeEff)
    ns <- length(Reg$res)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS2[m] <- logMargLike
  }else{
    logMargLikOLS2[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS2),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+TimeEff+Xm-1)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS2[m]-logMargLikOLS2[which.max(logMargLikOLS2)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLSfixedTim <- list(Top10 = Top10, PIP = PIP)
save(OLSfixedTim, file = "OLSfixedTimNewV2.RData")
# load(file = "OLSfixedTimNewV2.RData")

############################### No fixed effects OLS Reg ################################

logMargLikeFunctNormal <- function(Xr){
  Reg <- lm(gini~Xr)
  ns <- length(Reg$res)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}
logMargLikeFunctNormal(X)

logMargLikOLS3 <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~1)
    ns <- length(Reg$res)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS3[m] <- logMargLike
  }else{
    logMargLikOLS3[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS3),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~Xm)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS3[m]-logMargLikOLS3[which.max(logMargLikOLS3)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLS <- list(Top10 = Top10, PIP = PIP)
save(OLS, file = "OLSnewV2.RData")
# load(file = "OLSnewV2.RData")

#####Summary PIP ######
summary <- cbind(OLSfixedIndTim$PIP, OLSfixedInd$PIP, OLSfixedTim$PIP, OLS$PIP)


###### Which fixed effects are better? ########
logMargLikOLSall <- c(logMargLikOLS, logMargLikOLS1, logMargLikOLS2, logMargLikOLS3) 
OddsRatMaxALL <- NULL
for(m in 1:(dim(M)[1]*4)){
  OddsRatMaxALL[m] <- exp(logMargLikOLSall[m]-logMargLikOLSall[which.max(logMargLikOLSall)])
}
summary(OddsRatMaxALL)
Pmax <- 1/sum(OddsRatMaxALL)
ProbModel <- OddsRatMaxALL*Pmax
save(ProbModel, file = "PMPall.RData")
summary(ProbModel)
summary(ProbModel[1:2^16]) # Two-way
summary(ProbModel[(2^16+1):(2*2^16)]) # Individual
summary(ProbModel[(2*2^16+1):(3*2^16)]) # Time
summary(ProbModel[(3*2^16+1):(4*2^16)]) # None
sum(ProbModel[1:2^16])
sum(ProbModel[(2^16+1):(2*2^16)]) # Best models involve fixed effects, but not time effects 
sum(ProbModel[(2*2^16+1):(3*2^16)])
sum(ProbModel[(3*2^16+1):(4*2^16)])
plot(ProbModel)
idXs <- which(M[which.max(logMargLikOLSall)-2^16,] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+Xm)
summary(RegMax)

ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:10
Top10 <- cbind(ProbModelOrd[nbest], M[idBestMo[nbest],])
sum(ProbModelOrd[nbest])

################# Polynomial of degree 4 #########################

dataI <- read.csv("regionalKuznets.csv", sep = ",", header = TRUE)
attach(dataI)
dataI$lnGPDpcXfederal <- federal*lnGDPpc
summary(dataI)
data <- na.omit(dataI)
summary(data)
attach(data)

set.seed(010101)
FixedEff <- as.matrix(fastDummies::dummy_cols(country)[,-1])
TimeEff <- cbind(t1, t2, t3, t4)
X <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3, lnGDPpc4, trade, fdi, school, rents, gasoline,
           land, aid, areaXgasoline, ethnic_gini, lnGPDpcXfederal, polity2)
FixedReg <- lm(gini~FixedEff+TimeEff+X-1)
summary(FixedReg)
ResFixedReg <- summary(FixedReg) 
M <- expand.grid(c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0))
colnames(M) <- colnames(X)
head(M)
n <- length(gini)

####################################### Equal model probability #########################################
############################### Fixed effects OLS Reg: Individual + Time ################################
logMargLikeFunctNormal <- function(Xr){
  # Xr <- X[,1]
  Reg <- lm(gini~FixedEff+TimeEff+Xr)
  ns <- length(Reg$residuals)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}

logMargLikOLS <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~FixedEff+TimeEff)
    ns <- length(Reg$residuals)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS[m] <- logMargLike
  }else{
    logMargLikOLS[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+TimeEff+Xm-1)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS[m]-logMargLikOLS[which.max(logMargLikOLS)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
plot(ProbModelOrd)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLSfixedIndTim <- list(Top10 = Top10, PIP = PIP)
save(OLSfixedIndTim, file = "OLSfixedIndTimNewV2p4.RData")
# load(file = "OLSfixedIndTimNewV2p4.RData")

############################### Fixed effects OLS Reg: Individual ################################
logMargLikeFunctNormal <- function(Xr){
  Reg <- lm(gini~FixedEff+Xr)
  ns <- length(Reg$res)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}

logMargLikOLS1 <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~FixedEff)
    ns <- length(Reg$res)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS1[m] <- logMargLike
  }else{
    logMargLikOLS1[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS1),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+Xm-1)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS1[m]-logMargLikOLS1[which.max(logMargLikOLS1)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLSfixedInd <- list(Top10 = Top10, PIP = PIP)
save(OLSfixedInd, file = "OLSfixedIndNewV2p4.RData")
# load(file = "OLSfixedIndNewV2.RData")

############################### Fixed effects OLS Reg: Time ################################
logMargLikeFunctNormal <- function(Xr){
  Reg <- lm(gini~TimeEff+Xr)
  ns <- length(Reg$res)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}

logMargLikOLS2 <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~TimeEff)
    ns <- length(Reg$res)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS2[m] <- logMargLike
  }else{
    logMargLikOLS2[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS2),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+TimeEff+Xm-1)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS2[m]-logMargLikOLS2[which.max(logMargLikOLS2)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLSfixedTim <- list(Top10 = Top10, PIP = PIP)
save(OLSfixedTim, file = "OLSfixedTimNewV2p4.RData")
# load(file = "OLSfixedTimNewV2p4.RData")

############################### No fixed effects OLS Reg ################################

logMargLikeFunctNormal <- function(Xr){
  Reg <- lm(gini~Xr)
  ns <- length(Reg$res)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}

logMargLikOLS3 <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~1)
    ns <- length(Reg$res)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS3[m] <- logMargLike
  }else{
    logMargLikOLS3[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS3),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~Xm)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS3[m]-logMargLikOLS3[which.max(logMargLikOLS3)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLS <- list(Top10 = Top10, PIP = PIP)
save(OLS, file = "OLSnewV2p4.RData")
# load(file = "OLSnewV2p4.RData")

#####Summary PIP ######
summary <- cbind(OLSfixedIndTim$PIP, OLSfixedInd$PIP, OLSfixedTim$PIP, OLS$PIP)


###### Which fixed effects are better? ########
logMargLikOLSall <- c(logMargLikOLS, logMargLikOLS1, logMargLikOLS2, logMargLikOLS3) 
OddsRatMaxALL <- NULL
for(m in 1:(dim(M)[1]*4)){
  OddsRatMaxALL[m] <- exp(logMargLikOLSall[m]-logMargLikOLSall[which.max(logMargLikOLSall)])
}
summary(OddsRatMaxALL)
Pmax <- 1/sum(OddsRatMaxALL)
ProbModel <- OddsRatMaxALL*Pmax
save(ProbModel, file = "PMPallp4.RData")
summary(ProbModel)
summary(ProbModel[1:2^15]) # Two-way
summary(ProbModel[(2^15+1):(2*2^15)]) # Individual
summary(ProbModel[(2*2^15+1):(3*2^15)]) # Time
summary(ProbModel[(3*2^15+1):(4*2^15)]) # None
sum(ProbModel[1:2^15])
sum(ProbModel[(2^15+1):(2*2^15)]) # Best models involve fixed effects, but not time effects 
sum(ProbModel[(2*2^15+1):(3*2^15)])
sum(ProbModel[(3*2^15+1):(4*2^15)])
plot(ProbModel)
idXs <- which(M[which.max(logMargLikOLSall)-2^15,] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+Xm-1)
summary(RegMax)

ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:10
Top10 <- cbind(ProbModelOrd[nbest], M[idBestMo[nbest],])
sum(ProbModelOrd[nbest])


################# Polynomial of degree 3 #########################
rm(list=ls())
dataI <- read.csv("regionalKuznets.csv", sep = ",", header = TRUE)

attach(dataI)
dataI$lnGPDpcXfederal <- federal*lnGDPpc
summary(dataI)
data <- na.omit(dataI)
summary(data)
attach(data)

set.seed(010101)
FixedEff <- as.matrix(fastDummies::dummy_cols(country)[,-1])
TimeEff <- cbind(t1, t2, t3, t4)
X <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3, trade, fdi, school, rents, gasoline,
           land, aid, areaXgasoline, ethnic_gini, lnGPDpcXfederal, polity2)
FixedReg <- lm(gini~FixedEff+TimeEff+X-1)
summary(FixedReg)
ResFixedReg <- summary(FixedReg) 
M <- expand.grid(c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0))
colnames(M) <- colnames(X)
head(M)
n <- length(gini)

####################################### Equal model probability #########################################
############################### Fixed effects OLS Reg: Individual + Time ################################
logMargLikeFunctNormal <- function(Xr){
  # Xr <- X[,1]
  Reg <- lm(gini~FixedEff+TimeEff+Xr)
  ns <- length(Reg$residuals)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}

logMargLikOLS <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~FixedEff+TimeEff)
    ns <- length(Reg$residuals)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS[m] <- logMargLike
  }else{
    logMargLikOLS[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+TimeEff+Xm-1)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS[m]-logMargLikOLS[which.max(logMargLikOLS)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
plot(ProbModelOrd)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLSfixedIndTim <- list(Top10 = Top10, PIP = PIP)
save(OLSfixedIndTim, file = "OLSfixedIndTimNewV2p3.RData")
# load(file = "OLSfixedIndTimNewV2p3.RData")

############################### Fixed effects OLS Reg: Individual ################################
logMargLikeFunctNormal <- function(Xr){
  Reg <- lm(gini~FixedEff+Xr)
  ns <- length(Reg$res)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}

logMargLikOLS1 <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~FixedEff)
    ns <- length(Reg$res)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS1[m] <- logMargLike
  }else{
    logMargLikOLS1[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS1),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+Xm-1)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS1[m]-logMargLikOLS1[which.max(logMargLikOLS1)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLSfixedInd <- list(Top10 = Top10, PIP = PIP)
save(OLSfixedInd, file = "OLSfixedIndNewV2p3.RData")
# load(file = "OLSfixedIndNewV2p3.RData")

############################### Fixed effects OLS Reg: Time ################################
logMargLikeFunctNormal <- function(Xr){
  Reg <- lm(gini~TimeEff+Xr)
  ns <- length(Reg$res)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}

logMargLikOLS2 <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~TimeEff)
    ns <- length(Reg$res)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS2[m] <- logMargLike
  }else{
    logMargLikOLS2[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS2),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+TimeEff+Xm-1)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS2[m]-logMargLikOLS2[which.max(logMargLikOLS2)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLSfixedTim <- list(Top10 = Top10, PIP = PIP)
save(OLSfixedTim, file = "OLSfixedTimNewV2p3.RData")
# load(file = "OLSfixedTimNewV2p3.RData")

############################### No fixed effects OLS Reg ################################

logMargLikeFunctNormal <- function(Xr){
  Reg <- lm(gini~Xr)
  ns <- length(Reg$res)
  ResReg <- summary(Reg) 
  k <- length(ResReg[["coefficients"]][,1]) 
  BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
  logMargLike <- -BIC/2
  return(logMargLike)
}

logMargLikOLS3 <- NULL
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~1)
    ns <- length(Reg$res)
    ResReg <- summary(Reg) 
    k <- length(ResReg[["coefficients"]][,1]) 
    BIC <- k*log(ns)+ns*log(ResReg[["sigma"]]^2*(ns-k)/ns)
    logMargLike <- -BIC/2
    logMargLikOLS3[m] <- logMargLike
  }else{
    logMargLikOLS3[m] <- logMargLikeFunctNormal(Xm)
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS3),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~Xm)
summary(RegMax)

OddsRatMax <- NULL
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS3[m]-logMargLikOLS3[which.max(logMargLikOLS3)])
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:100
Top10 <- cbind(ProbModelOrd[nbest], M[unlist(idBestMo[nbest]),])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd[1:10000]*M[unlist(idBestMo[1:10000]),])
OLS <- list(Top10 = Top10, PIP = PIP)
save(OLS, file = "OLSnewV2p3.RData")
# load(file = "OLSnewV2p3.RData")

#####Summary PIP ######
summary <- cbind(OLSfixedIndTim$PIP, OLSfixedInd$PIP, OLSfixedTim$PIP, OLS$PIP)


###### Which fixed effects are better? ########
logMargLikOLSall <- c(logMargLikOLS, logMargLikOLS1, logMargLikOLS2, logMargLikOLS3) 
OddsRatMaxALL <- NULL
for(m in 1:(dim(M)[1]*4)){
  OddsRatMaxALL[m] <- exp(logMargLikOLSall[m]-logMargLikOLSall[which.max(logMargLikOLSall)])
}
summary(OddsRatMaxALL)
Pmax <- 1/sum(OddsRatMaxALL)
ProbModel <- OddsRatMaxALL*Pmax
save(ProbModel, file = "PMPallp3.RData")
summary(ProbModel)
summary(ProbModel[1:2^14]) # Two-way
summary(ProbModel[(2^14+1):(2*2^14)]) # Individual
summary(ProbModel[(2*2^14+1):(3*2^14)]) # Time
summary(ProbModel[(3*2^14+1):(4*2^14)]) # None
sum(ProbModel[1:2^14])
sum(ProbModel[(2^14+1):(2*2^14)]) # Best models involve fixed effects, but not time effects 
sum(ProbModel[(2*2^14+1):(3*2^14)])
sum(ProbModel[(3*2^14+1):(4*2^14)])
plot(ProbModel)
idXs <- which(M[which.max(logMargLikOLSall)-2^14,] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~FixedEff+Xm)
summary(RegMax)

ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:10
Top10 <- cbind(ProbModelOrd[nbest], M[idBestMo[nbest],])
sum(ProbModelOrd[nbest])

############################################################################
########### Figures #######################
rm(list=ls())
dataI <- read.csv("regionalKuznets.csv", sep = ",", header = TRUE)
attach(dataI)
FixedEff <- as.matrix(fastDummies::dummy_cols(country)[,-1])
TimeEff <- cbind(t1, t2, t3, t4)
Xgdp5 <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3, lnGDPpc4, lnGDPpc5)

form <- gini~Xgdp5+FixedEff+TimeEff
RegBestp5 <- lm(form)
Resultado <- summary(RegBestp5)
print(Resultado)

Xn <- cbind(1, Xgdp5, FixedEff, TimeEff)
IdVarno <- which(is.na(RegBestp5[["coefficients"]])==1)
Xno <- Xn[,-IdVarno]

RegMCMC <- MCMCpack::MCMCregress(gini~Xno-1, mcmc = 6000, burnin = 1001, thin = 1)
BayesReg <- summary(RegMCMC)
thetamean <- BayesReg$statistics[,1]
XnoBar <- colMeans(Xno)
XnoBar[2:6] <- 0
Cte <- c(XnoBar%*%thetamean[-length(thetamean)])

####### Level ########
IdOrd <- order(lnGDPpc)
lnGDPpcOrd <- lnGDPpc[IdOrd]
plot(lnGDPpcOrd)
funPred <- function(theta){Cte + theta[1]*lnGDPpcOrd +theta[2]*lnGDPpcOrd^2+theta[3]*lnGDPpcOrd^3+theta[4]*lnGDPpcOrd^4+theta[5]*lnGDPpcOrd^5}
theta0 <- thetamean[2:6]
MeanPred <- funPred(theta0)
plot(MeanPred)

require(ggplot2)
df <- data.frame(cbind(gini, lnGDPpcOrd, MeanPred))

fig1 <- ggplot(df, aes(lnGDPpcOrd)) +                    # basic graphical object
  geom_point(aes(y=gini), colour="black") +  # first layer
  geom_line(aes(y=MeanPred), colour="blue", linewidth = 1.2) +  # second layer
  xlab("lnGDP") + ylab("Gini")
fig1

Xgdp4 <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3, lnGDPpc4)

form <- gini~Xgdp4+FixedEff+TimeEff
RegBestp4 <- lm(form)
Resultado <- summary(RegBestp4)
print(Resultado)

Xn <- cbind(1, Xgdp4, FixedEff, TimeEff)
IdVarno <- which(is.na(RegBestp4[["coefficients"]])==1)
Xno <- Xn[,-IdVarno]

RegMCMC <- MCMCpack::MCMCregress(gini~Xno-1, mcmc = 6000, burnin = 1001, thin = 1)
BayesReg <- summary(RegMCMC)
thetamean <- BayesReg$statistics[,1]
XnoBar <- colMeans(Xno)
XnoBar[2:5] <- 0
Cte <- c(XnoBar%*%thetamean[-length(thetamean)])

####### Level ########
funPred <- function(theta){Cte + theta[1]*lnGDPpcOrd +theta[2]*lnGDPpcOrd^2+theta[3]*lnGDPpcOrd^3+theta[4]*lnGDPpcOrd^4}
theta0 <- thetamean[2:5]
MeanPred4 <- funPred(theta0)

df <- data.frame(cbind(gini, lnGDPpcOrd, MeanPred4))

fig2 <- ggplot(df, aes(lnGDPpcOrd)) +                    # basic graphical object
  geom_point(aes(y=gini), colour="black") +  # first layer
  geom_line(aes(y=MeanPred4), colour="blue", linewidth = 1.2) +  # second layer
  xlab("lnGDP") + ylab("Gini")
fig2

Xgdp3 <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3)

form <- gini~Xgdp3+FixedEff+TimeEff
RegBestp3 <- lm(form)
Resultado <- summary(RegBestp3)
print(Resultado)

Xn <- cbind(1, Xgdp3, FixedEff, TimeEff)
IdVarno <- which(is.na(RegBestp3[["coefficients"]])==1)
Xno <- Xn[,-IdVarno]

RegMCMC <- MCMCpack::MCMCregress(gini~Xno-1, mcmc = 6000, burnin = 1001, thin = 1)
BayesReg <- summary(RegMCMC)
thetamean <- BayesReg$statistics[,1]
XnoBar <- colMeans(Xno)
XnoBar[2:4] <- 0
Cte <- c(XnoBar%*%thetamean[-length(thetamean)])

####### Level ########
funPred <- function(theta){Cte + theta[1]*lnGDPpcOrd +theta[2]*lnGDPpcOrd^2+theta[3]*lnGDPpcOrd^3}
theta0 <- thetamean[2:4]
MeanPred3 <- funPred(theta0)

df <- data.frame(cbind(gini, lnGDPpcOrd, MeanPred3))

fig3 <- ggplot(df, aes(lnGDPpcOrd)) +                    # basic graphical object
  geom_point(aes(y=gini), colour="black") +  # first layer
  geom_line(aes(y=MeanPred3), colour="blue", linewidth = 1.2) +  # second layer
  xlab("lnGDP") + ylab("Gini")
fig3

Xgdp2 <- cbind(lnGDPpc, lnGDPpc2)

form <- gini~Xgdp2+FixedEff+TimeEff
RegBestp2 <- lm(form)
Resultado <- summary(RegBestp2)
print(Resultado)

Xn <- cbind(1, Xgdp2, FixedEff, TimeEff)
IdVarno <- which(is.na(RegBestp2[["coefficients"]])==1)
Xno <- Xn[,-IdVarno]

RegMCMC <- MCMCpack::MCMCregress(gini~Xno-1, mcmc = 6000, burnin = 1001, thin = 1)
BayesReg <- summary(RegMCMC)
thetamean <- BayesReg$statistics[,1]
XnoBar <- colMeans(Xno)
XnoBar[2:3] <- 0
Cte <- c(XnoBar%*%thetamean[-length(thetamean)])

####### Level ########
funPred <- function(theta){Cte + theta[1]*lnGDPpcOrd +theta[2]*lnGDPpcOrd^2}
theta0 <- thetamean[2:3]
MeanPred2 <- funPred(theta0)

df <- data.frame(cbind(gini, lnGDPpcOrd, MeanPred2))

fig4 <- ggplot(df, aes(lnGDPpcOrd)) +                    # basic graphical object
  geom_point(aes(y=gini), colour="black") +  # first layer
  geom_line(aes(y=MeanPred2), colour="blue", linewidth = 1.2) +  # second layer
  xlab("lnGDP") + ylab("Gini")
fig4

library(ggpubr)
Allfig <- ggarrange(fig4, fig3, fig2, fig1, 
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
Allfig

######################### Coefficients ############################

######### Polynomial 5 ############

data <- read.csv("regionalKuznets.csv", sep = ",", header = TRUE)
attach(data)

FixedEff <- as.matrix(fastDummies::dummy_cols(country)[,-1])
TimeEff <- cbind(t1, t2, t3, t4)
X <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3, lnGDPpc4, lnGDPpc5)
M <- expand.grid(c(1,0), c(1,0), c(1,0), c(1,0), c(1,0))


############################### Fixed effects OLS Reg: Individual + Time ################################
logMargLikeFunctNormalCoef <- function(Xr){
  Reg <- lm(gini~Xr+FixedEff+TimeEff-1)
  if(is.null(dim(Xr))==1){
    l <- 1
    n <- length(Xr)
  }else{
    l <- dim(Xr)[2]
    n <- dim(Xr)[1]
  }
  Coef <- Reg$coefficients[1:l]
  ResReg <- summary(Reg)
  VarCoef <- diag(vcov(Reg))[1:l]
  k <- length(ResReg[["coefficients"]][,1])
  BIC <- k*log(n)+n*log(ResReg[["sigma"]]^2*(n-k)/n)
  logMargLike <- -BIC/2
  Rest <- list(Coef = Coef, VarCoef = VarCoef, logMargLike = logMargLike)
  return(Rest)
}

logMargLikOLS <- NULL
List.Coef <- vector(mode='list', length=2^15)
List.Var <- vector(mode='list', length=2^15)
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~FixedEff+TimeEff)
    Coef <- 0
    ResReg <- summary(Reg)
    VarCoef <- 0
    k <- length(ResReg[["coefficients"]][,1]) 
    n <- dim(Xm)[1]
    BIC <- k*log(n)+n*log(ResReg[["sigma"]]^2*(n-k)/n)
    logMargLike <- -BIC/2
    logMargLikOLS[m] <- logMargLike
    List.Coef[[m]]<- 0
    List.Var[[m]] <- 0
  }else{
    Result <- logMargLikeFunctNormalCoef(Xm)
    logMargLikOLS[m] <- Result$logMargLike
    List.Coef[[m]] <- Result$Coef
    List.Var[[m]] <- Result$VarCoef
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~Xm+FixedEff+TimeEff-1)
summary(RegMax)

OddsRatMax <- NULL
MatCoef <- matrix(0,dim(M)[1],dim(M)[2])
MatVar <- matrix(0,dim(M)[1],dim(M)[2])
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS[m]-logMargLikOLS[which.max(logMargLikOLS)])
  idG <- which(M[m,]==1)
  MatCoef[m,idG] <- List.Coef[[m]]
  MatVar[m,idG] <- List.Var[[m]]
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:10
Top10 <- cbind(ProbModelOrd[nbest], M[idBestMo[nbest],])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd*M[idBestMo,])
CoefMean <- colSums(ProbModel*MatCoef)
CoefVari <- matrix(NA, dim(M)[1],dim(M)[2])
for(i in 1:dim(M)[1]){
  CoefVari[i, ] <- ProbModel[i]*MatVar[i,]+ProbModel[i]*((MatCoef[i,]-CoefMean)^2) 
}
CoefVar <- colSums(CoefVari)
CoefMean/CoefVar^0.5
OLSfixedIndTim <- list(Top10 = Top10, PIP = PIP)
cbind(PIP, CoefMean, CoefVar^0.5, c(Top10[1,-1]))
save(OLSfixedIndTim, file = "OLSfixedIndTimNewV1p5Coeff.RData")
# load(file = "OLSfixedIndTimNewV1p5Coeff.RData")

######### Polynomial 4 ############

data <- read.csv("regionalKuznets.csv", sep = ",", header = TRUE)
attach(data)

FixedEff <- as.matrix(fastDummies::dummy_cols(country)[,-1])
TimeEff <- cbind(t1, t2, t3, t4)
X <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3, lnGDPpc4)
M <- expand.grid(c(1,0), c(1,0), c(1,0), c(1,0))

############################### Fixed effects OLS Reg: Individual + Time ################################
logMargLikeFunctNormalCoef <- function(Xr){
  Reg <- lm(gini~Xr+FixedEff+TimeEff-1)
  if(is.null(dim(Xr))==1){
    l <- 1
    n <- length(Xr)
  }else{
    l <- dim(Xr)[2]
    n <- dim(Xr)[1]
  }
  Coef <- Reg$coefficients[1:l]
  ResReg <- summary(Reg)
  VarCoef <- diag(vcov(Reg))[1:l]
  k <- length(ResReg[["coefficients"]][,1])
  BIC <- k*log(n)+n*log(ResReg[["sigma"]]^2*(n-k)/n)
  logMargLike <- -BIC/2
  Rest <- list(Coef = Coef, VarCoef = VarCoef, logMargLike = logMargLike)
  return(Rest)
}

logMargLikOLS <- NULL
List.Coef <- vector(mode='list', length=2^15)
List.Var <- vector(mode='list', length=2^15)
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~FixedEff+TimeEff)
    Coef <- 0
    ResReg <- summary(Reg)
    VarCoef <- 0
    k <- length(ResReg[["coefficients"]][,1]) 
    n <- dim(Xm)[1]
    BIC <- k*log(n)+n*log(ResReg[["sigma"]]^2*(n-k)/n)
    logMargLike <- -BIC/2
    logMargLikOLS[m] <- logMargLike
    List.Coef[[m]]<- 0
    List.Var[[m]] <- 0
  }else{
    Result <- logMargLikeFunctNormalCoef(Xm)
    logMargLikOLS[m] <- Result$logMargLike
    List.Coef[[m]] <- Result$Coef
    List.Var[[m]] <- Result$VarCoef
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~Xm+FixedEff+TimeEff-1)
summary(RegMax)

OddsRatMax <- NULL
MatCoef <- matrix(0,dim(M)[1],dim(M)[2])
MatVar <- matrix(0,dim(M)[1],dim(M)[2])
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS[m]-logMargLikOLS[which.max(logMargLikOLS)])
  idG <- which(M[m,]==1)
  MatCoef[m,idG] <- List.Coef[[m]]
  MatVar[m,idG] <- List.Var[[m]]
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:10
Top10 <- cbind(ProbModelOrd[nbest], M[idBestMo[nbest],])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd*M[idBestMo,])
CoefMean <- colSums(ProbModel*MatCoef)
CoefVari <- matrix(NA, dim(M)[1],dim(M)[2])
for(i in 1:dim(M)[1]){
  CoefVari[i, ] <- ProbModel[i]*MatVar[i,]+ProbModel[i]*((MatCoef[i,]-CoefMean)^2) 
}
CoefVar <- colSums(CoefVari)
CoefMean/CoefVar^0.5
OLSfixedIndTim <- list(Top10 = Top10, PIP = PIP)
cbind(PIP, CoefMean, CoefVar^0.5, c(Top10[1,-1]))
save(OLSfixedIndTim, file = "OLSfixedIndTimNewV1p4Coeff.RData")
# load(file = "OLSfixedIndTimNewV1p5Coeff.RData")

######### Polynomial 3 ############

data <- read.csv("regionalKuznets.csv", sep = ",", header = TRUE)
attach(data)

FixedEff <- as.matrix(fastDummies::dummy_cols(country)[,-1])
TimeEff <- cbind(t1, t2, t3, t4)
X <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3)
M <- expand.grid(c(1,0), c(1,0), c(1,0))

############################### Fixed effects OLS Reg: Individual + Time ################################
logMargLikeFunctNormalCoef <- function(Xr){
  Reg <- lm(gini~Xr+FixedEff+TimeEff-1)
  if(is.null(dim(Xr))==1){
    l <- 1
    n <- length(Xr)
  }else{
    l <- dim(Xr)[2]
    n <- dim(Xr)[1]
  }
  Coef <- Reg$coefficients[1:l]
  ResReg <- summary(Reg)
  VarCoef <- diag(vcov(Reg))[1:l]
  k <- length(ResReg[["coefficients"]][,1])
  BIC <- k*log(n)+n*log(ResReg[["sigma"]]^2*(n-k)/n)
  logMargLike <- -BIC/2
  Rest <- list(Coef = Coef, VarCoef = VarCoef, logMargLike = logMargLike)
  return(Rest)
}
logMargLikeFunctNormalCoef(X)

logMargLikOLS <- NULL
List.Coef <- vector(mode='list', length=2^15)
List.Var <- vector(mode='list', length=2^15)
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~FixedEff+TimeEff)
    Coef <- 0
    ResReg <- summary(Reg)
    VarCoef <- 0
    k <- length(ResReg[["coefficients"]][,1]) 
    n <- dim(Xm)[1]
    BIC <- k*log(n)+n*log(ResReg[["sigma"]]^2*(n-k)/n)
    logMargLike <- -BIC/2
    logMargLikOLS[m] <- logMargLike
    List.Coef[[m]]<- 0
    List.Var[[m]] <- 0
  }else{
    Result <- logMargLikeFunctNormalCoef(Xm)
    logMargLikOLS[m] <- Result$logMargLike
    List.Coef[[m]] <- Result$Coef
    List.Var[[m]] <- Result$VarCoef
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~Xm+FixedEff-1)
summary(RegMax)

OddsRatMax <- NULL
MatCoef <- matrix(0,dim(M)[1],dim(M)[2])
MatVar <- matrix(0,dim(M)[1],dim(M)[2])
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS[m]-logMargLikOLS[which.max(logMargLikOLS)])
  idG <- which(M[m,]==1)
  MatCoef[m,idG] <- List.Coef[[m]]
  MatVar[m,idG] <- List.Var[[m]]
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:8
Top10 <- cbind(ProbModelOrd[nbest], M[idBestMo[nbest],])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd*M[idBestMo,])
CoefMean <- colSums(ProbModel*MatCoef)
CoefVari <- matrix(NA, dim(M)[1],dim(M)[2])
for(i in 1:dim(M)[1]){
  CoefVari[i, ] <- ProbModel[i]*MatVar[i,]+ProbModel[i]*((MatCoef[i,]-CoefMean)^2) 
}
CoefVar <- colSums(CoefVari)
CoefMean/CoefVar^0.5
OLSfixedIndTim <- list(Top10 = Top10, PIP = PIP)
cbind(PIP, CoefMean, CoefVar^0.5, c(Top10[1,-1]))
save(OLSfixedIndTim, file = "OLSfixedIndTimNewV1p3Coeff.RData")
# load(file = "OLSfixedIndTimNewV1p3Coeff.RData")

######### Polynomial 2 ############

data <- read.csv("regionalKuznets.csv", sep = ",", header = TRUE)
attach(data)

FixedEff <- as.matrix(fastDummies::dummy_cols(country)[,-1])
# TimeEffex <- as.matrix(fastDummies::dummy_cols(year)[,-1])
TimeEff <- cbind(t1, t2, t3, t4)
# DataF <- cbind(data, FixedEff)
# write.csv(DataF, file = "DataF.csv")
X <- cbind(lnGDPpc, lnGDPpc2)
M <- expand.grid(c(1,0), c(1,0))
# colnames(M) <- colnames(X)


############################### Fixed effects OLS Reg: Individual + Time ################################
logMargLikeFunctNormalCoef <- function(Xr){
  Reg <- lm(gini~Xr+FixedEff+TimeEff-1)
  if(is.null(dim(Xr))==1){
    l <- 1
    n <- length(Xr)
  }else{
    l <- dim(Xr)[2]
    n <- dim(Xr)[1]
  }
  Coef <- Reg$coefficients[1:l]
  ResReg <- summary(Reg)
  VarCoef <- diag(vcov(Reg))[1:l]
  k <- length(ResReg[["coefficients"]][,1])
  BIC <- k*log(n)+n*log(ResReg[["sigma"]]^2*(n-k)/n)
  logMargLike <- -BIC/2
  Rest <- list(Coef = Coef, VarCoef = VarCoef, logMargLike = logMargLike)
  return(Rest)
}
logMargLikeFunctNormalCoef(X)

logMargLikOLS <- NULL
List.Coef <- vector(mode='list', length=2^15)
List.Var <- vector(mode='list', length=2^15)
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~FixedEff+TimeEff)
    Coef <- 0
    ResReg <- summary(Reg)
    VarCoef <- 0
    k <- length(ResReg[["coefficients"]][,1]) 
    n <- dim(Xm)[1]
    BIC <- k*log(n)+n*log(ResReg[["sigma"]]^2*(n-k)/n)
    logMargLike <- -BIC/2
    logMargLikOLS[m] <- logMargLike
    List.Coef[[m]]<- 0
    List.Var[[m]] <- 0
  }else{
    Result <- logMargLikeFunctNormalCoef(Xm)
    logMargLikOLS[m] <- Result$logMargLike
    List.Coef[[m]] <- Result$Coef
    List.Var[[m]] <- Result$VarCoef
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~Xm+FixedEff-1)
summary(RegMax)

OddsRatMax <- NULL
MatCoef <- matrix(0,dim(M)[1],dim(M)[2])
MatVar <- matrix(0,dim(M)[1],dim(M)[2])
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS[m]-logMargLikOLS[which.max(logMargLikOLS)])
  idG <- which(M[m,]==1)
  MatCoef[m,idG] <- List.Coef[[m]]
  MatVar[m,idG] <- List.Var[[m]]
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:4
Top10 <- cbind(ProbModelOrd[nbest], M[idBestMo[nbest],])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd*M[idBestMo,])
CoefMean <- colSums(ProbModel*MatCoef)
CoefVari <- matrix(NA, dim(M)[1],dim(M)[2])
for(i in 1:dim(M)[1]){
  CoefVari[i, ] <- ProbModel[i]*MatVar[i,]+ProbModel[i]*((MatCoef[i,]-CoefMean)^2) 
}
CoefVar <- colSums(CoefVari)
CoefMean/CoefVar^0.5
OLSfixedIndTim <- list(Top10 = Top10, PIP = PIP)
cbind(PIP, CoefMean, CoefVar^0.5, c(Top10[1,-1]))
save(OLSfixedIndTim, file = "OLSfixedIndTimNewV1p2Coeff.RData")
# load(file = "OLSfixedIndTimNewV1p2Coeff.RData")


#######################################################
# load(file = "OLSfixedIndTimNewV2p3.RData")
dataI <- read.csv("regionalKuznets.csv", sep = ",", header = TRUE)
data <- na.omit(dataI)
attach(data)

FixedEff <- as.matrix(fastDummies::dummy_cols(country)[,-1])
# TimeEffex <- as.matrix(fastDummies::dummy_cols(year)[,-1])
TimeEff <- cbind(t1, t2, t3, t4)
# DataF <- cbind(data, FixedEff)
# write.csv(DataF, file = "DataF.csv")
X <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3, trade, fdi, school, rents, gasoline,
           land, aid, areaXgasoline, ethnic_gini, lnGPDpcXfederal, polity2)
M <- expand.grid(c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0), c(1,0))
colnames(M) <- colnames(X)


############################### Fixed effects OLS Reg: Individual + Time ################################
logMargLikeFunctNormalCoef <- function(Xr){
  Reg <- lm(gini~Xr+FixedEff+TimeEff-1)
  if(is.null(dim(Xr))==1){
    l <- 1
    n <- length(Xr)
  }else{
    l <- dim(Xr)[2]
    n <- dim(Xr)[1]
  }
  Coef <- Reg$coefficients[1:l]
  ResReg <- summary(Reg)
  VarCoef <- diag(vcov(Reg))[1:l]
  k <- length(ResReg[["coefficients"]][,1])
  BIC <- k*log(n)+n*log(ResReg[["sigma"]]^2*(n-k)/n)
  logMargLike <- -BIC/2
  Rest <- list(Coef = Coef, VarCoef = VarCoef, logMargLike = logMargLike)
  return(Rest)
}
logMargLikeFunctNormalCoef(X)

logMargLikOLS <- NULL
List.Coef <- vector(mode='list', length=2^15)
List.Var <- vector(mode='list', length=2^15)
for(m in 1:dim(M)[1]){
  idXs <- which(M[m,] == 1)
  Xm <- X[,idXs]
  if(length(idXs)==0){
    Reg <- lm(gini~FixedEff+TimeEff)
    Coef <- 0
    ResReg <- summary(Reg)
    VarCoef <- 0
    k <- length(ResReg[["coefficients"]][,1]) 
    n <- dim(Xm)[1]
    BIC <- k*log(n)+n*log(ResReg[["sigma"]]^2*(n-k)/n)
    logMargLike <- -BIC/2
    logMargLikOLS[m] <- logMargLike
    List.Coef[[m]]<- 0
    List.Var[[m]] <- 0
  }else{
    Result <- logMargLikeFunctNormalCoef(Xm)
    logMargLikOLS[m] <- Result$logMargLike
    List.Coef[[m]] <- Result$Coef
    List.Var[[m]] <- Result$VarCoef
  }
  print(m)
}

idXs <- which(M[which.max(logMargLikOLS),] == 1)
Xm <- X[,idXs]
RegMax <- lm(gini~Xm+FixedEff+TimeEff-1)
summary(RegMax)

OddsRatMax <- NULL
MatCoef <- matrix(0,dim(M)[1],dim(M)[2])
MatVar <- matrix(0,dim(M)[1],dim(M)[2])
for(m in 1:dim(M)[1]){
  OddsRatMax[m] <- exp(logMargLikOLS[m]-logMargLikOLS[which.max(logMargLikOLS)])
  idG <- which(M[m,]==1)
  MatCoef[m,idG] <- List.Coef[[m]]
  MatVar[m,idG] <- List.Var[[m]]
}
summary(OddsRatMax)
Pmax <- 1/sum(OddsRatMax)
ProbModel <- OddsRatMax*Pmax 
summary(ProbModel)
plot(ProbModel)
ProbModelOrd <- sort(ProbModel, decreasing = TRUE)
idBestMo <- sapply(1:length(ProbModelOrd), function(i){which(ProbModel == ProbModelOrd[i])})
nbest <- 1:10
Top10 <- cbind(ProbModelOrd[nbest], M[idBestMo[nbest],])
sum(ProbModelOrd[nbest])
PIP <- colSums(ProbModelOrd*M[idBestMo,])
CoefMean <- colSums(ProbModel*MatCoef)
CoefVari <- matrix(NA, dim(M)[1],dim(M)[2])
for(i in 1:dim(M)[1]){
  CoefVari[i, ] <- ProbModel[i]*MatVar[i,]+ProbModel[i]*((MatCoef[i,]-CoefMean)^2) 
}
CoefVar <- colSums(CoefVari)
CoefMean/CoefVar^0.5
OLSfixedIndTim <- list(Top10 = Top10, PIP = PIP)
cbind(PIP, CoefMean, CoefVar^0.5, c(Top10[1,-1]))
save(OLSfixedIndTim, file = "OLSfixedIndTimNewV1p3Coeff.RData")
# load(file = "OLSfixedIndTimNewV1p3Coeff.RData")

####### Jointness #######
Ms <- 1:2^14
Ind3 <- ifelse(M[,3]== 1, 1, 0)
PMPX3 <- sum(ProbModel[which(Ind3 == 1)])
PMPX3

Ind2n3 <- ifelse(M[,2]== 1 & M[,3]== 0, 1, 0)
PMPX2n3 <- sum(ProbModel[which(Ind2n3 == 1)])
PMPX2n3

Ind23 <- ifelse(M[,2]== 1 & M[,3]== 1, 1, 0)
PMPX23 <- sum(ProbModel[which(Ind23 == 1)])
PMPX23

Indn23 <- ifelse(M[,2]== 0 & M[,3]== 1, 1, 0)
PMPXn23 <- sum(ProbModel[which(Indn23 == 1)])
PMPXn23

Ind23 <- ifelse(M[,2]== 1 & M[,3]== 1, 1, 0)
PMPX23 <- sum(ProbModel[which(Ind23 == 1)])
PMPX23

Ind12 <- ifelse(M[,1] == 1 & M[,2]== 1, 1, 0)
PMPX12 <- sum(ProbModel[which(Ind12 == 1)])

Ind1n2 <- ifelse(M[,1] == 1 & M[,2]== 0, 1, 0)
PMPX1n2 <- sum(ProbModel[which(Ind1n2 == 1)])

Indn1n2 <- ifelse(M[,1] == 0 & M[,2]== 0, 1, 0)
PMPXn1n2 <- sum(ProbModel[which(Indn1n2 == 1)])

Indn12 <- ifelse(M[,1] == 0 & M[,2]== 1, 1, 0)
PMPXn12 <- sum(ProbModel[which(Indn12 == 1)])

Ind123 <- ifelse(M[,1] == 1 & M[,2]== 1 & M[,3] == 1, 1, 0)
Ind123T <- which(Ind123 == 1)
PPM123 <- ProbModel[Ind123T]
PMPX123 <- sum(PPM123)
PMPX123Max <- max(PPM123)

Ind12n3 <- ifelse(M[,1] == 1 & M[,2]== 1 & M[,3] == 0, 1, 0)
Ind12n3T <- which(Ind12n3 == 1)
PPM12n3 <- ProbModel[Ind12n3T]
PMPX12n3 <- sum(PPM12n3)
PMPX12n3Max <- max(PPM12n3)



######### Best model p3 ##########
dataI <- read.csv("regionalKuznets.csv", sep = ",", header = TRUE)
attach(dataI)

FixedEff <- as.matrix(fastDummies::dummy_cols(country)[,-1])
TimeEff <- cbind(t1, t2, t3, t4)
X <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3, trade, fdi, school, rents, gasoline,
           land, aid, areaXgasoline, ethnic_gini, lnGPDpcXfederal, polity2)

BestRegMPM <- which(OLSfixedIndTim[["PIP"]]>= 0.5)
XMPM <- cbind(lnGDPpc, lnGDPpc2, lnGDPpc3, rents, land, ethnic_gini)
BestRegHPM <- which(OLSfixedIndTim[["Top10"]][1,-1] == 1)
logMar3p <- logMargLikeFunctNormal(XMPM)
logMar2p <- logMargLikeFunctNormal(XMPM[,-3])
Oddsp2p3 <- exp(logMar2p-logMar3p)
Pp3 <- 1/sum(1+Oddsp2p3)
Pp2 <- Oddsp2p3/sum(1+Oddsp2p3)
Pp2; Pp3

RegXMPM <- lm(gini~XMPM+FixedEff+TimeEff-1)
RegXMPMf <- summary(RegXMPM)

Xno <- XMPM[-RegXMPM[["na.action"]],]
FixedEffno <- FixedEff[-RegXMPM[["na.action"]],]
TimeEffno <- TimeEff[-RegXMPM[["na.action"]],]

Xn <- cbind(Xno, FixedEffno, TimeEffno)
IdVarno <- which(is.na(RegXMPM[["coefficients"]])==1)
Xno <- Xn[,-IdVarno]
Yn <- gini[-RegXMPM[["na.action"]]]
solve(t(Xno)%*%Xno)%*%t(Xno)%*%Yn

RegMCMC <- MCMCpack::MCMCregress(Yn~Xno-1, mcmc = 6000, burnin = 1001, thin = 1)
BayesReg <- summary(RegMCMC)
plot(RegMCMC)
coda::raftery.diag(RegMCMC)

thetamean <- BayesReg[["statistics"]][,1]
XnoBar <- colMeans(Xno)
cbind(thetamean[-length(thetamean)], XnoBar)
XnoBar[1:3] <- 0
Cte <- c(XnoBar%*%thetamean[-length(thetamean)])


######## Level #######
IdOrd <- order(lnGDPpc)
lnfig <- c(lnGDPpc[IdOrd], seq(max(lnGDPpc[IdOrd]), 13, 0.1))
funPred <- function(theta){theta[1]+theta[2]*lnfig+theta[3]*lnfig^2+theta[4]*lnfig^3}
theta0 <- c(Cte, thetamean[1:3]) 
MeanPred <- funPred(theta0)
# plot(lnGDPpc[IdOrd], funPred(theta0), type = "l", col = "blue")

# theta0[2]+2*theta0[3]*8.655+3*theta0[4]*8.655^2
# theta0[1]+2*theta0[2]*10.381+3*theta0[3]*10.381^2
# 8.655 maximizes inequality 
# 10.381 minimizes inequality

######## Derivative #######
funPredDer <- function(theta){theta[1]+2*theta[2]*lnfig+3*theta[3]*lnfig^2}
theta0 <- thetamean[1:3]
MeanPredDer <- funPredDer(theta0)

funPredDer2 <- function(theta){2*theta[1]+6*theta[2]*lnfig}
theta0 <- thetamean[2:3]
MeanPredDer2 <- funPredDer2(theta0)

funPredDerPIB <- function(pib){theta0[1]+2*theta0[2]*pib+3*theta0[3]*pib^2}
theta0 <- thetamean[1:3]
pibMin <- (-2*theta0[2] - ((2*theta0[2])^2-4*3*theta0[3]*theta0[1])^0.5)/(2*3*theta0[3]) 
pibMax <- (-2*theta0[2] + ((2*theta0[2])^2-4*3*theta0[3]*theta0[1])^0.5)/(2*3*theta0[3]) 

######## Frequentist variability of Bayesian estimators #######
##### Optimal: Min and Max
fdelta <- function(B){
  a <- 3*B[3]
  b <- 2*B[2]
  c <- B[1]
  pibMin_f <- (- b - (b^2-4*a*c)^0.5)/(2*a)
  pibMax_f <- (- b + (b^2-4*a*c)^0.5)/(2*a)
  ResFun <- list(pibMin_f = pibMin_f, pibMax_f = pibMax_f)
  return(ResFun)
} 
theta0 <- thetamean[1:3]
fdelta(theta0)

fdeltaMin <- function(B){
  a <- 3*B[3]
  b <- 2*B[2]
  c <- B[1]
  pibMin_f <- (- b - (b^2-4*a*c)^0.5)/(2*a)
  return(pibMin_f)
} 
ttnewMin <- apply(RegMCMC[,1:3], 1, fdeltaMin)

BETA <- BayesReg[["statistics"]][,1][-length(BayesReg[["statistics"]][,1])] 
sigma2 <- BayesReg[["statistics"]][,1][length(BayesReg[["statistics"]][,1])] 
v <- dim(Xno)[1]-dim(Xno)[2]
falpha <- function(THETA){
  B <- THETA[-length(THETA)]
  S2 <- THETA[length(THETA)]
  aBETA <- (1/sigma2)*t(BETA - B)%*%(t(Xno)%*%Xno)
  aS2 <- -v/(2*sigma2)
  alpha <- c(aBETA, aS2)
  return(alpha)
}

aanew <- apply(RegMCMC, 1, falpha) 
dim(aanew)

SIGMABETA <- sigma2*solve(t(Xno)%*%Xno)
SIGMAS2 <- 2*sigma2^2/v
SIGMA <- Matrix::bdiag(SIGMABETA, SIGMAS2)
B <- 6000
ppnew <- rep(1,B)

fdeltaMax <- function(B){
  a <- 3*B[3]
  b <- 2*B[2]
  c <- B[1]
  pibMax_f <- (- b + (b^2-4*a*c)^0.5)/(2*a)
  return(pibMax_f)
} 
ttnewMax <- apply(RegMCMC[,1:3], 1, fdeltaMax)

source("freqacc.R")
fsdMin <- freqacc(tt = ttnewMin , aa = t(aanew), pp=ppnew, V = SIGMA, sw=0)
fsdMin
LimInfMin <- fsdMin$Ebayes-1.96*fsdMin$sdbayes
LimSupMin <- fsdMin$Ebayes+1.96*fsdMin$sdbayes
fsdMin$Ebayes
LimInfMin
LimSupMin


fsdMax <- freqacc(tt = ttnewMax , aa = t(aanew), pp=ppnew, V = SIGMA, sw=0)
fsdMax
LimInfMax <- fsdMax$Ebayes-1.96*fsdMax$sdbayes
LimSupMax <- fsdMax$Ebayes+1.96*fsdMax$sdbayes
fsdMax$Ebayes
LimInfMax
LimSupMax

##### Slope

funPredDerMin <- function(theta){theta[1]+2*theta[2]*min(lnGDPpc)+3*theta[3]*min(lnGDPpc)^2}
ttnewSlopeMin <- apply(RegMCMC[,1:3], 1, funPredDerMin)
fsdSlopeMin <- freqacc(tt = ttnewSlopeMin , aa = t(aanew), pp=ppnew, V = SIGMA, sw=0)
LimInfSlopeMin <- fsdSlopeMin$Ebayes-1.96*fsdSlopeMin$sdbayes
LimSupSlopeMin <- fsdSlopeMin$Ebayes+1.96*fsdSlopeMin$sdbayes
fsdSlopeMin$Ebayes
LimInfSlopeMin
LimSupSlopeMin

funPredDerMax <- function(theta){theta[1]+2*theta[2]*max(lnGDPpc)+3*theta[3]*max(lnGDPpc)^2}
ttnewSlopeMax <- apply(RegMCMC[,1:3], 1, funPredDerMax)
fsdSlopeMax <- freqacc(tt = ttnewSlopeMax, aa = t(aanew), pp=ppnew, V = SIGMA, sw=0)
LimInfSlopeMax <- fsdSlopeMax$Ebayes-1.96*fsdSlopeMax$sdbayes
LimSupSlopeMax <- fsdSlopeMax$Ebayes+1.96*fsdSlopeMax$sdbayes
fsdSlopeMax$Ebayes
LimInfSlopeMax
LimSupSlopeMax

# y2D0 <- 2*0.0203166287/(6*0.0006630751) # This is the value makes second derivative 0 
# y2D0 <- LimInfMin
# y2D0 <- LimSupMin
y2D0 <- 11.18 # This value makes LimSupSlopeMinLim with derivative equal 0
# funPredDer2 <- function(theta){2*theta[1]+6*theta[2]*y2D0}

funPredDerMinLim <- function(theta){theta[1]+2*theta[2]*y2D0+3*theta[3]*y2D0^2}
ttnewSlopeMinLim <- apply(RegMCMC[,1:3], 1, funPredDerMinLim)
fsdSlopeMinLim <- freqacc(tt = ttnewSlopeMinLim , aa = t(aanew), pp=ppnew, V = SIGMA, sw=0)
LimInfSlopeMinLim <- fsdSlopeMinLim$Ebayes-1.96*fsdSlopeMinLim$sdbayes
LimSupSlopeMinLim <- fsdSlopeMinLim$Ebayes+1.96*fsdSlopeMinLim$sdbayes
fsdSlopeMinLim$Ebayes
fsdSlopeMinLim$sdbayes
LimInfSlopeMinLim
LimSupSlopeMinLim


########## Figures #########
require(ggplot2)
le <- length(seq(max(lnGDPpc[IdOrd]), 13, 0.1))
giniNA <- c(gini, rep(NA, le))
df <- data.frame(cbind(giniNA, lnfig, MeanPred))

theta0 <- c(Cte, thetamean[1:3]) 
y1 <- theta0[1]+theta0[2]*fsdMin$Ebayes+theta0[3]*fsdMin$Ebayes^2+theta0[4]*fsdMin$Ebayes^3
y2 <- theta0[1]+theta0[2]*fsdMax$Ebayes+theta0[3]*fsdMax$Ebayes^2+theta0[4]*fsdMax$Ebayes^3

MeanPred <- funPred(theta0)
# colores <- c("Gini" = "black", "Prediction" = "blue", "Roots" = "red", "Extreme log(gdp)" = "green", "C.I. 95%" = "purple")

FigL <- ggplot(df, aes(lnfig)) +                    # basic graphical object
  geom_point(aes(y=giniNA), colour="black") +  # first layer
  geom_line(aes(y=MeanPred), colour="blue", linewidth = 1.2) +  # second layer
  labs(x = "lnGDP", y = "Gini") +
  geom_vline(xintercept=pibMin, colour = "red") +
  geom_vline(xintercept=pibMax, colour = "red") +
  geom_vline(xintercept=min(lnGDPpc), colour = "green") +
  geom_vline(xintercept=max(lnGDPpc), colour = "green") +
  geom_segment(aes(x = LimInfMin, y = y1, xend = LimSupMin, yend = y1), colour = "purple", linewidth = 3) 
# +
#   geom_segment(aes(x = LimInfMax, y = y2, xend = LimSupMax, yend = y2, colour = "credible set"), linewidth = 1.2)

FigL 

df1 <- data.frame(cbind(giniNA, lnfig, MeanPredDer))

# theta0 <- thetamean[1:3]
# y1 <- theta0[1]+2*theta0[2]*fsdSlopeMax$Ebayes+3*theta0[3]*fsdSlopeMax$Ebayes^2
# y11 <- theta0[1]+2*theta0[2]*LimInfSlopeMax+3*theta0[3]*LimInfSlopeMax^2
# y11 <- theta0[1]+2*theta0[2]*LimSupSlopeMax+3*theta0[3]*LimSupSlopeMax^2


FigD <- ggplot(df1, aes(lnfig)) +                    # basic graphical object
  geom_line(aes(y=MeanPredDer), colour="blue", linewidth = 1.2) +  # second layer
  xlab("lnGDP") + ylab("Derivative") +
  geom_vline(xintercept=pibMin, colour = "red") +
  geom_vline(xintercept=pibMax, colour = "red") +
  geom_hline(yintercept=0, colour = "red") +
  geom_vline(xintercept=min(lnGDPpc), colour = "green") +
  geom_vline(xintercept=max(lnGDPpc), colour = "green") +
  geom_segment(aes(x = min(lnGDPpc), y = LimInfSlopeMin, xend = min(lnGDPpc), yend = LimSupSlopeMin), colour = "purple", linewidth = 3) +
  geom_segment(aes(x = max(lnGDPpc), y = LimInfSlopeMax, xend = max(lnGDPpc), yend = LimSupSlopeMax), colour = "purple", linewidth = 3)  +
  geom_segment(aes(x = y2D0, y = LimInfSlopeMinLim, xend = y2D0, yend = LimSupSlopeMinLim), colour = "purple", linewidth = 3) 


FigD

df2 <- data.frame(cbind(giniNA, lnfig, MeanPredDer2))

FigD2 <- ggplot(df2, aes(lnfig)) +                    # basic graphical object
  geom_point(aes(y=giniNA), colour="black") +  # first layer
  geom_line(aes(y=MeanPredDer2), colour="blue", linewidth = 1.2) +  # second layer
  xlab("lnGDP") + ylab("Gini") +
  geom_vline(xintercept=pibMin, colour = "red") +
  geom_vline(xintercept=pibMax, colour = "red") +
  geom_hline(yintercept=0, colour = "red") +
  geom_vline(xintercept=min(lnGDPpc), colour = "green") +
  geom_vline(xintercept=max(lnGDPpc), colour = "green")



exp(pibMin) # around country position 80 by GDPpc
exp(pibMax) # around country position 160 by GDPpc


library(ggpubr)
ggarrange(FigL, FigD, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)


require(dplyr)
SumData <- dataI %>% 
  mutate(GDPpc = exp(lnGDPpc)) %>% 
  group_by(country) %>%  
  summarise(mean_GDPpc = mean(GDPpc), 
            .groups = "drop") 


###### Figures with Lessmann values ######

thetaL1 <- c(0.205, -0.022, 0.001) 
MeanPredDer <- funPredDer(thetaL1)
dfL1 <- data.frame(cbind(lnfig, MeanPredDer))
FigL1 <- ggplot(dfL1, aes(lnfig)) +                    # basic graphical object
  geom_line(aes(y=MeanPredDer), colour="blue", linewidth = 1.2) +  # second layer
  xlab("lnGDP") + ylab("Gini") +
  geom_vline(xintercept=pibMin, colour = "red") +
  geom_vline(xintercept=pibMax, colour = "red") +
  geom_hline(yintercept=0, colour = "red")

FigL1

