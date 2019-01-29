#---
# Loss Function given VaR and ES estimates
#---
#' @importFrom stats coefficients dnorm lm median pchisq pnorm qnorm
#' @importFrom utils head tail
#' @importFrom MASS ginv

#' @export
Loss = function(y, VaR, ES, alpha, choice){
  # Check input
  if(length(y) != length(VaR)){
    stop('Length of y and VaR must be the same')
  }
  if(length(y) != length(ES)){
    stop('If ES is provided, it should be at the same length as VaR')
  }
  if(sum(is.infinite(VaR)) > 0){ 
    stop('VaR can not be Inf - Recheck the VaR estimates')
  }
  if(sum(is.infinite(ES)) > 0){ 
    stop('ES can not be Inf - Recheck the ES estimates')
  }
  if(alpha < 0|alpha > 1){
    stop('Quantile level must be between 0 and 1')
  }
  if(!is.element(choice,c(1,2,3))){
    stop('Choice of loss function need to be either 1, 2 or 3')
  }
  if(choice == 1){
    loss = (y - VaR)*(alpha - (y <= VaR))
    score = mean(loss)
  }else if(choice == 2){
    loss1 = ((y <= VaR) - alpha)*VaR - (y <= VaR)*y
    loss2 = (exp(ES)/(1 + exp(ES)))*(ES - VaR + (y <= VaR)*((VaR-y)/alpha))
    loss3 = log(2/(1 + exp(ES)))
    loss = loss1 + loss2 + loss3
    score = mean(loss1 + loss2 + loss3)
  }else{
    loss = -log((alpha - 1)/ES) - (((y - VaR)*(alpha - (y <= VaR)))/(alpha*ES)) + y/ES;
    score = mean(loss);
  }
  return(list(loss = loss, score = score))
}

KupiecTest = function(Hit,alpha){
  N = length(Hit)
  x = sum(Hit)
  rate = x/N
  test = -2 * log(((1 - alpha)^(N - x) * alpha^x)/((1 - rate)^(N - x) * rate^x))
  if (is.nan(test)) 
    test = -2 * ((N - x) * log(1 - alpha) + x * log(alpha) - 
                   (N - x) * log(1 - rate) - x * log(rate))
  pvalue = 1 - stats::pchisq(test, df = 1)
  LRpof = c(test, pvalue)
  names(LRpof) = c("Test", "Pvalue")
  return(LRpof)  
}

DQtest_VaR = function(y,VaR,alpha,lags = 4){
  cT = length(y)
  vHit = numeric(cT)
  vHit[y <= VaR] = 1 - alpha
  vHit[y > VaR] = -alpha
  vConstant = rep(1, (cT - lags))
  vHIT = vHit[(lags + 1):cT]
  vVaRforecast = VaR[(lags + 1):cT]
  mZ = matrix(0, cT - lags, lags)
  vY2_lag = y[lags:(cT - 1)]^2
  for (st in 1:lags) {
    mZ[, st] = vHit[st:(cT - (lags + 1L - st))]
  }
  mX = cbind(vConstant, vVaRforecast, mZ)
  dDQstatOut = (t(vHIT) %*% mX %*% MASS::ginv(t(mX) %*% mX) %*% 
                  t(mX) %*% (vHIT))/(alpha * (1 - alpha))
  dDQpvalueOut = 1 - stats::pchisq(dDQstatOut, ncol(mX))
  return(dDQpvalueOut)
}

#' @export
var_backtest = function(realized,VaR,alpha,lags = 5){
  Hit = numeric(length(realized))
  Hit[which(realized <= VaR)] = 1L
  N = length(Hit)
  x = sum(Hit)
  rate = x/N
  AE = rate/alpha
  Kupiec = KupiecTest(Hit,alpha)
  DQVaR = DQtest_VaR(realized,VaR,alpha,lags)
  return(list(AE = AE,Kupiec = Kupiec,DQVaR = DQVaR))
}
# Define ES backtest function - McNeil 2000 + Patton 2017 -----
#' @export
es_McNeil = function(realized,VaR,ES,B=1000){
  # Extract standardized exceedances
  stdExceed = ((realized - ES)/VaR)[realized <= VaR]
  # Build the modified stdExceed that have zero mean to build the bootstrap population
  # Refer to Efron and Tibshirani (1993)
  newstdExceed = stdExceed - mean(stdExceed) # The newstdExceed has the mean = 1
  f = function(x)  mean(x) / (stats::sd(x) / sqrt(length(x)))
  t0 = f(stdExceed)
  t = c()
  for (i in 1:1000){ 
    newsample <- sample(newstdExceed, length(newstdExceed), replace=T)
    t <- c(t,f(newsample))
  }
  return(mean(abs(t) > abs(t0)))
}
#' @importFrom car linearHypothesis
es_Patton = function(y,VaR,ES,alpha,lags = 4){
  cT = length(y)
  eHit = numeric(cT)
  eHit[y <= VaR] = (1/alpha) * (y/ES)[y <= VaR] - 1
  eHit[y > VaR] = -1
  eConstant = rep(1, (cT - lags))
  eHIT = eHit[(lags + 1):cT]
  eESforecast = ES[(lags + 1):cT]
  mZ = matrix(0, cT - lags, lags)
  for (st in 1:lags) {
    mZ[, st] = eHit[st:(cT - (lags + 1L - st))]
  }
  reg = lm(eHIT ~ 0 + eConstant + eESforecast + mZ)
  Hnull = names(coefficients(reg))
  for(i in 1:length(Hnull)) Hnull[i] = paste(Hnull[i],"=0",sep = "")
  fTest = car::linearHypothesis(reg,Hnull)
  dDQpvalueOut = fTest$`Pr(>F)`[2]
  return(dDQpvalueOut)
}

# Define the ES backtest of Acerbi (2014) with regards to the Z2 test statistic
#' @export
es_Acerbi = function(y,VaR,ES,alpha,B = 10000){
  z2 = -((1/(length(y)*alpha)) * sum((y*(y <= VaR))/ES)) + 1
  stdResid <- (y - VaR)/abs(VaR)
  meanExceed = mean(stdResid[stdResid<0])
  z2_booted = vector()
  for(i in 1:B){
    stdResid_booted <- sample(stdResid,replace = T)
    y_booted = stdResid_booted*abs(VaR) + VaR
    exceeded_lot = which(y_booted < VaR)
    y_booted[exceeded_lot] = VaR[exceeded_lot] + stdResid_booted[exceeded_lot]*(ES[exceeded_lot] - VaR[exceeded_lot]) / meanExceed
    z2_booted[i] = -((1/(length(y_booted)*alpha)) * sum((y_booted*(y_booted < VaR))/ES)) + 1
  }
  return(c(z2,mean(z2_booted > z2)))
}

#' @export
es_backtest= function(realized,VaR,ES,alpha,lags = 4, B = 10000){
  UnC_backtest = es_McNeil(realized,VaR,ES,B)
  #Patton_backtest = es_Patton(realized,VaR,ES,alpha,lags)
  Acerbi_backtest = es_Acerbi(realized,VaR,ES,alpha,B)
  return(list(UnC = UnC_backtest, Acerbi = Acerbi_backtest))#, Patton = Patton_backtest))
}


