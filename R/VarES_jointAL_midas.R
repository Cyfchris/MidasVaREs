#--------------------------------
# MAIN FUNCTIONS: Estimate MIDAS quantile regression and ES using joint estimation with AL distribution
#--------------------------------
#' @importFrom forecast auto.arima Arima
#' @importFrom lmtest coeftest

#' @export VarEs_jointAL_midas
VarEs_jointAL_midas <- function(y,yDate,x = NULL, xDate = NULL, q = 0.01, armaOrder = c(1,0,0), horizon = 10, 
                                nlag = 100, ovlap = FALSE, numInitialsRand = 50000, numInitials = 20, 
                                GetSe = TRUE, GetSeSim = 200, startPars = NULL,constrained = TRUE,Params = NULL,
                                MainSolver = "bobyqa",SecondSolver = "ucminf",As = FALSE, fitcontrol = list(rep = 5),
                                beta2para = FALSE,warn = TRUE, simpleRet = FALSE,forecastLength = 0, multiSol = TRUE){
  #-- set up arguments ----
  if(length(yDate) != length(y))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
  y[is.na(y)] = mean(y,na.rm = TRUE)
  if(is.na(match(MainSolver,c("bobyqa","ucminf","neldermead","solnp","bfgs")))){
    stop("\nMidasQuantile-->error: available solvers are bobyqa, ucminf and neldermead... \n")
  }
  if(!is.null(SecondSolver)){
    if(is.na(match(SecondSolver,c("bobyqa","ucminf","neldermead","solnp","bfgs")))){
      stop("\nMidasQuantile-->error: Available solvers are bobyqa, ucminf and neldermead... \n")
    }
  }
  y = y[1:(length(y) - forecastLength)]
  yDate = yDate[1:(length(yDate) - forecastLength)]
  if(is.null(x)){
    x = abs(y) # If regressor is null, using absolute returns
  } else{
    if(length(x) != length(y))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
    x = x[1:(length(x) - forecastLength)]
  }
  if(is.null(xDate)){
    xDate = yDate
  } else{
    if(length(xDate) != length(yDate))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
    xDate = xDate[1:(length(xDate) - forecastLength)]
  }
  if(!beta2para){
    if(As){
      lb = c(-Inf,-Inf,-Inf,1,-Inf)
      ub = c(Inf, Inf, Inf, 400, Inf)
    } else {
      lb = c(-Inf,-Inf,1,-Inf)
      ub = c(Inf,Inf,400, Inf)
    }
  }else{
    if(As){
      lb = c(-Inf,-Inf,-Inf,1, 1,-Inf)
      ub = c(Inf, Inf, Inf, 400, 400, Inf)
    } else {
      lb = c(-Inf,-Inf,1, 1,-Inf)
      ub = c(Inf,Inf,400, 400, Inf)
    }
  }
  if(is.null(Params)){
    if(is.null(startPars)){
    UniQuantEst = try(MidasQuantile(y = y, yDate = yDate, x = x, xDate = xDate, q = q,constrained = constrained,
                                horizon = horizon, nlag = nlag, ovlap = ovlap, numInitialsRand = numInitialsRand,
                                numInitials = numInitials, GetSe = FALSE, GetSeSim = NULL, Params = NULL, multiSol = multiSol,
                                startPars = NULL, MainSolver = MainSolver, SecondSolver = SecondSolver, As = As,
                                fitcontrol = fitcontrol, beta2para = beta2para, warn = FALSE, simpleRet = simpleRet),silent = TRUE)
    if(inherits(UniQuantEst,"try-error")){
      warning("\nMidasQuantile -->error: The univariate quantile does not converge, try other solvers.\n")
      out = list(estPars = NA, pval = NA, y = y, yDate = yDate, 
                 condVaR = NA, condES = NA, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver),
                 fval = NA, conv = 1,multiSol = multiSol,forecastLength = forecastLength)
    } else{
      if(UniQuantEst$conv == 1){
      for(i in 1:3){
      UniQuantEst <- MidasQuantileResume(UniQuantEst)
        if(UniQuantEst$conv == 0){
          break
        }
      }
      if(UniQuantEst$conv == 1){
        warning("\nMidasQuantile -->error: The univariate quantile does not converge, try other solvers.\n")
        out = list(estPars = NA, pval = NA, y = y, yDate = yDate, 
                   condVaR = NA, condES = NA, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver),
                   fval = NA, conv = 1,multiSol = multiSol,forecastLength = forecastLength)
      }
      }
    }
  } 
  }
  #------ Get the mixed data to start estimation -----
  dataEst <- MixedFreqQuant(DataY = y,DataYdate = yDate,DataX = x,DataXdate = xDate,
                            xlag = nlag,period = horizon,ovlap = ovlap, simpleRet = simpleRet)
  dataHigh <- MixedFreqQuant(DataY = y,DataYdate = yDate,DataX = y,DataXdate = yDate,
                             xlag = nlag,period = horizon,ovlap = ovlap, simpleRet = simpleRet)
  y = dataEst$EstY;
  yDate <- dataEst$EstYdate
  x = dataEst$EstX
  xDate = dataEst$EstXdate
  yHigh = dataHigh$EstX
  x_neg = x
  x_pos = x
  x_neg[yHigh > 0] = 0
  x_pos[yHigh <= 0] = 0
  #------- Fit the conditional mean equation-------------
  if(is.null(armaOrder)){
    meanFit <- forecast::auto.arima(y, max.d = 0, max.D = 0)
    if(length(meanFit$coef) == 0){
      meanCoef = 0
      condMean = rep(0,length(y))
    } else{
      meanCoef <- lmtest::coeftest(meanFit)
      condMean = as.numeric(meanFit$fitted)
    }
  } else{
    meanFit = forecast::Arima(y,order = armaOrder)
    meanEst <- lmtest::coeftest(meanFit)
    meanCoef = unname(meanEst[,1])
    condMean = as.numeric(meanFit$fitted)
  }
  
  #----- Get the initial guess for the parameters-----
  if(is.null(Params)){
    if(is.null(startPars)){
    betaIni = GetIniParamsAL_midas(y = y, condMean = condMean, QuantEst = UniQuantEst$estPars, X = x, 
                                 X_neg = x_neg,X_pos = x_pos, q = q, numInitialsRand = numInitialsRand,
                                 numInitials = numInitials, beta2para = beta2para,As = As)
  }  else{
    betaIni = matrix(rep(startPars,numInitials),nrow = numInitials, byrow = TRUE)
  }
  #----- Estimate the paramters -----------
  sol = try(.sol(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFunAL_midas,
             condMean = condMean, y = y, x = x,x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para, 
             lb = lb, ub = ub, control = fitcontrol,warn = warn,multiSol = multiSol,As=As),silent = TRUE)
  if(inherits(sol,"try-error")){
    estPars = NA
    convergeFlag = 1
  } else{
  for(i in 1:3){
    if(sol$convergence == 0 && (sum(sol$par >= lb) + sum(sol$par <= ub)) == (2 * length(lb))){
      break
    } else{
      if(sol$convergence == 1){
        sol = .sol(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFunAL_midas,
                   condMean = condMean, y = y, x = x,x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para, 
                   lb = lb, ub = ub, control = fitcontrol,warn = warn,multiSol = multiSol,As=As)
      } else if(constrained){
        sol = .sol(MainSolver = "bobyqa",betaIni = betaIni,fun = objFunAL_midas,
                   condMean = condMean, y = y, x = x,x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para, 
                   lb = lb, ub = ub, control = fitcontrol,warn = warn,multiSol = FALSE,As=As)
      } else{

      }
    }
  }
    estPars = sol$par
    convergeFlag = sol$convergence
  } 
  } else{
    estPars = Params
    convergeFlag = 0
  }
  
  if(convergeFlag == 1){
    warnings("\nBoth Solvers failed to converge, try with other available solvers...\n")
    out = list(estPars = NA, pval = NA, y = y, yDate = yDate, 
               condVaR = NA, condES = NA, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver),
               conv = 1,multiSol = multiSol,forecastLength = forecastLength)
  } else if(convergeFlag == 0 && constrained && (sum(estPars >= lb)!=length(lb)||sum(estPars <= ub) != length(ub))){
    out = list(estPars = NA, pval = NA, y = y, yDate = yDate, 
               condVaR = NA, condES = NA, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver),
               conv = 1,multiSol = multiSol,forecastLength = forecastLength)
  } else{
    VaRES = condVaRES_midas(params = estPars, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos,beta2para = beta2para,As = As)
    betaIniSim = matrix(estPars,nrow = 1)
    condVaR = VaRES$VaR
    condES = VaRES$ES
    if((q < 0.3) && (any(condVaR > 0)||any(condES > 0))){
      constrainedLot = which(condVaR > 0)
      condVaR[constrainedLot] = -0.0001
      condES[constrainedLot] = condVaR[constrainedLot] * (1 + exp(tail(estPars,1)))
    }
    if(GetSe){
       if(beta2para){
        if(As){
          hypothesis = c(0,0,0,1,1,0)
        } else {
          hypothesis = c(0,0,1,1,0)
        }
      } else{
        if(As){
          hypothesis = c(0,0,0,1,0)
        } else {
          hypothesis = c(0,0,1,0)
        }
      }
      resid = (y - condVaR)/abs(condVaR)
      parSim = matrix(0,nrow = length(estPars), ncol = GetSeSim)
      for(i in 1:GetSeSim){
        residSim = sample(resid,size = length(resid),replace = TRUE)
        ySim = condVaR + residSim*abs(condVaR)
        ESsimMean = mean(residSim[residSim<0])
        Exceed = which(ySim < condVaR)
        ySim[Exceed] = condVaR[Exceed] + (residSim[Exceed]/ESsimMean)*(condES[Exceed]-condVaR[Exceed]);
        ySim_filter = forecast::Arima(ySim,order = meanFit$arma[1:3])
        condMeanSim <- as.numeric(ySim_filter$fitted)
       if(constrained){
          estSim = .solverSwitch(solver = "bobyqa", pars = betaIniSim, y = ySim,condMean = condMeanSim, x = x,x_neg = x_neg,
                                 x_pos = x_pos, fun = objFunAL_midas,q = q, beta2para = beta2para, lb = lb, ub = ub,
                                 control = fitcontrol, As = As)
        }else{
          estSim = .solverSwitch(solver = MainSolver, pars = betaIniSim, y = ySim,condMean = condMeanSim, x = x,x_neg = x_neg,
                                 x_pos = x_pos, fun = objFunAL_midas,q = q, beta2para = beta2para, lb = lb, ub = ub,
                                 control = fitcontrol, As = As)
        }
        parSim[,i] = estSim$par
      }
      pval = rowMeans(abs(parSim - rowMeans(parSim,na.rm = TRUE) + hypothesis) > matrix(rep(abs(estPars),GetSeSim),nrow = length(estPars)),na.rm = TRUE)
    } else{
      pval = rep(NA,length(estPars))
    }
    ViolateRate <- (sum(y < condVaR))/(length(y))
    out = list(estPars = estPars, pval = pval, yLowFreq = y, yDate = yDate, condVaR = condVaR, condES = condES, simpleRet = simpleRet,
               ViolateRate = ViolateRate, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver), As = As,
               horizon = horizon, conv = convergeFlag,forecastLength = forecastLength,ovlap = ovlap,nlag = nlag, meanFit = meanFit)
    }
  return(out)
}

#------------------------------------
# Forecasting function
#------------------------------------
#' @export VarEs_jointAL_midas_for

VarEs_jointAL_midas_for <- function(object, y, yDate, x = NULL, xDate = NULL){
  if(length(yDate) != length(y))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
  y[is.na(y)] = mean(y,na.rm = TRUE)
  if(is.null(x)){
    x = abs(y) # If regressor is null, using absolute returns
  } else{
    if(length(x) != length(y))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
  }
  if(is.null(xDate)){
    xDate = yDate
  } else{
    if(length(xDate) != length(yDate))  stop("\nMidasQuantile-->error: Length of yDate and xDate should be the same\n")
  }
  #----- Recall the in-sample argument ----
  estPars = object$estPars; simpleRet = object$simpleRet;q = object$quantile
  beta2para = object$beta2para; horizon = object$horizon
  ovlap = object$ovlap;nlag = object$nlag;As = object$As
  
  #------ Get the mixed data to start estimation -----
  dataEst <- MixedFreqQuant(DataY = y,DataYdate = yDate,DataX = x,DataXdate = xDate,
                            xlag = nlag,period = horizon,ovlap = ovlap, simpleRet = simpleRet)
  dataHigh <- MixedFreqQuant(DataY = y,DataYdate = yDate,DataX = y,DataXdate = yDate,
                             xlag = nlag,period = horizon,ovlap = ovlap, simpleRet = simpleRet)
  y = dataEst$EstY;
  yDate <- dataEst$EstYdate
  x = dataEst$EstX
  xDate = dataEst$EstXdate
  yHigh = dataHigh$EstX
  x_neg = x
  x_pos = x
  x_neg[yHigh > 0] = 0
  x_pos[yHigh <= 0] = 0
  VaRES = condVaRES_midas(params = estPars, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos,beta2para = beta2para,As = As)
  condVaR = VaRES$VaR
  condES = VaRES$ES
  ViolateRate <- (sum(y < condVaR))/(length(y))
  OSout = list(estPars = estPars, yLowFreq = y, yDate = yDate, condVaR = condVaR, condES = condES, simpleRet = simpleRet,
             ViolateRate = ViolateRate, quantile = q, beta2para = beta2para, horizon = horizon, ovlap = ovlap,nlag = nlag)
  out <- list(InSample = object, OutSample = OSout)
  return(out)
}
