#--------------------------------
# MAIN FUNCTIONS: Estimate MIDAS quantile regression
#--------------------------------

#' @export MidasQuantile
MidasQuantile <- function(y,yDate,x = NULL, xDate = NULL, q = 0.01,
                          horizon = 10, nlag = 100, ovlap = FALSE, numInitialsRand = 50000,constrained = TRUE,
                          numInitials = 20, GetSe = TRUE, GetSeSim = 200, Params = NULL, startPars = NULL,
                          MainSolver = "neldermead",SecondSolver = "ucminf",As = FALSE,forecastLength = 0,
                          fitcontrol = list(rep = 5),beta2para = NULL,warn = TRUE, simpleRet = FALSE,multiSol = TRUE){
  #-- set up arguments ----
  if(length(yDate) != length(y))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
  y[is.na(y)] = mean(y,na.rm = TRUE)
  if(is.na(match(MainSolver,c("bobyqa","ucminf","neldermead","solnp","bfgs")))){
    stop("\nMidasQuantile-->error: available solvers are bobyqa, ucminf, bfgs, solnp and neldermead... \n")
  }
  if(!is.null(SecondSolver)){
    if(is.na(match(SecondSolver,c("bobyqa","ucminf","neldermead","solnp","bfgs")))){
      stop("\nMidasQuantile-->error: Available solvers are bobyqa, ucminf, solnp, bfgs and neldermead... \n")
    }
  }
  y = y[1:(length(y) - forecastLength)]
  yDate = yDate[1:(length(yDate) - forecastLength)]
  if(is.null(x)){
    if(simpleRet){
      x = abs(exp(y) - 1)
    } else{
      x = abs(y) # If regressor is null, using absolute returns
    }
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
  if(is.null(beta2para)||!beta2para){
    beta2para = FALSE
    if(As){
      lb = c(-Inf,-Inf,-Inf,1)
      ub = c(Inf, Inf, Inf, 400)
    } else {
      lb = c(-Inf,-Inf,1)
      ub = c(Inf,Inf,400)
    }
  }else{
    if(As){
      lb = c(-Inf,-Inf,-Inf,1, 1)
      ub = c(Inf, Inf, Inf, 400, 400)
    } else {
      lb = c(-Inf,-Inf, 1, 1)
      ub = c(Inf,Inf,400, 400)
    }
  }
  #----- Get the initial guess for the parameters-----
  # It seems like the best way to find the initial beta guess in case of non-ovlapping dataset
  # is to run an overlapping estimation to ultilize the rich of data
  if(is.null(Params)){
  if(!is.null(startPars)){
    betaIni = matrix(rep(startPars,numInitials),nrow = numInitials, byrow = TRUE)
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
  if(is.null(Params)){
  if(is.null(startPars)){# && (horizon == 1 || (horizon > 1 && ovlap) || is.na(betaIni))){
    betaIni = GetIniParams_midas(y = y, X = x,X_neg = x_neg,X_pos = x_pos, q = q, numInitialsRand = numInitialsRand,
                                                          numInitials = numInitials, beta2para = beta2para,As = As)
  }
  #----- Estimate the paramters -----------
  sol = try(.sol(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFun_midas,
             y = y, x = x,x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para,
             lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As,multiSol = multiSol),silent = TRUE)
  if(inherits(sol,"try-error")){
    warnings("\nBoth Solvers failed to converge, resume with other available solvers...\n")
    out = list(estPars = NA, yLowFreq = y, yDate = yDate, betaIni = betaIni,x = x, simpleRet = simpleRet,fitcontrol = fitcontrol,betaIni = betaIni,
               x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para, nlag = nlag, horizon = horizon, GetSe = GetSe, GetSeSim = GetSeSim,
               lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As, conv = 1, forecastLength = forecastLength, constrained = constrained)
  }else{
  for(i in 1:3){
    if(sol$convergence == 0 && (sum(sol$par >= lb) + sum(sol$par <= ub)) == (2 * length(lb))){
      break
    } else{
      if(sol$convergence == 1){
        sol = .sol(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFun_midas,
                   y = y, x = x,x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para,
                   lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As,multiSol = multiSol)
      } else if(constrained){
        sol = .sol(MainSolver = "bobyqa", betaIni = betaIni,fun = objFun_midas,
                   y = y, x = x,x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para,
                   lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As,multiSol = FALSE)
      } else{
       # Do nothing
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
    warnings("\nBoth Solvers failed to converge, resume with other available solvers...\n")
    out = list(estPars = NA, yLowFreq = y, yDate = yDate, betaIni = betaIni,x = x, simpleRet = simpleRet,fitcontrol = fitcontrol,betaIni = betaIni,
               x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para, nlag = nlag, horizon = horizon, GetSe = GetSe, GetSeSim = GetSeSim,
               lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As, conv = 1, forecastLength = forecastLength, constrained = constrained)
  } else{
    condVaR = condVaR_midas(params = estPars,Xr = x,Xr_neg = x_neg,Xr_pos = x_pos,beta2para = beta2para,As = As)
    if(GetSe){
      betaIniSim = matrix(estPars,nrow = 1)
      if(beta2para){
        if(As){
          hypothesis = c(0,0,0,1,1)
        } else {
          hypothesis = c(0,0,1,1)
        }
      } else{
        if(As){
          hypothesis = c(0,0,0,1)
        } else {
          hypothesis = c(0,0,1)
        }
      }
      resid = (y - condVaR)/abs(condVaR)
      parSim = matrix(0,nrow = length(estPars), ncol = GetSeSim)
      for(i in 1:GetSeSim){
        residSim = sample(resid,size = length(resid),replace = TRUE)
        ySim = condVaR + residSim*abs(condVaR)
        if(constrained){
          estSim = .solverSwitch(solver = "bobyqa", pars = betaIniSim, y = ySim, x = x,x_neg = x_neg,
                                 x_pos = x_pos, fun = objFun_midas,q = q, beta2para = beta2para, lb = lb, ub = ub,
                                 control = fitcontrol, As = As)
        }else{
         estSim = .solverSwitch(solver = MainSolver, pars = betaIniSim, y = ySim, x = x,x_neg = x_neg,
                               x_pos = x_pos, fun = objFun_midas,q = q, beta2para = beta2para, lb = lb, ub = ub,
                               control = fitcontrol, As = As)
        }
        parSim[,i] = estSim$par
      }
      pval = rowMeans(abs(parSim - rowMeans(parSim,na.rm = TRUE) + hypothesis) > matrix(rep(abs(estPars),GetSeSim),nrow = length(estPars)),na.rm = TRUE)
    } else{
      pval = rep(NA,length(estPars))
    }
    ViolateRate <- (sum(y < condVaR))/(length(y))
    out = list(estPars = estPars, pval = pval, yLowFreq = y, yDate = yDate, x = x, condVaR = condVaR,simpleRet = simpleRet,As = As,
               ViolateRate = ViolateRate, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver),constrained = constrained,
               lb = lb, ub = ub, horizon = horizon, conv = convergeFlag,forecastLength = forecastLength,ovlap = ovlap,nlag = nlag)
    }
  return(out)
}


#-------------------------------------------------
# Forecasting function
#-------------------------------------------------
#' @export MidasQuantileFor

MidasQuantileFor <- function(object, y, yDate, x = NULL, xDate = NULL){
  #----- Recall the in-sample argument ----
  estPars = object$estPars; beta2para = object$beta2para; As = object$As
  horizon = object$horizon; ovlap = object$ovlap; simpleRet = object$simpleRet; nlag = object$nlag

  #-- set up arguments ----
  if(length(yDate) != length(y))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
  y[is.na(y)] = mean(y,na.rm = TRUE)
  if(is.null(x)){
    if(simpleRet){
      x = abs(exp(y) - 1)
    } else{
      x = abs(y) # If regressor is null, using absolute returns
    }
  } else{
    if(length(x) != length(y))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
  }
  if(is.null(xDate)){
    xDate = yDate
  } else{
    if(length(xDate) != length(yDate))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
  }

   if(!beta2para){
    if(As){
      lb = c(-Inf,-Inf,-Inf,1)
      ub = c(Inf, Inf, Inf, 200)
    } else {
      lb = c(-Inf,-Inf,1)
      ub = c(Inf,Inf,200)
    }
  }else{
    if(As){
      lb = c(-Inf,-Inf,-Inf,1,1)
      ub = c(Inf, Inf, Inf, 200, 200)
    } else {
      lb = c(-Inf,-Inf,1, 1)
      ub = c(Inf,Inf,200, 200)
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

  #------ Compute condVaR and return -----
  condVaR = condVaR_midas(params = estPars,Xr = x,Xr_neg = x_neg,Xr_pos = x_pos,beta2para = beta2para,As = As)
  ViolateRate <- (sum(y < condVaR))/(length(y))
  OSout = list(estPars = estPars, yLowFreq = y, yDate = yDate, condVaR = condVaR,simpleRet = simpleRet,
             ViolateRate = ViolateRate, quantile = object$quantile, beta2para = beta2para,
             horizon = horizon, ovlap = ovlap)
  out = list(InSample = object, OutSample = OSout)
  return(out)
}

#-------------------------------------------------
# Resume function
#-------------------------------------------------
#' @export MidasQuantileResume

MidasQuantileResume <- function(object, MainSolver = "neldermead", SecondSolver = "bobyqa", multiSol = FALSE){
  if(is.na(match(MainSolver,c("bobyqa","ucminf","neldermead","solnp","bfgs")))){
    stop("\nMidasQuantile-->error: available solvers are bobyqa, ucminf and neldermead... \n")
  }
  if(!is.null(SecondSolver)){
    if(is.na(match(SecondSolver,c("bobyqa","ucminf","neldermead","solnp","bfgs")))){
      stop("\nMidasQuantile-->error: Available solvers are bobyqa, ucminf, solnp, bfgs and neldermead... \n")
    }
  }
  #----- Recall the argument -----
  y = object$yLowFreq; yDate = object$yDate;startPars = NULL;
  x = object$x; x_neg = object$x_neg; x_pos = object$x_pos; q = object$q
  beta2para = object$beta2para; lb = object$lb; ub = object$ub; simpleRet = object$simpleRet
  fitcontrol = object$control; warn = object$warn; As= object$As; ovlap = object$ovlap
  nlag = object$nlag; forecastLength = object$forecastLength; horizon = object$horizon
  GetSe = object$GetSe; GetSeSim = object$GetSeSim;constrained = object$constrained;
  if(is.null(startPars)){
    betaIni = GetIniParams_midas(y = y, X = x,X_neg = x_neg,X_pos = x_pos, q = q, numInitialsRand = 50000,
                                 numInitials = 20, beta2para = beta2para,As = As)
  } else{
    betaIni = matrix(rep(startPars,numInitials),nrow = numInitials, byrow = TRUE)
  }
  #----- Estimate the paramters -----------
  sol = try(.sol(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFun_midas,
             y = y, x = x,x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para,
             lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As,multiSol = multiSol),silent = TRUE)
  if(inherits(sol,"try-error")){
    warnings("\nBoth Solvers failed to converge, resume with other available solvers...\n")
    out = list(estPars = NA, yLowFreq = y, yDate = yDate, betaIni = betaIni,x = x, simpleRet = simpleRet,fitcontrol = fitcontrol,
               x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para, nlag = nlag, horizon = horizon, GetSe = GetSe, GetSeSim = GetSeSim,
               lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As, conv = convergeFlag, forecastLength = forecastLength, constrained = constrained)
  }else{
  for(i in 1:3){
    if(sol$convergence == 0 && (sum(sol$par >= lb) + sum(sol$par <= ub)) == (2 * length(lb))){
      break
    } else{
      if(sol$convergence == 1){
        sol = .sol(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFun_midas,
                   y = y, x = x,x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para,
                   lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As,multiSol = multiSol)
      } else if(constrained){
        sol = .sol(MainSolver = "bobyqa", betaIni = betaIni,fun = objFun_midas,
                   y = y, x = x,x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para,
                   lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As,multiSol = FALSE)
      } else{
        break
      }
    }
  }
  estPars = sol$par
  fval = sol$value
  convergeFlag = sol$convergence
  if(convergeFlag == 1){
    warnings("\nBoth Solvers failed to converge, resume with other available solvers...\n")
    out = list(estPars = NA, yLowFreq = y, yDate = yDate, betaIni = betaIni,x = x, simpleRet = simpleRet,fitcontrol = fitcontrol,
               x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para, nlag = nlag, horizon = horizon, GetSe = GetSe, GetSeSim = GetSeSim,
               lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As, conv = convergeFlag, forecastLength = forecastLength, constrained = constrained)
  } else{
    condVaR = condVaR_midas(params = estPars,Xr = x,Xr_neg = x_neg,Xr_pos = x_pos,beta2para = beta2para,As = As)
    if(q < 0.5 && any(condVaR > 0)){
      out = list(estPars = NA, yLowFreq = y, yDate = yDate, betaIni = betaIni,x = x, simpleRet = simpleRet,fitcontrol = fitcontrol,
                 x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para, nlag = nlag, horizon = horizon, GetSe = GetSe, GetSeSim = GetSeSim,
                 lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As, conv = 1, forecastLength = forecastLength, constrained = constrained)
    }else{
      betaIniSim = matrix(estPars,nrow = 1)
      if(beta2para){
        if(As){
          hypothesis = c(0,0,0,1,1)
        } else {
          hypothesis = c(0,0,1,1)
        }
      } else{
        if(As){
          hypothesis = c(0,0,0,1)
        } else {
          hypothesis = c(0,0,1)
        }
      }
      if(GetSe){
        resid = (y - condVaR)/abs(condVaR)
        parSim = matrix(0,nrow = length(estPars), ncol = GetSeSim)
        for(i in 1:GetSeSim){
          residSim = sample(resid,size = length(resid),replace = TRUE)
          ySim = condVaR + residSim*abs(condVaR)
          if(constrained){
            estSim = .solverSwitch(solver = "bobyqa", pars = betaIniSim, y = ySim, x = x,x_neg = x_neg,
                                   x_pos = x_pos, fun = objFun_midas,q = q, beta2para = beta2para, lb = lb, ub = ub,
                                   control = fitcontrol, As = As)
          }else{
            estSim = .solverSwitch(solver = MainSolver, pars = betaIniSim, y = ySim, x = x,x_neg = x_neg,
                                   x_pos = x_pos, fun = objFun_midas,q = q, beta2para = beta2para, lb = lb, ub = ub,
                                   control = fitcontrol, As = As)
          }
          parSim[,i] = estSim$par
        }
        pval = rowMeans(abs(parSim - rowMeans(parSim,na.rm = TRUE) + hypothesis) > matrix(rep(abs(estPars),GetSeSim),nrow = length(estPars)),na.rm = TRUE)
      } else{
        pval = rep(NA,length(estPars))
      }
      ViolateRate <- (sum(y < condVaR))/(length(y))
      out = list(estPars = estPars, pval = pval, yLowFreq = y, yDate = yDate, condVaR = condVaR,simpleRet = simpleRet,As = As,
                 ViolateRate = ViolateRate, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver),constrained = constrained,
                 horizon = horizon, fval = fval, conv = convergeFlag,forecastLength = forecastLength,ovlap = ovlap,nlag = nlag)
    }
  }
  }
  return(out)
}
