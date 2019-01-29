#--------------------------------
# MAIN FUNCTIONS: Estimate MIDAS quantile regression
#--------------------------------
#' @importFrom fExtremes gpdFit
#' @export VarEs_evt
VarEs_evt <- function(y,yDate,x = NULL, xDate = NULL, q = 0.01,qThreshold = 0.075,multiSol = TRUE,Params = NULL,
                      horizon = 10, nlag = 100, ovlap = FALSE, numInitialsRand = 50000, constrained = FALSE,
                      numInitials = 20, GetSe = FALSE, GetSeSim = 200, startPars = NULL, forecastLength = 0,
                      MainSolver = "bobyqa",SecondSolver = "ucminf",quantEst = "midas", As = FALSE, empQuant = NULL,
                      fitcontrol = list(rep = 5),beta2para = FALSE,warn = TRUE, simpleRet = FALSE,Uni = TRUE){
  #-- set up arguments ----
  if(length(yDate) != length(y))  stop("\nMidasQuantile-->error: Length of y and yDate should be the same\n")
  if(is.na(match(quantEst,c("midas","cav")))){
    stop("\nError: the quantile regression should be either midas or cav..\n")
  }
  if(is.na(match(MainSolver,c("bobyqa","ucminf","neldermead","solnp","bfgs")))){
    stop("\nMidasQuantile-->error: available solvers are bobyqa, ucminf and neldermead... \n")
  }
  if(!is.null(SecondSolver)){
    if(is.na(match(SecondSolver,c("bobyqa","ucminf","neldermead","solnp","bfgs")))){
      stop("\nMidasQuantile-->error: Available solvers are bobyqa, ucminf and neldermead... \n")
    }
  }
   #------ Get the threshold for the extreme value theory-----
    VaRthreshold <- switch (quantEst,
                            cav = CAViaR(y = y, yDate = yDate, x = x, q = qThreshold, horizon = horizon, ovlap = ovlap, numInitialsRand = numInitialsRand,
                                         numInitials = numInitials, empQuant = empQuant, GetSe = GetSe, startPars = startPars,multiSol = multiSol,
                                         MainSolver = MainSolver,SecondSolver = SecondSolver, As = As, fitcontrol = list(rep = 5),Params = Params,
                                         warn = warn, simpleRet = simpleRet, Uni = Uni, forecastLength = forecastLength, constrained = constrained), 
                            midas = MidasQuantile(y = y, yDate = yDate, x = x, xDate = xDate, q = qThreshold, horizon = horizon, nlag = nlag, ovlap = ovlap,
                                                  numInitialsRand = numInitialsRand, numInitials = numInitials, GetSe = GetSe, GetSeSim = GetSeSim, Params = Params,
                                                  startPars = startPars, MainSolver = MainSolver, SecondSolver = SecondSolver, As = As, fitcontrol = fitcontrol,
                                                  beta2para = beta2para, warn = warn, constrained = constrained, simpleRet = simpleRet, forecastLength = forecastLength,multiSol = multiSol)
    )
    
    #------ Estimate the General Pareto Distribution parameters----
    if(VaRthreshold$conv == 1){
      warning("\nThe quantile regression for threshold level is not converged, refit the model with other solver...\n")
      out = list(estPars = NA,thresholdEst = NA, yLowFreq = y, yDate = yDate, condVaR = NA,condES = NA,
                 GPDest = NA, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver), simpleRet = simpleRet,
                 As = As, Uni = Uni, forecastLength = forecastLength,empQuant = empQuant, quantEst = quantEst, conv = 1)
      
    } else{
      condThres = VaRthreshold$condVaR
      y = VaRthreshold$yLowFreq
      yDate = VaRthreshold$yDate
      StdExceed = y/condThres - 1
      GPDfit = try(fExtremes::gpdFit(x = StdExceed,u = 0),silent = TRUE)
      if(inherits(GPDfit,"try-error")){
        warning("\n Unrealistic estimate of GPD paramters or positive condVaR, reestimate with another threshold level..\n")
        out = list(estPars = NA,thresholdEst = VaRthreshold, yLowFreq = y, yDate = yDate, condVaR = NA,condES = NA,
                   GPDest = NA, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver), simpleRet = simpleRet,
                   As = As, Uni = Uni,ovlap = ovlap, forecastLength = forecastLength,empQuant = empQuant, quantEst = quantEst, conv = 1, horizon = horizon)
      }else{
        GPDpars = GPDfit@fit$par.ests
        gamma = GPDpars[1]; beta = GPDpars[2]
        Threshold_ViolateRate <- VaRthreshold$ViolateRate
        StdQuantile = ((q*(1/Threshold_ViolateRate))^(-gamma)-1)*(beta/gamma)
        condVaR = condThres * (1 + StdQuantile)
        ExpectedMean_StdQuant = StdQuantile * (1/(1-gamma) + beta/((1-gamma)*StdQuantile))
        condES = condThres * (1 + ExpectedMean_StdQuant)
        if(gamma >= 1 || any(condVaR >= 0)){
          warning("\n Unrealistic estimate of GPD paramters or positive condVaR, reestimate with another threshold level..\n")
          out = list(estPars = c(VaRthreshold$estPars,unname(GPDfit@fit$par.ests)),thresholdEst = VaRthreshold, yLowFreq = y, yDate = yDate, condVaR = condVaR,condES = condES,
                   GPDest = GPDfit, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver), simpleRet = simpleRet,
                   As = As, Uni = Uni,ovlap = ovlap, forecastLength = forecastLength,empQuant = empQuant, quantEst = quantEst, conv = 1, horizon = horizon)
        } else{
        ViolateRate <- (sum(y < condVaR))/(length(y))
        out = list(estPars = c(VaRthreshold$estPars,unname(GPDfit@fit$par.ests)),thresholdEst = VaRthreshold, yLowFreq = y, yDate = yDate, condVaR = condVaR,condES = condES,
                   GPDest = GPDfit, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver), simpleRet = simpleRet,
                   As = As, Uni = Uni,ovlap = ovlap, forecastLength = forecastLength,empQuant = empQuant, quantEst = quantEst, horizon = horizon, ViolateRate = ViolateRate, conv = 0)
        }
      }
    }
  return(out)
}

#---------------------------------
# Forecasting function
#----------------------------------
#' @export VarEs_evt_for

VarEs_evt_for <- function(object, y, yDate, x = NULL, xDate = NULL){
  #--------- Recall InSample arguments-----
  quantEst <- object$quantEst;  estPars = object$estPars; ISthresholdEst = object$thresholdEst; q = object$quantile
  beta2para = object$beta2para; simpleRet = object$simpleRet; As = object$As; Uni = object$Uni;
  quantEst = object$quantEst; horizon = object$horizon;ovlap = object$ovlap
  #----- Get the condVaR and condEs -----
  n = length(estPars)
  gamma = estPars[n-1]; beta = estPars[n]
  VaRthreshold <- switch (quantEst,
                          cav = CAViaRfor(object = ISthresholdEst, y = y, yDate = yDate, x = x, xDate = xDate),
                          midas = MidasQuantileFor(object = ISthresholdEst, y = y, yDate = yDate, x = x, xDate = xDate))
  VaRthreshold = VaRthreshold$OutSample
  condThres = VaRthreshold$condVaR
  y = VaRthreshold$yLowFreq
  yDate = VaRthreshold$yDate
  StdExceed = y/condThres - 1
  IS_Threshold_ViolateRate <- ISthresholdEst$ViolateRate
  StdQuantile = ((q*(1/IS_Threshold_ViolateRate))^(-gamma)-1)*(beta/gamma)
  condVaR = condThres * (1 + StdQuantile)
  ExpectedMean_StdQuant = StdQuantile * (1/(1-gamma) + beta/((1-gamma)*StdQuantile))
  condES = condThres * (1 + ExpectedMean_StdQuant)
  ViolateRate <- (sum(y < condVaR))/(length(y))
  OSout = list(yLowFreq = y, yDate = yDate, condVaR = condVaR,condES = condES,VaRthreshold = VaRthreshold,
             quantile = q, beta2para = beta2para, simpleRet = simpleRet,
             As = As, Uni = Uni,quantEst = quantEst, ViolateRate = ViolateRate)
  out <- list(InSample = object, OutSample = OSout)
  return(out)
}

