#--------------------------------
# MAIN FUNCTIONS: Estimate MIDAS quantile regression
#--------------------------------
#' @importFrom forecast auto.arima Arima
#' @importFrom lmtest coeftest

#' @export VarEs_jointAL_cav
  VarEs_jointAL_cav <- function(y,yDate,x = NULL, xDate = NULL, q = 0.01, armaOrder = c(1,0,0), horizon = 10, 
                                ovlap = FALSE, numInitialsRand = 50000, numInitials = 20, GetSe = TRUE, 
                                GetSeSim = 200, startPars = NULL,forecastLength = 0,multiSol = TRUE,Params = NULL,
                                MainSolver = "neldermead",SecondSolver = "ucminf",As= FALSE,Uni = TRUE,
                                empQuant = NULL, fitcontrol = list(rep = 5),warn = TRUE, simpleRet = FALSE, constrained = FALSE){
    if(is.na(match(MainSolver,c("bobyqa","ucminf","neldermead","solnp","bfgs")))){
      stop("\nMidasQuantile-->error: available solvers are bobyqa, ucminf and neldermead... \n")
    }
    if(!is.null(SecondSolver)){
      if(is.na(match(SecondSolver,c("bobyqa","ucminf","neldermead","solnp","bfgs")))){
        stop("\nMidasQuantile-->error: Available solvers are bobyqa, ucminf and neldermead... \n")
      }
    }
    #--- Form the single period setting for the return series
    y = y[1:(length(y) - forecastLength)]
    yDate = yDate[1:(length(y) - forecastLength)]
    
    nobs <- length(y)
    nobsShort = nobs-horizon+1
    yLowFreq = matrix(NA,ncol = 1, nrow = nobsShort)
    yDateLowFreq = matrix(NA,ncol = 1, nrow = nobsShort)
    for(t in 1:(nobs - horizon + 1)){
      yLowFreq[t,1] = sum(y[t:(t+horizon-1)]);
      yDateLowFreq[t,1] = yDate[t];
    }
    if(!ovlap){
      yLowFreq = yLowFreq[seq(1,nobsShort,horizon),1]
      yDateLowFreq = yDateLowFreq[seq(1,nobsShort,horizon),1]
    }
    if(simpleRet){
      y = exp(yLowFreq) - 1
    } else{
      y = yLowFreq
    }
    yDate = as.Date.numeric(yDateLowFreq,origin = "1970-01-01")
    if(is.null(x)){
      x = abs(y)
    }else {
      if(length(x) != length(y)){
        stop("\nCAviaR-->error: The explanatory variable should have the same length as the objective variable... \n")
      }
    }
    if(is.null(xDate)) xDate = yDate
    #-----Check the solvers input------
    # The CAViaR model is sensitive to the choice of the empirical quantile to start the dynamics. Here, we use the empirical
    # quantile of the first 10% of the data sample to start the quantile dynamics.
    if(is.null(empQuant)){
      if(horizon == 1){
        empQuant = unname(quantile(y[1:300],q,type = 5))
      } else if(horizon > 1 && horizon <= 5){
        empQuant = unname(quantile(y[1:100],q,type = 5))
      } else{
        empQuant = unname(quantile(y[1:50],q,type = 5))
      }
    }
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
    # Set the bounds for the parameters. The autoregressive paramter should be between 0 and 1?
    tol = 1e-10
    if(is.null(Params)){
    if(is.null(startPars)){
      UniQuantEst = CAViaR(y = y,yDate = yDate,x = x,xDate = xDate,q = q,horizon = 1,ovlap = FALSE,numInitialsRand = numInitialsRand,
                           numInitials = numInitials,empQuant = empQuant, GetSe = FALSE, startPars = NULL,MainSolver = MainSolver,multiSol = multiSol,
                           SecondSolver = SecondSolver,As = As,fitcontrol = fitcontrol,warn = FALSE,simpleRet = FALSE,Uni = Uni, constrained = constrained)
      if(UniQuantEst$conv == 1){
        UniQuantEst <- CAViaRresume(UniQuantEst)
        if(UniQuantEst$conv == 1){
        warning("\nCAViaR -->error: The univariate quantile does not converge, try other solvers.\n")
          out = list(estPars = NA,simpleRet = simpleRet,horizon = horizon,
                     control = fitcontrol,y = y,x = x,As = As,q = q,empQuant = empQuant,ovlap = ovlap,
                     Uni = Uni,condMean = NA, meanFit = NA, conv = 1, multiSol = multiSol)
        }
      }
    }
    #----- Get the initial guess for the parameters-----
    if(is.null(startPars)){
      betaIni = GetIniParamsAL_cav(y = y, condMean = condMean, QuantEst = UniQuantEst$estPars, X = x, 
                                 q = q, numInitialsRand = numInitialsRand, As = As, empQuant = empQuant,
                                 Uni = Uni)
    } else{
      betaIni = matrix(rep(startPars,numInitials),nrow = numInitials, byrow = TRUE)
    }
    #----- Estimate the paramters -----------
    if(!Uni){
      if(As){
        lb = c(-Inf,-Inf,-Inf,-Inf,tol,-Inf)
        ub = c(Inf, Inf, Inf, Inf, 1-tol, Inf)
      } else{
        lb = c(-Inf, -Inf, -Inf, tol, -Inf)
        ub = c(Inf, Inf, Inf, 1-tol, Inf)
      }
    } else{
      if(As){
        lb = c(-Inf,-Inf,-Inf,tol, -Inf)
        ub = c(Inf, Inf, Inf, 1-tol, Inf)
      } else{
        lb = c(-Inf, -Inf, tol, -Inf)
        ub = c(Inf, Inf, 1-tol, Inf)
      }
    }
    sol = .sol_cav(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFunAL_cav,
                   control = fitcontrol,y = y,x = x,As = As,q = q,empQuant = empQuant,
                   Uni = Uni,condMean = condMean, lb = lb, ub = ub, multiSol = multiSol)
      estPars = sol$par
      convergeFlag = sol$convergence
    } else{
      estPars = Params
      convergeFlag = 0
    }
    if(convergeFlag == 1){
      warnings("\nBoth Solvers failed to converge, try with other available solvers...\n")
      out = list(estPars = estPars, betaIni = betaIni,fun = objFunAL_cav,simpleRet = simpleRet,horizon = horizon,
                 control = fitcontrol,y = y,x = x,As = As,q = q,empQuant = empQuant,ovlap = ovlap,
                 Uni = Uni,condMean = condMean, meanFit = meanFit, conv = 1, multiSol = multiSol)
    } else{
      VaRES = condVaRES_cav(params = estPars, yr = y, Xr = x, empQuant = empQuant, As = As, Uni = Uni)
      betaIniSim = matrix(estPars,nrow = 1)
      condVaR = VaRES$VaR
      condES = VaRES$ES
      if((q < 0.25) && (any(condVaR > 0)||any(condES > 0))){
        # In case we have condVaR for q < 0.5, using this method, condES will be larger than 0, which does not make sence.
        # Hence, in this case, we constrained VaR to be slightly less than 0
        constrainedLot = which(condVaR > 0)
        condVaR[constrainedLot] = -0.0001
        condES[constrainedLot] = condVaR[constrainedLot] * (1 + exp(tail(estPars,1)))
      }
      if(Uni){
        if(!As){
          hypothesis = c(0,0,0,0) # intercept, absolute lag returns, autoregressive, phi
        } else{
          hypothesis = c(0,0,0,0,0)
        }
      } else{
        if(!As){
          hypothesis = c(0,0,0,0,0)
        } else{
          hypothesis = c(0,0,0,0,0,0)
        }
      }
      
      if(GetSe){
        resid = (y - condVaR)/abs(condVaR)
        parSim = matrix(0,nrow = length(estPars), ncol = GetSeSim)
        for(i in 1:GetSeSim){
          residSim = sample(resid,size = length(resid),replace = TRUE)
          xSim = sample(x, size = length(x),replace = TRUE)
          NegResidMean = mean(residSim[residSim < 0])
          ySim = cavSim(estPars, residSim, NegResidMean, y, condVaR, condES, As, Uni, xSim)
          ySim_filter = forecast::Arima(ySim,order = meanFit$arma[1:3])
          condMeanSim <- as.numeric(ySim_filter$fitted)
          estSim = .solverSwitch_cav(solver = MainSolver, pars = estPars, fun = objFunAL_cav, control = fitcontrol,
                                     As = As, y = ySim, x = xSim, empQuant = empQuant, Uni = Uni,
                                     condMean = condMeanSim, q = q,lb = lb, ub = ub)
          parSim[,i] = estSim$par
        }
        pval = rowMeans(abs(parSim - rowMeans(parSim,na.rm = TRUE) + hypothesis) > matrix(rep(abs(estPars),GetSeSim),nrow = length(estPars)),na.rm = TRUE)
      } else{
        pval = rep(NA,length(estPars))
      }
      ViolateRate <- (sum(y < condVaR))/(length(y))
      out = list(estPars = estPars, pval = pval, yLowFreq = y, yDate = yDate, condES = condES,ovlap = ovlap,
                 condVaR = condVaR, quantile = q, Uni = Uni, Solvers = c(MainSolver,SecondSolver),
                 conv = 0, simpleRet = simpleRet,As = As, meanFit = meanFit,
                 empQuant = empQuant, multiSol = multiSol, forecastLength = forecastLength,horizon = horizon,ViolateRate = ViolateRate)
      }
    return(out)
  }
    
  #--------------------------------------
  # Forecasting function
  #--------------------------------------
  #' @export VarEs_jointAL_cav_for
  VarEs_jointAL_cav_for <- function(object, y, yDate, x = NULL, xDate = NULL){
    #----- Recall the arguments ----
    estPars = object$estPars; pval = object$pval; q = object$q; Uni = object$Uni 
    simpleRet = object$simpleRet;As = object$As; meanFit = object$meanFit
    empQuant = object$empQuant; forecastLength = object$forecastLength;ovlap = object$ovlap
    horizon = object$horizon;
    #--- Form the single period setting for the return series
    nobs <- length(y)
    nobsShort = nobs-horizon+1
    yLowFreq = matrix(NA,ncol = 1, nrow = nobsShort)
    yDateLowFreq = matrix(NA,ncol = 1, nrow = nobsShort)
    for(t in 1:(nobs - horizon + 1)){
      yLowFreq[t,1] = sum(y[t:(t+horizon-1)]);
      yDateLowFreq[t,1] = yDate[t];
    }
    if(!ovlap){
      yLowFreq = yLowFreq[seq(1,nobsShort,horizon),1]
      yDateLowFreq = yDateLowFreq[seq(1,nobsShort,horizon),1]
    }
    if(simpleRet){
      y = exp(yLowFreq) - 1
    } else{
      y = yLowFreq
    }
    yDate = as.Date.numeric(yDateLowFreq,origin = "1970-01-01")
    if(is.null(x)){
      x = abs(y)
    }else {
      if(length(x) != length(y)){
        stop("\nCAviaR-->error: The explanatory variable should have the same length as the objective variable... \n")
      }
    }
    if(is.null(xDate)) xDate = yDate
    #----- Compute condVaR and condES and return
    VaRES = condVaRES_cav(params = estPars, yr = y, Xr = x, empQuant = empQuant, As = As, Uni = Uni)
    condVaR = VaRES$VaR
    condES = VaRES$ES
    ViolateRate <- (sum(y < condVaR))/(length(y))
    OSout = list(estPars = estPars, yLowFreq = y, yDate = yDate, condES = condES,
               condVaR = condVaR, quantile = q, Uni = Uni,
               simpleRet = simpleRet,As = As, empQuant = empQuant, forecastLength = forecastLength,
               horizon = horizon, ViolateRate = ViolateRate)
    out = list(InSample = object, OutSample = OSout)
  return(out)
}
