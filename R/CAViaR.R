#--------------------------------
# MAIN FUNCTIONS: Estimate MIDAS quantile regression
#--------------------------------
#' @importFrom stats quantile
#' @export CAViaR
CAViaR <- function(y,yDate,x = NULL, xDate = NULL,q = 0.01,horizon = 10, ovlap = FALSE, numInitialsRand = 50000, 
                   empQuant = NULL,numInitials = 20, GetSe = TRUE, startPars = NULL,Params = NULL,
                   MainSolver = "neldermead",SecondSolver = "ucminf",As = FALSE, fitcontrol = list(rep = 10),
                   warn = TRUE, simpleRet = FALSE, Uni = TRUE, forecastLength = 0, constrained = FALSE, multiSol = TRUE){
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
  nobs = length(y)
  nobsShort = nobs-horizon+1
  yLowFreq = matrix(NA,ncol = 1, nrow = nobsShort)
  yDateLowFreq = matrix(NA,ncol = 1, nrow = nobsShort)
  for(t in 1:(nobs - horizon + 1)){
    yLowFreq[t,1] = sum(y[t:(t+horizon-1)])
    yDateLowFreq[t,1] = yDate[t]
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
  
  tol = 1e-10
  #----- Get the initial guess for the parameters-----
  if(is.null(Params)){
  if(is.null(startPars)){
    betaIni = GetIniParams_cav(y = y, x = x, q = q, As = As, empQuant = empQuant, Uni = Uni, 
                             numInitialsRand = numInitialsRand, numInitials = numInitials)
  } else{
    betaIni = matrix(rep(startPars,numInitials),nrow = numInitials, byrow = TRUE)
  }
  #----- Estimate the paramters -----------
  if(!Uni){
    if(As){
      lb = c(-Inf,-Inf,-Inf,-Inf,tol)
      ub = c(Inf, Inf, Inf, Inf, 1-tol)
    } else{
      lb = c(-Inf, -Inf, -Inf, tol)
      ub = c(Inf, Inf, Inf, 1-tol)
    }
  } else{
    if(As){
      lb = c(-Inf,-Inf,-Inf,tol)
      ub = c(Inf, Inf, Inf, 1-tol)
    } else{
      lb = c(-Inf, -Inf, tol)
      ub = c(Inf, Inf, 1-tol)
    }
  }
  
  sol = .sol_cav(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFun_cav,
                 y = y, x = x,As = As,empQuant = empQuant,Uni = Uni, q = q, lb = lb, ub = ub,
                 control = fitcontrol,warn = warn,multiSol = multiSol)
  for(i in 1:3){
    if(sol$convergence == 0 && (sum(sol$par >= lb) + sum(sol$par <= ub)) == (2 * length(lb))){
      break
    } else{
      if(sol$convergence == 1){
      sol = .sol_cav(MainSolver = SecondSolver,SecondSolver = MainSolver,betaIni = betaIni,fun = objFun_cav,
                     y = y, x = x,As = As,empQuant = empQuant,Uni = Uni, q = q, lb = lb, ub = ub,
                     control = fitcontrol,warn = warn,multiSol = multiSol)
      } else if(constrained){
        sol = .sol_cav(MainSolver = "bobyqa",SecondSolver = SecondSolver,betaIni = betaIni,fun = objFun_cav,
                       y = y, x = x,As = As,empQuant = empQuant,Uni = Uni, q = q, lb = lb, ub = ub,
                       control = fitcontrol,warn = warn,multiSol = FALSE)
      } else{
        break
      }
    }
  }
  #----- Preparing outputs and computing standard errors for estimated paramters-----
  # For the CAViaR model, the standard errors are calculated using the code from Engle and Manganelli
  
  estPars = sol$par
  convergeFlag = sol$convergence
  } else{
    estPars = Params
    convergeFlag = 0
  }
  if(convergeFlag == 1){
    warning("\nBoth Solvers failed to converge, try with other available solvers...\n")
    out = list(estPars = NaN, betaIni = betaIni, y = y, x = x,As = As,empQuant = empQuant,Uni = Uni, q = q,lb = lb, ub = ub,multiSol = multiSol,
               control = fitcontrol,warn = warn, empQuant = empQuant, simpleRet = simpleRet, forecastLength = forecastLength,
               ovlap = ovlap, conv = 1, GetSe = GetSe, fitcontrol = fitcontrol, yDate = yDate, horizon = horizon,constrained = constrained)
  } else{
    condVaR = condVaR_cav(params = estPars,yr = y,Xr = x,As = As,empQuant = empQuant,Uni = Uni)
    if(q < 0.5 && any(condVaR > 0)){
      # In case we have condVaR for q < 0.5, using this method, condES will be larger than 0, which does not make sence.
      # Hence, in this case, we constrained VaR to be slightly less than 0
      constrainedLot = which(condVaR > 0)
      condVaR[constrainedLot] = -0.0001
    }
    if(GetSe){
      VarCov = .VarCovarCaviaR(pars = estPars,As = As,y = y,x = x,q = q,condVaR = condVaR,Uni = Uni)
      if (inherits(VarCov, "try-error")) {
        warning("\nCAViaR : Cannot calculate VarCov...\n")
        stdErr = rep(NA,1,length(estPars))
        pval = rep(NA,1,length(estPars))
      } else{
        stdErr = sqrt(diag(VarCov))
        pval = pnorm(-abs(estPars)/stdErr)
      }
    } else{
      stdErr = rep(NA,1,length(estPars))
      pval = rep(NA,1,length(estPars))
    }
    ViolateRate <- (sum(y < condVaR))/(length(y))
    out = list(estPars = estPars, pval = pval, yLowFreq = y, yDate = yDate, condVaR = condVaR,ViolateRate = ViolateRate,
               quantile = q, As = As,Uni = Uni, simpleRet = simpleRet, Solvers = c(MainSolver,SecondSolver),horizon = horizon,
               multiSol = multiSol,conv = 0, empQuant = empQuant, forecastLength = forecastLength,ovlap = ovlap,constrained = constrained)  
    }
  return(out)
}

#--------------------------------------
# Forecasting function
#--------------------------------------
#' @export CAViaRfor

CAViaRfor <- function(object,y,yDate,x = NULL, xDate = NULL){
  #---- Get in-sample arguments ----
  estPars = object$estPars; horizon = object$horizon
  ovlap = object$ovlap; As = object$As
  Uni = object$Uni; simpleRet = object$simpleRet
  empQuant = object$empQuant
  #---- Form the single period setting for the return series
  nobs = length(y)
  nobsShort = nobs-horizon+1
  yLowFreq = matrix(NA,ncol = 1, nrow = nobsShort)
  yDateLowFreq = matrix(NA,ncol = 1, nrow = nobsShort)
  for(t in 1:(nobs - horizon + 1)){
    yLowFreq[t,1] = sum(y[t:(t+horizon-1)])
    yDateLowFreq[t,1] = yDate[t]
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

  #---- Get the Forecasts ----
  condVaR = condVaR_cav(params = estPars,yr = y,Xr = x,As = As,empQuant = empQuant,Uni = Uni)
  ViolateRate <- (sum(y < condVaR))/(length(y))
  OutSample = list(estPars = estPars, yLowFreq = y, yDate = yDate, condVaR = condVaR, ViolateRate = ViolateRate,
             quantile = object$quantile, As = As, Uni = Uni, simpleRet = simpleRet, horizon = horizon,
             empQuant = empQuant, ovlap = ovlap)
  out = list(InSample = object, OutSample = OutSample)
  return(out)
}

.VarCovarCaviaR <- function(pars, As, y, x, q, condVaR,Uni){
  T = length(y)
  resid <- y - condVaR
  SortedRes <- sort(abs(resid))
  if(q == 0.01){
    k = 40; bandwidth = SortedRes[40]
  }else if(q == 0.025){
    k = 50; bandwidth = SortedRes[50]
  } else if(q == 0.05){
    k = 60; bandwidth = SortedRes[60]
  } else{
    kk = median(abs(resid - median(resid)))
    hh = (T^(-1/3))*(qnorm(1 - 0.05/2)^2/3)*(((1.5*((dnorm(qnorm(q)))^2))/(2*(qnorm(q))^2 + 1)))^(1/3)
    bandwidth = kk * (qnorm(q + hh) - qnorm(q - hh))
  }
  D = matrix(0,ncol = length(pars), nrow = length(pars))
  A = D
  t = 0
  derivative1 = matrix(0,nrow = T, ncol = 1)
  derivative2 = matrix(0,nrow = T, ncol = 1)
  derivative3 = matrix(0,nrow = T, ncol = 1)
  derivative4 = matrix(0,nrow = T, ncol = 1)
  if(Uni){
    if(!As){
      gradient = matrix(0,nrow = T, ncol = length(pars))
      for(i in 2:T){
        derivative1[i] = 1 + pars[1] * derivative1[i-1]
        derivative2[i] = pars[2] * derivative2[i-1] + abs(y[i-1])
        derivative3[i] = pars[3] * derivative3[i-1] + condVaR[i-1]
        gradient[i,] = c(derivative1[i],derivative2[i],derivative3[i])
        A = A + gradient[i,] %*% t(gradient[i,])
        if(abs(resid[i]) <= bandwidth){
          t = t + 1
          D = D + gradient[i,] %*% t(gradient[i,])
        }
      }
    } else{
      gradient = matrix(0,nrow = T, ncol = length(pars))
      for(i in 2:T){
        derivative1[i] = 1 + pars[1] * derivative1[i-1]
        derivative2[i] = pars[2] * derivative2[i-1] + abs(y[i-1]) * (y[i-1] <= 0)
        derivative3[i] = pars[3] * derivative3[i-1] + abs(y[i-1]) * (y[i-1] > 0)
        derivative4[i] = pars[4] * derivative4[i-1] + condVaR[i-1]
        gradient[i,] = c(derivative1[i],derivative2[i],derivative3[i],derivative4[i])
        A = A + gradient[i,] %*% t(gradient[i,])
        if(abs(resid[i]) <= bandwidth){
          t = t + 1
        D = D + gradient[i,] %*% t(gradient[i,])
        }
      }
    }
  } else{
    derivative5 = matrix(0,nrow = T, ncol = 1)
    if(!As){
      gradient = matrix(0,nrow = T, ncol = length(pars))
      for(i in 2:T){
        derivative1[i] = 1 + pars[1] * derivative1[i-1]
        derivative2[i] = pars[2] * derivative2[i-1] + abs(y[i-1])
        derivative3[i] = pars[3] * derivative3[i-1] + x[i-1]
        derivative4[i] = pars[4] * derivative4[i-1] + condVaR[i-1]
        gradient[i,] = c(derivative1[i],derivative2[i],derivative3[i],derivative4[i])
        A = A + gradient[i,] %*% t(gradient[i,])
        if(abs(resid[i]) <= bandwidth){
          t = t + 1
          D = D + gradient[i,] %*% t(gradient[i,])
        }
      }
    } else{
      gradient = matrix(0,nrow = T, ncol = length(pars))
      for(i in 2:T){
        derivative1[i] = 1 + pars[1] * derivative1[i-1]
        derivative2[i] = pars[2] * derivative2[i-1] + abs(y[i-1]) * (y[i-1] <= 0)
        derivative3[i] = pars[3] * derivative3[i-1] + abs(y[i-1]) * (y[i-1] > 0)
        derivative4[i] = pars[4] * derivative4[i-1] + x[i-1]
        derivative5[i] = pars[5] * derivative5[i-1] + condVaR[i-1]
        gradient[i,] = c(derivative1[i],derivative2[i],derivative3[i],derivative4[i],derivative5[i])
        A = A + gradient[i,] %*% t(gradient[i,])
        if(abs(resid[i]) <= bandwidth){
          t = t + 1
          D = D + gradient[i,] %*% t(gradient[i,])
        }
      }
    }
  }
  tStdError = t
  A = A/T
  D = D/(2*bandwidth*T)
  VCmatrix = try((q * (1-q)) * (solve(D) %*% A %*% solve(D)) * 1/T,silent = TRUE)
  return(VCmatrix)
}

#--------------------------------------
# Resume function
#--------------------------------------
#' @export CAViaRresume

CAViaRresume <- function(object, MainSolver = "neldermead", SecondSolver = "bobyqa",multiSol = FALSE){
  if(is.na(match(MainSolver,c("bobyqa","ucminf","neldermead","solnp")))){
    stop("\nMidasQuantile-->error: available solvers are bobyqa, ucminf and neldermead... \n")
  }
  if(!is.null(SecondSolver)){
    if(is.na(match(SecondSolver,c("bobyqa","ucminf","neldermead")))){
      stop("\nMidasQuantile-->error: Available solvers are bobyqa, ucminf and neldermead... \n")
    }
  }
  #--------- Recall the arguments -----
  betaIni = object$betaIni; y = object$y; x = object$x ; As = object$As 
  empQuant = object$empQuant; Uni = object$Uni; q = object$q; yDate = object$yDate; horizon = object$horizon
  control = object$fitcontrol; warn = object$warn; empQuant = object$empQuant
  simpleRet = object$simpleRet; forecastLength = object$forecastLength
  ovlap = object$ovlap; GetSe = object$GetSe; fitcontrol = object$fitcontrol; constrained = object$constrained
  lb = object$lb ; ub = object$ub;
  #----- Estimate the paramters -----------
  sol = .sol_cav(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFun_cav,
                 y = y, x = x,As = As,empQuant = empQuant,Uni = Uni, q = q,
                 control = fitcontrol,warn = warn, multiSol = multiSol)
  for(i in 1:3){
    if(sol$convergence == 0 && (sum(sol$par >= lb) + sum(sol$par <= ub)) == (2 * length(lb))){
      break
    } else{
      if(sol$convergence == 1){
        sol = .sol_cav(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFun_cav,
                       y = y, x = x,As = As,empQuant = empQuant,Uni = Uni, q = q, lb = lb, ub = ub,
                       control = fitcontrol,warn = warn, multiSol = multiSol)
      } else if(constrained){
        sol = .sol_cav(MainSolver = "bobyqa",SecondSolver = SecondSolver,betaIni = betaIni,fun = objFun_cav,
                       y = y, x = x,As = As,empQuant = empQuant,Uni = Uni, q = q, lb = lb, ub = ub,
                       control = fitcontrol,warn = warn, multiSol = FALSE)
      } else{
        break
      }
    }
  }
  #----- Preparing outputs and computing standard errors for estimated paramters-----
  # For the CAViaR model, the standard errors are calculated using the code from Engle and Manganelli
  
  estPars = sol$par
  fval = sol$value
  convergeFlag = sol$convergence
  if(convergeFlag == 1){
    warning("\nBoth Solvers failed to converge, try with other available solvers...\n")
    out = list(estPars = NaN, betaIni = betaIni, y = y, x = x,As = As,empQuant = empQuant,Uni = Uni, q = q,
               control = fitcontrol,warn = warn, empQuant = empQuant, simpleRet = simpleRet, forecastLength = forecastLength,
               ovlap = ovlap, multiSol = multiSol, conv = 1)
  } else{
    condVaR = condVaR_cav(params = estPars,yr = y,Xr = x,As = As,empQuant = empQuant,Uni = Uni)
    if(GetSe){
      VarCov = .VarCovarCaviaR(pars = estPars,As = As,y = y,x = x,q = q,condVaR = condVaR,Uni = Uni)
      if (inherits(VarCov, "try-error")) {
        warning("\nCAViaR : Cannot calculate VarCov...\n")
        stdErr = rep(NA,1,length(estPars))
        pval = rep(NA,1,length(estPars))
      } else{
        stdErr = sqrt(diag(VarCov))
        pval = pnorm(-abs(estPars)/stdErr)
      }
    } else{
      stdErr = rep(NA,1,length(estPars))
      pval = rep(NA,1,length(estPars))
    }
    ViolateRate <- (sum(y < condVaR))/(length(y))
    out = list(estPars = estPars, pval = pval, yLowFreq = y, yDate = yDate, condVaR = condVaR,ViolateRate = ViolateRate,
               quantile = q, As = As,Uni = Uni, simpleRet = simpleRet, Solvers = c(MainSolver,SecondSolver),horizon = horizon,
               fval = fval, multiSol = multiSol, conv = 0, empQuant = empQuant, forecastLength = forecastLength,ovlap = ovlap)  
  }
  return(out)
  
}