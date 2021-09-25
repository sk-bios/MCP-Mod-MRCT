library(DoseFinding)
library(gridExtra)

# MCP MoD example of Thomas 2017
dose = c( 0, 0.025, 0.125, 0.250, 0.5, 1.0)
maxeff = 0.4
alpha = 0.05 # significance level of one-sided upper

DFMRCT = function(rep,rep2,sigma,tau1,tau2,tau3,dose,maxeff,s,alpha,pi = 0.5){
  ## define candidate dose-response shapes
  models <- Mods(linear = NULL, emax = c(0.1,0.2,0.3),
                 quadratic = c(-0.833, -0.714),
                 doses = dose,maxEff = maxeff)
  altmodels <- Mods(emax = c(0.13),doses = dose,maxEff = maxeff)

  theta = altmodels$emax
  e0 = theta[1]
  eMax = theta[2]
  ed50 = theta[3]
  
  r = length(dose)
  
  n0 = c(80,40,40,40,40,80) # Total sample size by dose level
  n1 = c(60,60,60,60,60,60) # Total sample size by dose level
  
  f1 = rep(1,s)
  f1 = f1/sum(f1)
  f2 = c(rep(0.75,s/2),rep(1,s/2))
  f2 = f2/sum(f2)
  f3 = c(rep(2.5,s/2),rep(1,s/2))
  f3 = f3/sum(f3)

  if (rep == 1){n = n0}
  if (rep == 2){n = n1}

  if (rep2 == 1){f = f1}
  if (rep2 == 2){f = f2}
  if (rep2 == 3){f = f3}

  # Matrix Sigma
  A = array(0,dim = c(r,r,s))
  for (i in 1:s){
    for (j1 in 1:r){
      for (j2 in 1:r){
        if (j1 == j2){tmp = sigma^2 / n[j1] / f[i] + tau1^2 + (dose[j1] / (ed50 + dose[j1]))^2 * tau2^2 + (eMax * ed50 *dose[j1] / (ed50 + dose[j1])^2)^2 * tau3^2
        }else{tmp = tau1^2 + (dose[j1] / (ed50 + dose[j1])) * (dose[j2] / (ed50 + dose[j2]))* tau2^2 +
          (eMax*ed50 *dose[j1] / (ed50 + dose[j1])^2) * (eMax*ed50 *dose[j2] / (ed50 + dose[j2])^2) * tau3^2}
        
        A[j1,j2,i] = tmp
      }
    }
    
    if (i == 1){Sprime = f[i]^2 * A[,,i]
    }else{Sprime = Sprime + f[i]^2 * A[,,i]}
  }

  S = diag(sigma^2 / n)
  
  LL = optContr(models, w = n)
  contMat = LL$contMat
  
  muMat <- getResp(altmodels)
  covMat <- t(contMat) %*% S %*% contMat
  den <- sqrt(diag(covMat))
  corMat <- cov2cor(covMat)
  deltaMat <- t(contMat) %*% muMat
  deltaMat <- deltaMat/den
  critV = qmvnorm(1 - alpha,tail = c("lower.tail"),corr = corMat)$quantile
  power = powMCT(contMat = LL, alpha = alpha, altModels = altmodels, S = S, df = Inf, alternative = c("one.sided"))

  Is <- diag(sqrt(1/diag(covMat)))
  Vprime = Is %*% t(contMat) %*% Sprime %*% contMat %*% Is
  powerr = (1 - pmvt(lower = rep(-Inf,nrow(Vprime)),upper = rep(critV,nrow(Vprime)),delta = as.vector(deltaMat),sigma = Vprime,df = Inf)) * 100

  # ED 0.8
  muMat = getResp(altmodels)
  set.seed(1)
  rdat = rmvnorm(10000,muMat,Sprime)
  for (i in 1:nrow(rdat)){
    fitmodel = fitMod(dose,rdat[i,],S = S, model = "emax",bnds = c(0.001, 1.5), type = "general")
    yhat = emax(dose,coef(fitmodel)[1],coef(fitmodel)[2],coef(fitmodel)[3])
    ytrue = emax(dose,altmodels$emax[1],altmodels$emax[2],altmodels$emax[3])
    tmpED80 = ED(fitmodel,0.8)
    tmpsqdifm = mean((yhat - ytrue)^2)
    if (i == 1){ED80 = tmpED80; sqdifm = tmpsqdifm
    }else{ED80 = c(ED80,tmpED80); sqdifm = c(sqdifm,tmpsqdifm)}
  }
  ED80true = ED(altmodels,0.8)
  MSE = mean((ED80 - ED80true)^2)
  sqdifmm = mean(sqrt(sqdifm))

  # Consistency assessment
  muMat <- getResp(altmodels)
  num.m1 = (1 - pi) * t(LL$contMat) %*% muMat
  den.m1.1 = (1 - pi * f[1])^2 * diag(t(LL$contMat) %*% A[,,1] %*% LL$contMat)
  den.m1.2 = rep(0,length(den.m1.1))
  for (i in 2:s){
    den.m1.2 = den.m1.2 + pi^2 * f[i]^2 * diag(t(LL$contMat) %*% A[,,i] %*% LL$contMat)
  }
  m1= pnorm(num.m1 / sqrt(den.m1.1 + den.m1.2)) * 100
  
  m2 = 1
  for (i in 1:s){
    m2 = m2 * pnorm(t(LL$contMat) %*% muMat / sqrt(diag(t(LL$contMat) %*% A[,,i] %*% LL$contMat)))    
  }
  m2 = m2 *100
  
  out = data.frame(
    ncase = rep,
    fcase = rep2,
    sigma = sigma,
    tau1 = tau1,
    tau2 = tau2,
    tau3 = tau3,
    maxeff = maxeff,
    s = s,
    alpha = alpha,
    power = powerr,
    m1 = t(m1),
    m1.ave = mean(m1),
    m2 = t(m2),
    m2.ave = mean(m2),
    MSEED80 = MSE,
    sqdifmm
  )
  return(out)
}

slist = c(2,4,8)
replist = c(1,2)
rep2list = c(1,2,3)
for (i in 1:length(slist)){
  for (j in 1:length(replist)){
    for (k in 1:length(rep2list)){
      s = slist[i]
      rep = replist[j]
      rep2 = rep2list[k]
      
      out0 = DFMRCT(rep,rep2,1,0,0,0,dose,0.4,s,0.05)
      
      out11 = DFMRCT(rep,rep2,1,0.1,0,0,dose,0.4,s,0.05)
      out12 = DFMRCT(rep,rep2,1,0.2,0,0,dose,0.4,s,0.05)
      out13 = DFMRCT(rep,rep2,1,0.3,0,0,dose,0.4,s,0.05)
      out14 = DFMRCT(rep,rep2,1,0.4,0,0,dose,0.4,s,0.05)
      
      out211 = DFMRCT(rep,rep2,1,0,0.2,0,dose,0.4,s,0.05)
      out212 = DFMRCT(rep,rep2,1,0,0.2,0.4,dose,0.4,s,0.05)
      out213 = DFMRCT(rep,rep2,1,0,0.2,0.8,dose,0.4,s,0.05)
      out221 = DFMRCT(rep,rep2,1,0,0.4,0,dose,0.4,s,0.05)
      out222 = DFMRCT(rep,rep2,1,0,0.4,0.4,dose,0.4,s,0.05)
      out223 = DFMRCT(rep,rep2,1,0,0.4,0.8,dose,0.4,s,0.05)
      out231 = DFMRCT(rep,rep2,1,0,0.6,0,dose,0.4,s,0.05)
      out232 = DFMRCT(rep,rep2,1,0,0.6,0.4,dose,0.4,s,0.05)
      out233 = DFMRCT(rep,rep2,1,0,0.6,0.8,dose,0.4,s,0.05)
      out241 = DFMRCT(rep,rep2,1,0,0.8,0,dose,0.4,s,0.05)
      out242 = DFMRCT(rep,rep2,1,0,0.8,0.4,dose,0.4,s,0.05)
      out243 = DFMRCT(rep,rep2,1,0,0.8,0.8,dose,0.4,s,0.05)
      out2 = rbind(out211,out212,out213,out221,out222,out223,
                   out231,out232,out233,out241,out242,out243)
      
      out311 = DFMRCT(rep,rep2,1,0,0,0.2,dose,0.4,s,0.05)
      out312 = DFMRCT(rep,rep2,1,0,0.4,0.2,dose,0.4,s,0.05)
      out313 = DFMRCT(rep,rep2,1,0,0.8,0.2,dose,0.4,s,0.05)
      out321 = DFMRCT(rep,rep2,1,0,0,0.4,dose,0.4,s,0.05)
      out331 = DFMRCT(rep,rep2,1,0,0,0.6,dose,0.4,s,0.05)
      out332 = DFMRCT(rep,rep2,1,0,0.4,0.6,dose,0.4,s,0.05)
      out333 = DFMRCT(rep,rep2,1,0,0.8,0.6,dose,0.4,s,0.05)
      out341 = DFMRCT(rep,rep2,1,0,0,0.8,dose,0.4,s,0.05)
      out3 = rbind(out311,out312,out313,out321,
                   out331,out332,out333,out341)
      
      out = rbind(out0,out11,out12,out13,out14,
                  out2,out3)
      if (i == 1 & j == 1 & k == 1){output = out
      }else{output = rbind(output,out)}
    }
  }
}

output
