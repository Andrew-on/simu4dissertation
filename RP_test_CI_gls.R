RP_test_CI_gls <- function(reg.m, reg.y, data, npermute,m,tc1,tc2,tc3) {

  ar=1; ma=0; ordr=ar+ma
  ## Pre-whitening & Back-transformation
  ### reg.y
  e.y <- resid(reg.y)#, type = "normalized")  # normalized residuals (already whitened)
  rho.y <- coef(reg.y$modelStruct$corStruct, unconstrained = FALSE)  # estimated AR(1) coefficient
  n.y <- length(e.y)  # number of residuals
  y_e_white <- numeric(n.y)
  y_e_white[1] <- sqrt(1 - rho.y^2) * e.y[1]
  for (t in 2:n.y) {
    y_e_white[t] <- e.y[t] - rho.y * e.y[t - 1]
  }
  ### reg.m
  e.m <- resid(reg.m)#, type = "normalized")  # normalized residuals (already whitened)
  rho.m <- coef(reg.m$modelStruct$corStruct, unconstrained = FALSE)  # estimated AR(1) coefficient
  n.m <- length(e.m)  # number of residuals
  m_e_white <- numeric(n.m)
  m_e_white[1] <- sqrt(1 - rho.m^2) * e.m[1]
  for (t in 2:n.m) {
    m_e_white[t] <- e.m[t] - rho.m * e.m[t - 1]
  }
  
  MY <- t(sapply(1:npermute, function(i) {
    # Randomization of Residuals : yres & mres
    y_e_white_star <- sample(y_e_white, length(y_e_white), replace = FALSE)
    m_e_white_star <- sample(m_e_white, length(m_e_white), replace = FALSE)
    ye_star <- numeric(n.y)
    ye_star[1] <- y_e_white_star[1] / sqrt(1 - rho.y^2)  # y* first residual
    for (t in 2:n.y) {
      ye_star[t] <- rho.y * ye_star[t - 1] + y_e_white_star[t]  # recursive reconstruction
    }
    me_star <- numeric(n.m)
    me_star[1] <- m_e_white_star[1] / sqrt(1 - rho.m^2)  # m* first residual
    for (t in 2:n.m) {
      me_star[t] <- rho.m * me_star[t - 1] + m_e_white_star[t]  # recursive reconstruction
    }

    # New data generating : ystar & mstar
    data$ystar <- predict(reg.y) + ye_star
    data$mstar <- predict(reg.m) + me_star
    
    reg.m_perm <- try(nlme::gls(mstar ~ phase1+phase2+phase3, data=data, correlation = corARMA(form = ~ 1,p=ar,q=ma),method="REML"),silent = T)    
    #m_perm_coef <- coef(m_perm)
    reg.y_perm <- try(nlme::gls(ystar ~ phase1+phase2+phase3+m+phase1m+phase2m+phase3m, data=data, correlation = corARMA(form = ~ 1,p=ar,q=ma),method="REML"),silent = T)
    #y_perm_coef <- coef(y_perm)
    if (length(reg.m_perm)>1 & length(reg.y_perm)>1) {
      mcoef=reg.m_perm$coefficients;m_coef=c(mcoef[1],0, mcoef[2],0, mcoef[3],0, mcoef[4],0)
      ycoef=reg.y_perm$coefficients;y_coef=c(ycoef[1],0, ycoef[2],0, ycoef[3],0, ycoef[4],0, ycoef[5:8])
      c(m_coef,y_coef)
      }else{ c(rep(NA,20)) }
  }))
  colnames(MY)<-c('im', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7','iy', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'b', 'h1', 'h2', 'h3')
  PNIE.B1 = MY[,'b']*(MY[,'a2']+MY[,'a3']*tc1)
  PNIE.A2 = (MY[,'b']+MY[,'h2'])*((MY[,'a2']+MY[,'a3']*tc1)-(MY[,'a4']+MY[,'a5']*tc2))
  PNIE.B2 = (MY[,'b']+MY[,'h2'])*((MY[,'a6']+MY[,'a7']*tc3)-(MY[,'a4']+MY[,'a5']*tc2))
  TNIE.B1 = (MY[,'b']+MY[,'h1'])*(MY[,'a2']+MY[,'a3']*tc1)
  TNIE.A2 = (MY[,'b']+MY[,'h1'])*((MY[,'a2']+MY[,'a3']*tc1)-(MY[,'a4']+MY[,'a5']*tc2))
  TNIE.B2 = (MY[,'b']+MY[,'h3'])*((MY[,'a6']+MY[,'a7']*tc3)-(MY[,'a4']+MY[,'a5']*tc2))
  PNDE.B1 = (MY[,'c2']+MY[,'c3']*tc1)+MY[,'h1']*(MY[,'im']+MY[,'a1']*t1)
  PNDE.A2 = (MY[,'c2']+MY[,'c3']*tc1)-(MY[,'c4']+MY[,'c5']*tc2)+(MY[,'h1']-MY[,'h2'])*((MY[,'im']+MY[,'a1']*t2)+(MY[,'a4']+MY[,'a5']*tc2))
  PNDE.B2 = (MY[,'c6']+MY[,'c7']*tc3)-(MY[,'c4']+MY[,'c5']*tc2)-(MY[,'h2']-MY[,'h3'])*((MY[,'im']+MY[,'a1']*t2)+(MY[,'a4']+MY[,'a5']*tc2))
  TNDE.B1 = (MY[,'c2']+MY[,'c3']*tc1)+MY[,'h1']*(MY[,'im']+MY[,'a1']*t1+MY[,'a2']+MY[,'a3']*tc1)
  TNDE.A2 = (MY[,'c2']+MY[,'c3']*tc1)-(MY[,'c4']+MY[,'c5']*tc2)+(MY[,'h1']-MY[,'h2'])*((MY[,'im']+MY[,'a1']*t1)+(MY[,'a2']+MY[,'a3']*tc1))
  TNDE.B2 = (MY[,'c6']+MY[,'c7']*tc3)-(MY[,'c4']+MY[,'c5']*tc2)-(MY[,'h2']-MY[,'h3'])*((MY[,'im']+MY[,'a1']*t3)+(MY[,'a6']+MY[,'a7']*tc3))
  CDE.B1 =  (MY[,'c2']+MY[,'c3']*tc1)+MY[,'h1']*m
  CDE.A2 =  (MY[,'c2']+MY[,'c3']*tc1)-(MY[,'c4']+MY[,'c5']*tc2)+MY[,'h1']*m-MY[,'h2']*m
  CDE.B2 =  (MY[,'c6']+MY[,'c7']*tc3)-(MY[,'c4']+MY[,'c5']*tc2)-MY[,'h2']*m+MY[,'h3']*m
  TE.B1 = PNIE.B1 + TNDE.B1
  TE.A2 = PNIE.A2 + TNDE.A2
  TE.B2 = PNIE.B2 + TNDE.B2
  RP.estimates = cbind(PNIE.B1, PNIE.A2, PNIE.B2,
                       TNIE.B1, TNIE.A2, TNIE.B2,
                       PNDE.B1, PNDE.A2, PNDE.B2,
                       TNDE.B1, TNDE.A2, TNDE.B2,
                       CDE.B1, CDE.A2, CDE.B2,
                       TE.B1, TE.A2, TE.B2)
  RP.CI=t(apply(RP.estimates, 2, FUN=function(effect_perm){
    quantile(effect_perm, probs = c(0.025, 0.975),na.rm = TRUE)
  }))
  colnames(RP.CI) <- c('RP95LCL','RP95UCL')
  RP.CI=as.data.frame(RP.CI)
}
