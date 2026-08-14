RP_test_CI_orcutt <- function(reg.m, reg.y, data, npermute,m,tc1,tc2,tc3) {
  #tc1 = 0; tc2 = 0; tc3 = 0; m=mean(data$m)
  library(orcutt)
  y_pred=as.numeric(reg.y$fitted.values)
  m_pred=as.numeric(reg.m$fitted.values)
  ## Pre-whitening & Back-transformation
  rho.y <- reg.y$rho
  rho.m <- reg.m$rho
  
  eps_y <- residuals(reg.y)
  eps_m <- residuals(reg.m)
  
  u_y <- eps_y[-1] - rho.y * eps_y[-length(eps_y)]
  u_m <- eps_m[-1] - rho.m * eps_m[-length(eps_m)]
  
  MY <- t(sapply(1:npermute, function(i) {
    # Randomization of Residuals : yres & mres
    u_y_perm <- sample(u_y, length(u_y), replace = FALSE)
    u_m_perm <- sample(u_m, length(u_m), replace = FALSE)
    
    eps_y_perm <- numeric(length(eps_y))
    eps_m_perm <- numeric(length(eps_m))
    
    # initial value
    eps_y_perm[1] <- eps_y[1]
    eps_m_perm[1] <- eps_m[1]
    
    for (t in 2:length(eps_y)) {
      eps_y_perm[t] <- rho.y * eps_y_perm[t-1] + u_y_perm[t-1]
      eps_m_perm[t] <- rho.m * eps_m_perm[t-1] + u_m_perm[t-1]
    }
    
    y_residuals <- eps_y_perm
    m_residuals <- eps_m_perm
    # New data generating : ystar & mstar
    data$ystar <- y_pred +y_residuals
    data$mstar <- m_pred +m_residuals
    
    reg.m_perm <- try(orcutt::cochrane.orcutt(lm(mstar ~ phase1+phase2+phase3, data = data), max.iter=1e06),silent = T)
    reg.y_perm <- try(orcutt::cochrane.orcutt(lm(ystar ~ phase1+phase2+phase3+m+phase1m+phase2m+phase3m, data = data), max.iter=1e06),silent = T)
    #y_perm_coef <- coef(y_perm)
    if (length(reg.m_perm)>1 & length(reg.y_perm)>1) {
      mcoef = c(reg.m_perm$coefficients);m_coef=c(mcoef[1],0, mcoef[2],0, mcoef[3],0, mcoef[4],0)
      ycoef = c(reg.y_perm$coefficients);y_coef=c(ycoef[1],0, ycoef[2],0, ycoef[3],0, ycoef[4],0, ycoef[5:8])
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
