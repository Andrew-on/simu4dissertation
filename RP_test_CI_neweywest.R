RP_test_CI_neweywest <- function(reg.m, reg.y, data, npermute,m,tc1,tc2,tc3) {
  #tc1 = 0; tc2 = 0; tc3 = 0; m=mean(data$m)
  y_pred=as.numeric(predict(reg.y))
  m_pred=as.numeric(predict(reg.m))
  MY <- t(sapply(1:npermute, function(i) {
    # Randomization of Residuals : yres & mres
    y_residuals <- sample(residuals(reg.y), length(residuals(reg.y)), replace = FALSE)
    m_residuals <- sample(residuals(reg.m), length(residuals(reg.m)), replace = FALSE)
    
    # New data generating : ystar & mstar
    data$ystar <- y_pred +y_residuals
    data$mstar <- m_pred +m_residuals
    
    reg.m_perm <- try(lm(mstar ~ phase1+phase2+phase3, data = data),silent = T)

    reg.y_perm <- try(lm(ystar ~ phase1+phase2+phase3+m+phase1m+phase2m+phase3m, data = data),silent = T)
    
    if (length(reg.m_perm)>1 & length(reg.y_perm)>1) {
      mvcv <- NeweyWest(reg.m_perm, lag = 1, prewhite = FALSE)
      mcoef <- coeftest(reg.m_perm, data = data, vcov=mvcv);m_coef=c(mcoef[1,1],0, mcoef[2,1],0, mcoef[3,1],0, mcoef[4,1],0)
      yvcv <- NeweyWest(reg.y_perm, lag = 1, prewhite = FALSE)
      ycoef <- coeftest(reg.y_perm, data = data, vcov=yvcv);y_coef=c(ycoef[1,1],0, ycoef[2,1],0, ycoef[3,1],0, ycoef[4,1],0, ycoef[5:8,1])
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
