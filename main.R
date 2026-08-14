load("./SIMU_dataset.Rdata")
## Main Simulation
library(pbapply)
t0=Sys.time()
system.time(
  simu_summary <- pblapply(simu_dataset, function(SIMU_dataset){
    coefs.m = c('ar1','im', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7')
    coefs.y = c('ar1','iy', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'b', 'h1', 'h2', 'h3')
    
    ###-----------------------------------------------------------###
    
    library(parallel)
    num_cores <- parallel::detectCores(logical = FALSE)
    cl <- parallel::makeCluster(num_cores, type="PSOCK")
    
    ###-----------------------------------------------------------###
    ## Functions loading ##
    source("MCDLT_CI.R")          # MC & DLT
    source("RP_test_CI_gls.R")       # RP(1)
    source("RP_test_CI_arima.R")     # RP(2)
    source("RP_test_CI_dbfit.R")     # RP(3)
    source("RP_test_CI_neweywest.R") # RP(4)
    source("RP_test_CI_ols.R")       # RP(5)
    source("RP_test_CI_orcutt.R")    # RP(6)
    clusterExport(cl, varlist = c("MCDLT_CI", "RP_test_CI_gls", "RP_test_CI_arima", "RP_test_CI_dbfit", "RP_test_CI_neweywest", "RP_test_CI_ols", "RP_test_CI_orcutt"))
    
    ### 1. nlme::gls ## REML ###---------------------------------------###
    clusterEvalQ(cl,library(nlme))
    ar=1; ma=0; ordr=ar+ma
    gls.REML <- parSapply(cl,SIMU_dataset,function(dat){
      reg.m <- try(nlme::gls(m ~ time+phase1+phase1time+phase2+phase2time+phase3+phase3time, data=dat, correlation = corARMA(form = ~ 1,p=ar,q=ma),method="REML"),silent = T)
      reg.y <- try(nlme::gls(y ~ time+phase1+phase1time+phase2+phase2time+phase3+phase3time+m+phase1m+phase2m+phase3m, data=dat, correlation = corARMA(form = ~ 1,p=ar,q=ma),method="REML"),silent = T)
      
      if (length(reg.m)>1 & length(reg.y)>1) {
        m.coef = c(coef(reg.m$modelStruct$corStruct,unconstrained=FALSE), reg.m$coefficients)
        m.covm = reg.m$varBeta
        y.coef = c(coef(reg.y$modelStruct$corStruct,unconstrained=FALSE), reg.y$coefficients)
        y.covm = reg.y$varBeta
        ## CIs
        m <<- mean(dat[c(1:(which(dat$phase1==1)[1]-1)),'m'])
        tc1<<-0; t1 <<- min(which(dat$phase1 != tc1))
        tc2<<-0; t2 <<- min(which(dat$phase2 != tc2))
        tc3<<-0; t3 <<- min(which(dat$phase3 != tc3))
        MCDLT.CI <- MCDLT_CI(m.coef=m.coef,m.covm=m.covm, y.coef=y.coef, y.covm=y.covm, m=m,tc1=tc1,tc2=tc2,tc3=tc3,t1=t1,t2=t2,t3=t3)
        RP.CI <- RP_test_CI_gls(reg.m=reg.m, reg.y=reg.y,data=dat, npermute=500, m=m,tc1=tc1,tc2=tc2,tc3=tc3)
        ## Output format
        m.coef <- as.data.frame(t(m.coef));colnames(m.coef) <- c(coefs.m)
        y.coef <- as.data.frame(t(y.coef));colnames(y.coef) <- c(coefs.y)
        m.covm <- as.data.frame(m.covm)
        y.covm <- as.data.frame(y.covm)
        return( c(list(m.coef),list(m.covm),list(y.coef),list(y.covm),list(cbind(MCDLT.CI,RP.CI))) )
        
      }else{c(list(rep("Error",2))) }
    }, simplify = F)
    
    ### 2. stats::arima ## ML ###---------------------------------------###
    clusterEvalQ(cl,{
      library(stats)
      library(forecast)})
    
    arima.ML <- parSapply(cl,SIMU_dataset,function(dat){
      reg.m <- try(stats::arima(dat$m,order=c(1,0,0),xreg=dat[,c('time','phase1','phase1time','phase2','phase2time','phase3','phase3time')],method=c("ML")),silent = T)
      reg.y <- try(stats::arima(dat$y,order=c(1,0,0),xreg=dat[,c('time','phase1','phase1time','phase2','phase2time','phase3','phase3time','m','phase1m','phase2m','phase3m')],method=c("ML")),silent = T)
      
      if (length(reg.m)>1 & length(reg.y)>1) {
        m.coef = c(reg.m$coef)
        m.covm = reg.m$var.coef[-1,-1]
        y.coef = c(reg.y$coef)
        y.covm = reg.y$var.coef[-1,-1]
        ## CIs
        m <<- mean(dat[c(1:(which(dat$phase1==1)[1]-1)),'m'])
        tc1<<-0; t1 <<- min(which(dat$phase1 != tc1))
        tc2<<-0; t2 <<- min(which(dat$phase2 != tc2))
        tc3<<-0; t3 <<- min(which(dat$phase3 != tc3))
        #library(forecast) ## fitted() ##
        y_pred=as.numeric(fitted(reg.y))
        m_pred=as.numeric(fitted(reg.m))
        MCDLT.CI <- MCDLT_CI(m.coef=m.coef,m.covm=m.covm, y.coef=y.coef, y.covm=y.covm, m=m,tc1=tc1,tc2=tc2,tc3=tc3,t1=t1,t2=t2,t3=t3)
        RP.CI <- RP_test_CI_arima(reg.m=reg.m, reg.y=reg.y,data=dat, npermute=500, m=m,tc1=tc1,tc2=tc2,tc3=tc3, m_pred, y_pred)
        ## Output format
        m.coef <- as.data.frame(t(m.coef));colnames(m.coef) <- c(coefs.m)
        y.coef <- as.data.frame(t(y.coef));colnames(y.coef) <- c(coefs.y)
        m.covm <- as.data.frame(m.covm)
        y.covm <- as.data.frame(y.covm)
        return( c(list(m.coef),list(m.covm),list(y.coef),list(y.covm),list(cbind(MCDLT.CI,RP.CI))) )
      }else{c(list(rep("Error",2))) }
    }, simplify = F)
    
    ### 3. DBfit::dbfit ## OLS ###-------------------------------------###
    clusterEvalQ(cl,library(DBfit))
    
    dbfit.OLS <- parSapply(cl,SIMU_dataset,function(dat){
      set.seed(33620)
      dat$ones=1
      reg.m <- try(DBfit::dbfit(x=dat[,c('ones','time','phase1','phase1time','phase2','phase2time','phase3','phase3time')], y=dat[,"m"], arp=1, method="OLS"),silent = T)
      reg.y <- try(DBfit::dbfit(x=dat[,c('ones','time','phase1','phase1time','phase2','phase2time','phase3','phase3time','m','phase1m','phase2m','phase3m')], y=dat[,"y"], arp=1, method="OLS"),silent = T)
      if (length(reg.m)>1 & length(reg.y)>1) {
        m.coef = c(reg.m$adjar,reg.m$coef)
        m.covm = reg.m$betacov
        y.coef = c(reg.y$adjar,reg.y$coef)
        y.covm = reg.y$betacov
        ## CIs
        m <<- mean(dat[c(1:(which(dat$phase1==1)[1]-1)),'m'])
        tc1<<-0; t1 <<- min(which(dat$phase1 != tc1))
        tc2<<-0; t2 <<- min(which(dat$phase2 != tc2))
        tc3<<-0; t3 <<- min(which(dat$phase3 != tc3))
        MCDLT.CI <- MCDLT_CI(m.coef=m.coef,m.covm=m.covm, y.coef=y.coef, y.covm=y.covm, m=m,tc1=tc1,tc2=tc2,tc3=tc3,t1=t1,t2=t2,t3=t3)
        RP.CI <- RP_test_CI_dbfit(reg.m=reg.m, reg.y=reg.y,data=dat, npermute=500, m=m,tc1=tc1,tc2=tc2,tc3=tc3)
        ## Output format
        m.coef <- as.data.frame(t(m.coef));colnames(m.coef) <- c(coefs.m)
        y.coef <- as.data.frame(t(y.coef));colnames(y.coef) <- c(coefs.y)
        m.covm <- as.data.frame(m.covm)
        y.covm <- as.data.frame(y.covm)
        return( c(list(m.coef),list(m.covm),list(y.coef),list(y.covm),list(cbind(MCDLT.CI,RP.CI))) )
      }else{c(list(rep("Error",2))) }
    }, simplify = F)
    
    ### 4.Sandwich::NeweyWest ## OLS ###-----------------------------------###
    clusterEvalQ(cl,{
      library(stats)
      library(sandwich)
      library(lmtest)})
    neweywest.OLS <- parSapply(cl,SIMU_dataset,function(dat){
      reg.m <- try(lm(m ~ time+phase1+phase1time+phase2+phase2time+phase3+phase3time, data = dat ),silent = T)
      reg.y <- try(lm(y ~ time+phase1+phase1time+phase2+phase2time+phase3+phase3time+m+phase1m+phase2m+phase3m, data = dat ),silent = T)
      if (length(reg.m)>1 & length(reg.y)>1) {
        mvcv <- NeweyWest(reg.m, lag = 1, prewhite = FALSE)
        mcoef <- coeftest(reg.m, data = dat, vcov=mvcv)
        yvcv <- NeweyWest(reg.y, lag = 1, prewhite = FALSE)
        ycoef <- coeftest(reg.y, data = dat, vcov=yvcv)
        
        m.coef = c(NA, mcoef[,1])
        m.covm = mvcv
        y.coef = c(NA, ycoef[,1])
        y.covm = yvcv
        ## CIs
        m <<- mean(dat[c(1:(which(dat$phase1==1)[1]-1)),'m'])
        tc1<<-0; t1 <<- min(which(dat$phase1 != tc1))
        tc2<<-0; t2 <<- min(which(dat$phase2 != tc2))
        tc3<<-0; t3 <<- min(which(dat$phase3 != tc3))
        MCDLT.CI <- MCDLT_CI(m.coef=m.coef,m.covm=m.covm, y.coef=y.coef, y.covm=y.covm, m=m,tc1=tc1,tc2=tc2,tc3=tc3,t1=t1,t2=t2,t3=t3)
        RP.CI <- RP_test_CI_neweywest(reg.m=reg.m, reg.y=reg.y,data=dat, npermute=500, m=m,tc1=tc1,tc2=tc2,tc3=tc3)
        ## Output format
        m.coef <- as.data.frame(t(m.coef));colnames(m.coef) <- c(coefs.m)
        y.coef <- as.data.frame(t(y.coef));colnames(y.coef) <- c(coefs.y)
        m.covm <- as.data.frame(m.covm)
        y.covm <- as.data.frame(y.covm)
        return( c(list(m.coef),list(m.covm),list(y.coef),list(y.covm),list(cbind(MCDLT.CI,RP.CI)) ) ) 
      }else{c(list(rep("Error",2))) }
    }, simplify = F)  
    
    ### 5.stats::lm ## OLS ###---------------------------------------###
    clusterEvalQ(cl,library(stats))
    OLS <- parSapply(cl,SIMU_dataset,function(dat){
      reg.m <- try(lm(m ~ time+phase1+phase1time+phase2+phase2time+phase3+phase3time, data = dat),silent = T)
      reg.y <- try(lm(y ~ time+phase1+phase1time+phase2+phase2time+phase3+phase3time+m+phase1m+phase2m+phase3m, data = dat),silent = T)
      if (length(reg.m)>1 & length(reg.y)>1) {
        m.coef = c(NA, reg.m$coefficients); m.covm = vcov(reg.m)
        y.coef = c(NA, reg.y$coefficients); y.covm = vcov(reg.y)
        ## CIs
        m <<- mean(dat[c(1:(which(dat$phase1==1)[1]-1)),'m'])
        tc1<<-0; t1 <<- min(which(dat$phase1 != tc1))
        tc2<<-0; t2 <<- min(which(dat$phase2 != tc2))
        tc3<<-0; t3 <<- min(which(dat$phase3 != tc3))
        MCDLT.CI <- MCDLT_CI(m.coef=m.coef,m.covm=m.covm, y.coef=y.coef, y.covm=y.covm, m=m,tc1=tc1,tc2=tc2,tc3=tc3,t1=t1,t2=t2,t3=t3)
        RP.CI <- RP_test_CI_ols(reg.m=reg.m, reg.y=reg.y,data=dat, npermute=500, m=m,tc1=tc1,tc2=tc2,tc3=tc3)
        ## Output format
        m.coef <- as.data.frame(t(m.coef));colnames(m.coef) <- c(coefs.m)
        y.coef <- as.data.frame(t(y.coef));colnames(y.coef) <- c(coefs.y)
        m.covm <- as.data.frame(m.covm)
        y.covm <- as.data.frame(y.covm)
        return( c(list(m.coef),list(m.covm),list(y.coef),list(y.covm),list(cbind(MCDLT.CI,RP.CI))) )
      }else{c(list(rep("Error",2))) }
    }, simplify = F) 
    
    ## 6.orcutt::cochrane.orcutt ## FGLS ###------------------------------------###
    clusterEvalQ(cl,library(orcutt))
    orcutt.FGLS <- parSapply(cl,SIMU_dataset,function(dat){
      reg.m <- try(orcutt::cochrane.orcutt(lm(m ~ time+phase1+phase1time+phase2+phase2time+phase3+phase3time, data = dat), max.iter=1e06),silent = T)
      reg.y <- try(orcutt::cochrane.orcutt(lm(y ~ time+phase1+phase1time+phase2+phase2time+phase3+phase3time+m+phase1m+phase2m+phase3m, data = dat), max.iter=1e06),silent = T)
      if (length(reg.m)>1 & length(reg.y)>1) {
        m.coef = c(reg.m$rho, reg.m$coefficients)
        m.covm = diag(reg.m$std.error**2)
        y.coef = c(reg.y$rho, reg.y$coefficients)
        y.covm = diag(reg.y$std.error**2)
        ## CIs
        m <<- mean(dat[c(1:(which(dat$phase1==1)[1]-1)),'m'])
        tc1<<-0; t1 <<- min(which(dat$phase1 != tc1))
        tc2<<-0; t2 <<- min(which(dat$phase2 != tc2))
        tc3<<-0; t3 <<- min(which(dat$phase3 != tc3))
        MCDLT.CI <- MCDLT_CI(m.coef=m.coef,m.covm=m.covm, y.coef=y.coef, y.covm=y.covm, m=m,tc1=tc1,tc2=tc2,tc3=tc3,t1=t1,t2=t2,t3=t3)
        RP.CI <- RP_test_CI_orcutt(reg.m=reg.m, reg.y=reg.y,data=dat, npermute=500, m=m,tc1=tc1,tc2=tc2,tc3=tc3)
        ## Output format
        m.coef <- as.data.frame(t(m.coef));colnames(m.coef) <- c(coefs.m)
        y.coef <- as.data.frame(t(y.coef));colnames(y.coef) <- c(coefs.y)
        m.covm <- as.data.frame(m.covm)
        y.covm <- as.data.frame(y.covm)
        return( c(list(m.coef),list(m.covm),list(y.coef),list(y.covm),list(cbind(MCDLT.CI,RP.CI))) )
      }else{c(list(rep("Error",2))) }
    }, simplify = F) 
    
    stopCluster(cl)
    
    SIMU.summary <- list(gls = gls.REML, 
                         arima = arima.ML, 
                         dbfit = dbfit.OLS, 
                         nw = neweywest.OLS,
                         ols= OLS,
                         orcutt= orcutt.FGLS)
  })
)
Sys.time()-t0
save(simu_summary, file = "./SIMU_Summary.Rdata")
#load("./SIMU_Summary.Rdata")