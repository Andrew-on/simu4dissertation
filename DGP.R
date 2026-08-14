### DGP Function ###
SIMU_data_generator <- function(nsim, ss, ar, a1,a2,a3,a4,a5,a6,a7, b, c1,c2,c3,c4,c5,c6,c7, h1,h2,h3){
  library(stats)    
  #var_name <- c("time","phase1","timec1","phase2","timec2","phase3","timec3","em","m","ey","y")
  var_name <- c("time","phase1","phase1time","phase2","phase2time","phase3","phase3time","em","m","ey","y")
  dt = data.frame(matrix(ncol = length(var_name), nrow = ss ))
  colnames(dt) <- var_name
  
  dt[,'time'] = c(1:ss)
  ## Phase B1
  dt[,'phase1'] = ifelse(dt$time <= ss*2/4 & dt$time > ss/4, 1, 0)
  dt[,'phase1time'] = ifelse((dt$phase1*dt$time - floor(ss/4)) <= 0,0,dt$phase1*dt$time - floor(ss/4))
  ## Phase A2
  dt[,'phase2'] = ifelse(dt$time <= ss*3/4 & dt$time > ss/2, 1, 0)
  dt[,'phase2time'] = ifelse((dt$phase2*dt$time - floor(ss/2)) <= 0,0,dt$phase2*dt$time - floor(ss/2))
  ## Phase B2
  dt[,'phase3'] = ifelse(dt$time > ss*3/4, 1, 0)
  dt[,'phase3time'] = ifelse((dt$phase3*dt$time - floor(ss*3/4)) <=0,0,dt$phase3*dt$time - floor(ss*3/4))
  
  lapply(c(1:nsim),function(n){
    seed = sample.int(.Machine$integer.max, 1 , replace = FALSE)
    set.seed(seed)
    dt[,'em'] = stats::arima.sim(n = ss, list(ar = ar), sd = 1)
    dt[,'ey'] = stats::arima.sim(n = ss, list(ar = ar), sd = 1)
    dt[,'m'] = a1*dt$time + (a2 + a3*dt$phase1time)*dt$phase1 + (a4 + a5*dt$phase2time)*dt$phase2 + (a6 + a7*dt$phase3time)*dt$phase3 + dt$em
    dt[,'y'] = c1*dt$time + (c2 + c3*dt$phase1time)*dt$phase1 + (c4 + c5*dt$phase2time)*dt$phase2 + (c6 + c7*dt$phase3time)*dt$phase3 + (b + h1*dt$phase1 + h2*dt$phase2 + h3*dt$phase3)*dt$m + dt$ey
    dt[,'phase1m'] = dt$phase1*dt$m
    dt[,'phase2m'] = dt$phase2*dt$m
    dt[,'phase3m'] = dt$phase3*dt$m
    dt$seed = seed
    return(dt)
  })
  
}

### Parameter Pool ###
load("./SIMU_param_pool.Rdata")

### Data-Generating Process ###
simu_dataset <- lapply(c(1:dim(param.pool)[1]),function(pc){
  data.set <- SIMU_data_generator(
    nsim=param.pool[pc,'nsim'],
    ss = param.pool[pc,'ss'], 
    ar = param.pool[pc,'ar'], 
    a1 = param.pool[pc,'a1'], 
    a2 = param.pool[pc,'a2'], 
    a3 = param.pool[pc,'a3'], 
    a4 = param.pool[pc,'a4'], 
    a5 = param.pool[pc,'a5'], 
    a6 = param.pool[pc,'a6'], 
    a7 = param.pool[pc,'a7'], 
    b  = param.pool[pc,'b'], 
    c1 = param.pool[pc,'c1'], 
    c2 = param.pool[pc,'c2'], 
    c3 = param.pool[pc,'c3'], 
    c4 = param.pool[pc,'c4'], 
    c5 = param.pool[pc,'c5'], 
    c6 = param.pool[pc,'c6'], 
    c7 = param.pool[pc,'c7'],
    h1 = param.pool[pc,'h1'], 
    h2 = param.pool[pc,'h2'], 
    h3 = param.pool[pc,'h3']    )
})

save(simu_dataset, file = "./SIMU_dataset.Rdata")