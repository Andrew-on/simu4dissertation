## MC+DLT Testing Procedures
MCDLT_CI<-function(m.coef,m.covm, y.coef, y.covm, m, tc1, tc2, tc3,t1,t2,t3  ){
ar=1; ma=0; ordr=ar+ma
## Model M
coefs = c('im', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7')
m.coef=c(m.coef[1:2],0, m.coef[3],0, m.coef[4],0, m.coef[5],0)
for (i in 2:(length(coefs)+1) ) {
  var.name<-paste0('',coefs[i-1])
  assign(var.name, as.numeric(m.coef[i]))
}

rownames(m.covm) <- c("im", "a2", "a4", "a6")
colnames(m.covm) <- c("im", "a2", "a4", "a6")
### Zero-injection Model M
vcm <- matrix(0,nrow = length(coefs),ncol = length(coefs),
  dimnames = list(coefs, coefs)
)
vcm[rownames(m.covm), colnames(m.covm)] <- as.matrix(m.covm)
m.covm <- as.data.frame(vcm)

## Model Y
coefs = c('iy', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'b', 'h1', 'h2', 'h3')
y.coef=c(y.coef[1:2],0, y.coef[3],0, y.coef[4],0, y.coef[5],0, y.coef[6:9])
for (i in 2:(length(coefs)+1) ) {
  var.name<-paste0('',coefs[i-1])
  assign(var.name, as.numeric(y.coef[i]))
}

rownames(y.covm) <- c('iy', 'c2', 'c4', 'c6', 'b', 'h1', 'h2', 'h3')
colnames(y.covm) <- c('iy', 'c2', 'c4', 'c6', 'b', 'h1', 'h2', 'h3')
### Zero-injection Model Y
vcy <- matrix(0,nrow = length(coefs),ncol = length(coefs),
             dimnames = list(coefs, coefs)
)
vcy[rownames(y.covm), colnames(y.covm)] <- as.matrix(y.covm)
y.covm <- as.data.frame(vcy)

#tc1 = 0; tc2 = 0; tc3 = 0; 
#m = mean(dat$m)
#t1 = min(which(dat$phase1 != tc1))
#t2 = min(which(dat$phase2 != tc2)) 
#t3 = min(which(dat$phase3 != tc3))

PNIE.B1 = b*(a2+a3*tc1)
PNIE.A2 = (b+h2)*((a2+a3*tc1)-(a4+a5*tc2))
PNIE.B2 = (b+h2)*((a6+a7*tc3)-(a4+a5*tc2))

TNIE.B1 = (b+h1)*(a2+a3*tc1)
TNIE.A2 = (b+h1)*((a2+a3*tc1)-(a4+a5*tc2))
TNIE.B2 = (b+h3)*((a6+a7*tc3)-(a4+a5*tc2))

PNDE.B1 = (c2+c3*tc1)+h1*(im+a1*t1)
PNDE.A2 = (c2+c3*tc1)-(c4+c5*tc2)+(h1-h2)*((im+a1*t2)+(a4+a5*tc2))
PNDE.B2 = (c6+c7*tc3)-(c4+c5*tc2)-(h2-h3)*((im+a1*t2)+(a4+a5*tc2))

TNDE.B1 = (c2+c3*tc1)+h1*(im+a1*t1+a2+a3*tc1)
TNDE.A2 = (c2+c3*tc1)-(c4+c5*tc2)+(h1-h2)*((im+a1*t1)+(a2+a3*tc1))
TNDE.B2 = (c6+c7*tc3)-(c4+c5*tc2)-(h2-h3)*((im+a1*t3)+(a6+a7*tc3))

CDE.B1 =  (c2+c3*tc1)+h1*m
CDE.A2 =  (c2+c3*tc1)-(c4+c5*tc2)+h1*m-h2*m
CDE.B2 =  (c6+c7*tc3)-(c4+c5*tc2)-h2*m+h3*m

TE.B1 = PNIE.B1 + TNDE.B1
TE.A2 = PNIE.A2 + TNDE.A2
TE.B2 = PNIE.B2 + TNDE.B2

Effect = c(PNIE.B1, PNIE.A2, PNIE.B2, 
           TNIE.B1, TNIE.A2, TNIE.B2, 
           PNDE.B1, PNDE.A2, PNDE.B2,
           TNDE.B1, TNDE.A2, TNDE.B2,
           CDE.B1, CDE.A2, CDE.B2,
           TE.B1, TE.A2, TE.B2)

library(MASS)
Num.Samples = 20000
M.Mean = matrix(c(a1, a2, a3, a4, a5, a6, a7, im), nrow = 1, ncol = 8, byrow = TRUE)
M.Cov  = as.matrix(m.covm[c(2,3,4,5,6,7,8,1),c(2,3,4,5,6,7,8,1)])

set.seed(33612)
M = mvrnorm(n = Num.Samples, mu = M.Mean, Sigma = M.Cov, empirical = FALSE)
M.coeffs = c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'im')
colnames(M) <- M.coeffs

Y.Mean = matrix(c(iy, b, c1, c2, c3, c4, c5, c6, c7, h1, h2, h3), nrow = 1, ncol = 12, byrow = TRUE)
Y.Cov  = as.matrix(y.covm[c(1,9,2,3,4,5,6,7,8,10,11,12),c(1,9,2,3,4,5,6,7,8,10,11,12)])

set.seed(33620)
Y = mvrnorm(n = Num.Samples, mu = Y.Mean, Sigma = Y.Cov, empirical = FALSE)
Y.coeffs = c('iy', 'b', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'h1', 'h2', 'h3')
colnames(Y) <- Y.coeffs

PNIE.B1 = Y[,'b']*(M[,'a2']+M[,'a3']*tc1)
PNIE.A2 = (Y[,'b']+Y[,'h2'])*((M[,'a2']+M[,'a3']*tc1)-(M[,'a4']+M[,'a5']*tc2))
PNIE.B2 = (Y[,'b']+Y[,'h2'])*((M[,'a6']+M[,'a7']*tc3)-(M[,'a4']+M[,'a5']*tc2))

TNIE.B1 = (Y[,'b']+Y[,'h1'])*(M[,'a2']+M[,'a3']*tc1)
TNIE.A2 = (Y[,'b']+Y[,'h1'])*((M[,'a2']+M[,'a3']*tc1)-(M[,'a4']+M[,'a5']*tc2))
TNIE.B2 = (Y[,'b']+Y[,'h3'])*((M[,'a6']+M[,'a7']*tc3)-(M[,'a4']+M[,'a5']*tc2))

PNDE.B1 = (Y[,'c2']+Y[,'c3']*tc1)+Y[,'h1']*(M[,'im']+M[,'a1']*t1)
PNDE.A2 = (Y[,'c2']+Y[,'c3']*tc1)-(Y[,'c4']+Y[,'c5']*tc2)+(Y[,'h1']-Y[,'h2'])*((M[,'im']+M[,'a1']*t2)+(M[,'a4']+M[,'a5']*tc2))
PNDE.B2 = (Y[,'c6']+Y[,'c7']*tc3)-(Y[,'c4']+Y[,'c5']*tc2)-(Y[,'h2']-Y[,'h3'])*((M[,'im']+M[,'a1']*t2)+(M[,'a4']+M[,'a5']*tc2))

TNDE.B1 = (Y[,'c2']+Y[,'c3']*tc1)+Y[,'h1']*(M[,'im']+M[,'a1']*t1+M[,'a2']+M[,'a3']*tc1)
TNDE.A2 = (Y[,'c2']+Y[,'c3']*tc1)-(Y[,'c4']+Y[,'c5']*tc2)+(Y[,'h1']-Y[,'h2'])*((M[,'im']+M[,'a1']*t1)+(M[,'a2']+M[,'a3']*tc1))
TNDE.B2 = (Y[,'c6']+Y[,'c7']*tc3)-(Y[,'c4']+Y[,'c5']*tc2)-(Y[,'h2']-Y[,'h3'])*((M[,'im']+M[,'a1']*t3)+(M[,'a6']+M[,'a7']*tc3))

CDE.B1 =  (Y[,'c2']+Y[,'c3']*tc1)+Y[,'h1']*m
CDE.A2 =  (Y[,'c2']+Y[,'c3']*tc1)-(Y[,'c4']+Y[,'c5']*tc2)+Y[,'h1']*m-Y[,'h2']*m
CDE.B2 =  (Y[,'c6']+Y[,'c7']*tc3)-(Y[,'c4']+Y[,'c5']*tc2)-Y[,'h2']*m+Y[,'h3']*m

TE.B1 = PNIE.B1 + TNDE.B1
TE.A2 = PNIE.A2 + TNDE.A2
TE.B2 = PNIE.B2 + TNDE.B2

MC.estimates = cbind(PNIE.B1, PNIE.A2, PNIE.B2, 
                     TNIE.B1, TNIE.A2, TNIE.B2, 
                     PNDE.B1, PNDE.A2, PNDE.B2, 
                     TNDE.B1, TNDE.A2, TNDE.B2, 
                     CDE.B1, CDE.A2, CDE.B2, 
                     TE.B1, TE.A2, TE.B2)

## 1.Monte Carlo (MC) Method
CI.MC <- matrix(ncol = 2, nrow = ncol(MC.estimates))
rownames(CI.MC) <- colnames(MC.estimates)
colnames(CI.MC) <- c('MC95LCL','MC95UCL')

alpha = 0.05

for (effect in colnames(MC.estimates)) 
{
  CI.MC[effect,] = quantile(MC.estimates[,effect], probs = c(alpha/2, 1-alpha/2))
}

#--------------------------------------------------------------#

## 2.Delta (DLT) Method
coef = c(M.Mean, Y.Mean)
cov  = matrix(0,nrow = nrow(M.Cov)+nrow(Y.Cov), ncol = ncol(M.Cov)+ncol(Y.Cov))
cov[1:nrow(M.Cov), 1:ncol(M.Cov)] <- M.Cov
cov[(1+nrow(M.Cov)):(nrow(M.Cov)+ncol(Y.Cov)), (1+nrow(M.Cov)):(ncol(M.Cov)+ncol(Y.Cov))] <- Y.Cov

library(msm)
z = qnorm(1-alpha/2,mean = 0,sd = 1)
PNIE.B1.se.delta = deltamethod(~ x10*(x2+x3*tc1), coef,cov)
PNIE.A2.se.delta = deltamethod(~ (x10+x19)*((x2+x3*tc1)-(x4+x5*tc2)), coef,cov)
PNIE.B2.se.delta = deltamethod(~ (x10+x19)*((x6+x7*tc3)-(x4+x5*tc2)), coef,cov)

TNIE.B1.se.delta = deltamethod(~ (x10+x18)*(x2+x3*tc1), coef,cov)
TNIE.A2.se.delta = deltamethod(~ (x10+x18)*((x2+x3*tc1)-(x4+x5*tc2)), coef,cov)
TNIE.B2.se.delta = deltamethod(~ (x10+x20)*((x6+x7*tc3)-(x4+x5*tc2)), coef,cov)

PNDE.B1.se.delta = deltamethod(~ (x12+x13*tc1)+x18*(x8+x1*t1), coef,cov)
PNDE.A2.se.delta = deltamethod(~ (x12+x13*tc1)-(x14+x15*tc2)+(x18-x19)*((x8+x1*t2)+(x4+x5*tc2)), coef,cov)
PNDE.B2.se.delta = deltamethod(~ (x16+x17*tc3)-(x14+x15*tc2)-(x19-x20)*((x8+x1*t2)+(x4+x5*tc2)), coef,cov)

TNDE.B1.se.delta = deltamethod(~ (x12+x13*tc1)+x18*(x8+x1*t1+x2+x3*tc1), coef,cov)
TNDE.A2.se.delta = deltamethod(~ (x12+x13*tc1)-(x14+x15*tc2)+(x18-x19)*((x8+x1*t1)+(x2+x3*tc1)), coef,cov)
TNDE.B2.se.delta = deltamethod(~ (x16+x17*tc3)-(x14+x15*tc2)-(x19-x20)*((x8+x1*t3)+(x6+x7*tc3)), coef,cov)

CDE.B1.se.delta = deltamethod(~ (x12+x13*tc1)+x18*m, coef,cov)
CDE.A2.se.delta = deltamethod(~ (x12+x13*tc1)-(x14+x15*tc2)+x18*m-x19*m, coef,cov)
CDE.B2.se.delta = deltamethod(~ (x16+x17*tc3)-(x14+x15*tc2)-x19*m+x20*m, coef,cov)

TE.B1.se.delta = deltamethod(~ x10*(x2+x3*tc1) + (x12+x13*tc1)+x18*(x8+x1*t1+x2+x3*tc1), coef,cov)
TE.A2.se.delta = deltamethod(~ (x10+x19)*((x2+x3*tc1)-(x4+x5*tc2)) + (x12+x13*tc1)-(x14+x15*tc2)+(x18-x19)*((x8+x1*t1)+(x2+x3*tc1)), coef,cov)
TE.B2.se.delta = deltamethod(~ (x10+x19)*((x6+x7*tc3)-(x4+x5*tc2)) + (x16+x17*tc3)-(x14+x15*tc2)-(x19-x20)*((x8+x1*t3)+(x6+x7*tc3)), coef,cov)

SE.delta = c(PNIE.B1.se.delta, PNIE.A2.se.delta, PNIE.B2.se.delta,
             TNIE.B1.se.delta, TNIE.A2.se.delta, TNIE.B2.se.delta, 
             PNDE.B1.se.delta, PNDE.A2.se.delta, PNDE.B2.se.delta,
             TNDE.B1.se.delta, TNDE.A2.se.delta, TNDE.B2.se.delta,
             CDE.B1.se.delta, CDE.A2.se.delta, CDE.B2.se.delta,
             TE.B1.se.delta, TE.A2.se.delta, TE.B2.se.delta)
CI.delta <- matrix(ncol = 2, nrow = ncol(MC.estimates))
rownames(CI.delta) <- colnames(MC.estimates)
colnames(CI.delta) <- c('DLT95LCL','DLT95UCL')
for (effect in 1:ncol(MC.estimates)) 
{
  CI.delta[effect,1] = Effect[effect] - z*SE.delta[effect]
  CI.delta[effect,2] = Effect[effect] + z*SE.delta[effect]
}
M.CI=as.data.frame(cbind(Effect,CI.MC,CI.delta))
}
