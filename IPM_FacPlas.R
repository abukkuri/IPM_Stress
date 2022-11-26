#Code written by Anuraag Bukkuri, expanded from code provided in "Advancing population ecology with integral projection models: a practical guide" - Merow et al 2013

params = data.frame(
  txon = FALSE,
  varG = .05, #.01,.05,1
  varF = .05, #.01,.05,.1
  endtx = 150)

lexp = 1.5 #1.5,2
beta = 1.5 #2
evolv = 1 #1,4
natd = .01*evolv
bex = -.6 #-.6,-.8 or -1

vpar1 = .01 
vpar2 = .2
vpar3 = .5
plotstop = 0
starttx=0
stressG = TRUE
stressF = FALSE
min = 0
max = 1
n = 100
tol = 1e-8
bc = min + c(0:n) * (max - min)/n
y = 0.5 * (bc[1:n] + bc[2:(n+1)])
h = y[2] - y[1]

Nt = rnorm(n,.5,0)
sim_steps = params$endtx-1
v = 0

popsize = vector(,sim_steps)
evoldyn = vector(,sim_steps)
avgtemp = vector(,length(y))
avgplot1 = vector(,sim_steps)
avgvar = vector(,length(y))
avgplot2 = vector(,sim_steps)

for (i in 1:sim_steps) {
  if (params$txon==TRUE) {m=1} #1.5
  else {m=0}
  
  total_cells = sum(Nt)
  Ntnorm = Nt / total_cells
  popsize[i] = total_cells
  evoldyn[i] = v
  for (j in 1:length(y)) {
    avgtemp[j] = y[j] * Ntnorm[j]
    x = m/(exp(lexp*y[j]) + beta*v)
    avgvar[j] = Ntnorm[j]*(vpar1+vpar2*exp(-vpar3/x))  
  }
  avgcell2 = sum(avgvar)
  avgcell = sum(avgtemp)
  avgplot1[i] = avgcell
  avgplot2[i] = avgcell2
  
  if (total_cells>1e6 && params$txon==FALSE){starttx=i
  params$txon=TRUE}
  
  #if (i==47){plotstop=i #47,53
  #break}
  
  if ((total_cells>1e6 && params$txon==TRUE && i>starttx+3)||total_cells<10){plotstop=i
  break}
  
  # if ((total_cells>popsize[i-1] && params$txon==TRUE && i>starttx+2||total_cells<10)){plotstop=i
  # break}
  
  l.x = function(x, params) {exp(lexp*x)}
  
  b.x = function(x, params) {.3*exp(bex*x)} #.4,.5
  
  s.x = function(x, params) {
    return(1-natd-m/(l.x(x,params=params)+beta*v))}
  
  s.x2 = function(x, params) {
    return(1-natd-m/(l.x(x,params=params)+beta*(v+.00001)))}
  
  g.yx = function(xp, x, params) {
    drugd = m/(l.x(x, params=params) + beta*v)
    varG = vpar1+vpar2*exp(-vpar3/drugd)
    if (stressG==FALSE){return (dnorm(xp,mean=x,sd= params$varG))}
    else if (stressG==TRUE){return (dnorm(xp,mean=x,sd= varG))}
  }
  
  g.yx2 = function(xp, x, params) {
    drugd = m/(l.x(x, params=params) + beta*(v+.00001))
    varG = vpar1+vpar2*exp(-vpar3/drugd)
    if (stressG==FALSE){return (dnorm(xp,mean=x,sd= params$varG))}
    else if (stressG==TRUE){return (dnorm(xp,mean=x,sd= varG))}
  }
  
  f.yx = function(xp, x, params) {
    drugd = m/(l.x(x,params=params) + beta*v)
    varF = vpar1+vpar2*exp(-vpar3/drugd)
    if (stressF==FALSE){return (dnorm(xp,mean=x,sd=params$varF))}
    else if (stressF==TRUE){return(dnorm(xp,mean=x,sd= varF))}
  }
  
  f.yx2 = function(xp, x, params) {
    drugd = m/(l.x(x,params=params) + beta*(v+.00001))
    varF = vpar1+vpar2*exp(-vpar3/drugd)
    if (stressF==FALSE){return (dnorm(xp,mean=x,sd=params$varF))}
    else if (stressF==TRUE){return(dnorm(xp,mean=x,sd= varF))}
  }
  
  dsdv <-outer(y, y, function(xp,x,params) {beta*m/(beta*v+exp(lexp*x)) ^ 2}, params=params)
  
  dg1dv <-outer(y,y,g.yx,params=params)
  dg2dv <-outer(y,y,g.yx2,params=params)
  dgdv <- outer(y, y, function(xp,x,params) {(dg2dv-dg1dv)/.00001},params=params)
  
  df1dv <-outer(y,y,f.yx,params=params)
  df2dv <-outer(y,y,f.yx2,params=params)
  dfdv <- outer(y, y, function(xp,x,params) {(df2dv-df1dv)/.00001},params=params)
  
  dkds <-outer(y, y, function(xp,x,params) {g.yx(xp,x,params)},params=params)
  
  dkdg <-outer(y, y, function(xp,x,params) {s.x(x,params)},params=params)
  
  dkdf <-outer(y, y, function(xp,x,params) {b.x(x,params)},params=params)
  
  #print(s.x(0))
  #print(b.x(0))
  
  F = h*outer(y, y, f.yx, params=params)
  B = b.x(y, params=params)
  G = h*outer(y, y, g.yx, params=params)
  S = s.x(y, params=params)
  Q = F
  P = G
  for (j in 1:n) {
    Q[, j] = F[, j] * B[j]
    P[, j] = G[, j] * S[j]
  }
  
  #ADRN Eviction
  for (j in 1:(n / 2)) {
    F[1, j] <- F[1, j] + 1 - sum(F[, j])
    Q[, j] <- F[, j] * B[j]
    G[1, j] <- G[1, j] + 1 - sum(G[, j])
    P[, j] <- G[, j] * S[j]
  }
  
  #MES Eviction
  for (j in (n / 2 + 1):n) {
    F[n, j] <- F[n, j] + 1 - sum(F[, j])
    Q[, j] <- F[, j] * B[j]
    G[n, j] <- G[n, j] + 1 - sum(G[, j])
    P[, j] <- G[, j] * S[j]
  }
  
  K = P + Q
  Nt1 = K %*% Nt
  Nt = Nt1
  
  lam <- abs(eigen(K)$values[1])
  #damp = abs(eigen(K)$values[1])/abs(eigen(K)$values[2]) # damping ratio
  w.eigen <- abs(eigen(K)$vectors[, 1])
  stable.dist <- w.eigen/sum(w.eigen)
  v.eigen <- abs(eigen(t(K))$vectors[, 1])
  repro.val <- v.eigen/v.eigen[1]
  v.dot.w = sum(stable.dist * repro.val) * h
  sens = outer(repro.val, stable.dist)/v.dot.w
  #elas = matrix(as.vector(sens) * as.vector(K)/lam, nrow = n)
  
  conts = dkds * dsdv
  if (stressG==FALSE||m==0){contg = 0}
  else if (stressG==TRUE){contg = dkdg*dgdv}
  if (stressF==FALSE||m==0){contf = 0}
  else if (stressF==TRUE){contf = dkdf*dfdv}
  sensv <- sum(sens*(conts+contg+contf)) * h ^ 2
  v <- v + evolv * sum(sensv) / lam
  
  #print(c(sum(sens*(conts)),i))
  #print(c(sum(sens),i))
  #print(c(sum(conts),i))
  
  #print(c(sum(Nt),i))
  #print(c(avgcell,i))
  #print(c(lam,i))
  
  # jpeg(file = paste("/Users/an7531bu/Downloads/",i))
  # plot(y,Nt,xlab = "State",ylab = "Population Density",type = 'l',main = paste("State distribution: Time ", i),family='serif')
  # dev.off()
}

message('Drug Resistance Level:',v)
message('Average Cell State:',avgcell)
message('Nadir of Population:',min(popsize[-c(1:starttx,plotstop+1:params$endtx)]))
message('Tx Start:',starttx)
message('Time to Progression:',plotstop-starttx)

par(mfrow = c(1, 4))#,oma=c(1,0,2,0))
plot(avgplot1[-c(plotstop+1:params$endtx)],xlab = "Time Steps",ylab='Cell State',type = 'l',main = 'Average Cell State over Time',family='serif')
plot(popsize[-c(plotstop+1:params$endtx)],xlab = "Time Steps",ylab = 'Population Size',main = 'Population Size over Time',type = 'l',family='serif')
abline(h=1e6,col='red',lty=3)
plot(evoldyn[-c(plotstop+1:params$endtx)],xlab = "Time Steps",ylab = 'Trait Value',main = 'Drug Resistance over Time',type = 'l',family='serif')
plot(avgplot2[-c(plotstop+1:params$endtx)],xlab = "Time Steps",ylab = 'Growth Variance',main = 'Average Growth Variance over Time',type = 'l',family='serif',cex.main=1)

# par(mfrow = c(1,2))
# plot(y,s.x(y,params),ylim=c(0.6,1),xlab="State",type="l",ylab="Survival Probability",lwd=12,main = 'Survival w/Correction: Pre-Tx',family='serif') #
# points(y,apply(P,2,sum),col="brown1",lwd=3,cex=.1,pch=19)
# plot(y,b.x(y,params),xlab="State",type="l",ylab="Birth Probability",lwd=12,main = 'Birth w/Correction: Pre-Tx',family='serif')
# points(y,apply(Q,2,sum),col="cyan",lwd=3,cex=.1,pch=19)

# plot(y,Nt,xlab="State",ylab='Pop Density',type='l',main='Population Distribution',family='serif')
# image(y,y,t(elas),xlab = "State (t)",ylab = "State (t+1)",main = "Elasticity",family='serif')
# image(y,y,t(sens),xlab = "State (t)",ylab = "State (t+1)",main = "Sensitivity",family='serif')
#mtext("Low Evolvability", line=0, side=3, outer=TRUE, cex=1.5,family='serif')
#plot(y,stable.dist,xlab = "State",type = "l",main = "Stable state distribution",family='serif')
#plot(y,repro.val,xlab="State",type="l",main="Reproductive values",family='serif')
#image(y,y,t(K), xlab="State (t)",ylab="State (t+1)",col=topo.colors(100), main="IPM matrix",family='serif')
#contour(y,y,t(K), add = TRUE, drawlabels = TRUE)