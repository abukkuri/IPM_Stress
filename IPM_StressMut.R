#Code written by Anuraag Bukkuri, expanded from code provided in "Advancing population ecology with integral projection models: a practical guide" - Merow et al 2013

params = data.frame(
  txon = FALSE,
  varG = .05,
  varF = .05,
  endtx = 150)

lexp = 1.5
beta = 1.5
bex = -.6

dn = .01
vpar1 = 0 
vpar2 = 5
vpar3 = .1
plotstop = 0
starttx=0
stressG = FALSE
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
avgplot = vector(,sim_steps)
avgplot2 = vector(,sim_steps)
avgevolv = vector(,length(y))

for (i in 1:sim_steps) {
  
  if (params$txon==TRUE) {m=1}
  else {m=0}
  
  total_cells = sum(Nt)
  Ntnorm = Nt / total_cells
  popsize[i] = total_cells
  evoldyn[i] = v
  
  kevol = function(v1,params){
    x = m/(exp(lexp*y[j]) + beta*v1)
    avgevolv[j] = Ntnorm[j]*(vpar1+vpar2*exp(-vpar3/x))  
    return (sum(avgevolv))
  }
  
  for (j in 1:length(y)) {avgtemp[j] = y[j] * Ntnorm[j]
  x = m/(exp(lexp*y[j]) + beta*v)
  avgevolv[j] = Ntnorm[j]*(vpar1+vpar2*exp(-vpar3/x))
  }

  avgcell2 = kevol(v,params=params) 
  avgplot2[i] = avgcell2
  
  avgcell = sum(avgtemp)
  avgplot[i] = avgcell
  
  if (total_cells>1e6 && params$txon==FALSE){starttx=i
  params$txon=TRUE}
  
  if ((total_cells>1e6 && params$txon==TRUE && i>starttx)||total_cells<10){plotstop=i
  break} #cancer progression
  
  if ((total_cells>popsize[i-1] && params$txon==TRUE && i>starttx+2||total_cells<10)){plotstop=i
  break} #tx failure
  
  l.x = function(x, params) {exp(lexp*x)}
  
  b.x = function(x, params) {.3*exp(bex*x)}
  
  s.x = function(x, params) {
    return(1-dn*kevol(v,params)-m/(l.x(x,params=params)+beta*v))}
  
  s.x2 = function(x, params) {
    return(1-dn*kevol(v+.0000001,params)-m/(l.x(x,params=params)+beta*(v+.0000001)))}
  
  g.yx = function(xp, x, params) {dnorm(xp,mean=x,sd= params$varG)}
  f.yx = function(xp, x, params) {dnorm(xp,mean=x,sd= params$varF)}
  
  dsdv <-outer(y,y,function(xp,x,params) {(s.x2(x,params)-s.x(x,params))/.0000001}, params=params)
  dkds <-outer(y,y,function(xp,x,params) {g.yx(xp,x,params)},params=params)
  
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
  
  sensv <- sum(sens*(dkds * dsdv)) * h ^ 2
  v <- v + sum(avgevolv) * sum(sensv) / lam
  
  #print(sum(avgevolv))
  
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
plot(avgplot[-c(plotstop+1:params$endtx)],xlab = "Time Steps",ylab='Cell State',type = 'l',main = 'Average Cell State over Time',family='serif')
plot(popsize[-c(plotstop+1:params$endtx)],xlab = "Time Steps",ylab = 'Population Size',main = 'Population Size over Time',type = 'l',family='serif')
abline(h=1e6,col='red',lty=3)
plot(evoldyn[-c(plotstop+1:params$endtx)],xlab = "Time Steps",ylab = 'Trait Value',main = 'Drug Resistance over Time',type = 'l',family='serif')
plot(avgplot2[-c(plotstop+1:params$endtx)],xlab = "Time Steps",ylab='Evolvability',type = 'l',main = 'Average Evolvability over Time',family='serif')

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