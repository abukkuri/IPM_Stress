#Code written by Anuraag Bukkuri, expanded from code provided in "Advancing population ecology with integral projection models: a practical guide" - Merow et al 2013

params = data.frame(
  txon = FALSE,
  endtx = 150)

varG = .05
varF = .05
lexp = 1.5
beta = 1.5 
evolv = 1
natd = .01*evolv
bex = -.6

vpar1 = .01
vpar2 = .3
vpar3 = .4
plotstop = 0
starttx=0
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

for (i in 1:sim_steps) {

  total_cells = sum(Nt)
  Ntnorm = Nt / total_cells
  popsize[i] = total_cells
  evoldyn[i] = v
  for (j in 1:length(y)) {avgtemp[j] = y[j] * Ntnorm[j]}
  avgcell = sum(avgtemp)
  avgplot[i] = avgcell
  
  if (total_cells>1e6 && params$txon==FALSE){starttx=i
  params$txon=TRUE}
  
  if ((total_cells>1e6 && params$txon==TRUE && i>starttx+3)||total_cells<10){plotstop=i
  break}
  
  if ((total_cells>popsize[i-1] && params$txon==TRUE && i>starttx||total_cells<10)){plotstop=i
  break} #time to failure
  
  if (params$txon==TRUE) {m=1}
  else {m=0}
  
  l.x = function(x, params) {exp(lexp*x)}
  
  b.x = function(x, params) {.3*exp(bex*x)}
  
  s.x = function(x, params) {
    return(1-natd-m/(l.x(x,params=params)+beta*v))}
  
  g.yx = function(xp, x, params) {dnorm(xp,mean=x,sd= varG)}
  
  f.yx = function(xp, x, params) {dnorm(xp,mean=x,sd= varF)}
  
  dsdv <-outer(y, y, function(xp,x,params) {beta*m/(beta*v+exp(lexp*x))^2}, params=params)
  
  dkds <-outer(y, y, function(xp,x,params) {g.yx(xp,x,params)},params=params)
  
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
  
  #Low state eviction
  for (j in 1:(n / 2)) {
    F[1, j] <- F[1, j] + 1 - sum(F[, j])
    Q[, j] <- F[, j] * B[j]
    G[1, j] <- G[1, j] + 1 - sum(G[, j])
    P[, j] <- G[, j] * S[j]
  }

  #High state eviction
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
  w.eigen <- abs(eigen(K)$vectors[, 1])
  stable.dist <- w.eigen/sum(w.eigen)
  v.eigen <- abs(eigen(t(K))$vectors[, 1])
  repro.val <- v.eigen/v.eigen[1]
  v.dot.w = sum(stable.dist * repro.val) * h
  sens = outer(repro.val, stable.dist)/v.dot.w
  
  sensv <- sum(sens*(dkds * dsdv)) * h ^ 2
  v <- v + evolv * sum(sensv) / lam
}

message('Drug Resistance Level:',v)
message('Average Cell State:',avgcell)
message('Nadir of Population:',min(popsize[-c(1:starttx,plotstop+1:params$endtx)]))
message('Tx Start:',starttx)
message('Time to Progression:',plotstop-starttx)

par(mfrow = c(1, 3))
plot(avgplot[-c(plotstop+1:params$endtx)],xlab = "Time Steps",ylab='Cell State',type = 'l',main = 'Average Cell State over Time',family='serif')
plot(popsize[-c(plotstop+1:params$endtx)],xlab = "Time Steps",ylab = 'Population Size',main = 'Population Size over Time',type = 'l',family='serif')
abline(h=1e6,col='red',lty=3)
plot(evoldyn[-c(plotstop+1:params$endtx)],xlab = "Time Steps",ylab = 'Trait Value',main = 'Drug Resistance over Time',type = 'l',family='serif')
