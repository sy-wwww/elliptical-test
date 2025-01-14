library(mvtnorm)
library(mnormt)
library(Rlab)
library(gTests)
library(expm)
library(MASS)
library(normwhn.test)
library(psych)
library(tilting)
library(Matrix)
library(foreach)

library(wordspace)


file = getwd()
source(paste0(file,"/funcs.r"))
source(paste0(file,"/funcs2.R"))


normality_test = function(X){
  m = dim(X)[1]
  d = dim(X)[2]
  
  BB = 500  #repeat times for Step 2 in Algo 1
  L = 1  #estimate p-values L times
  
  
  # PseudoGaussian_re = matrix(0,B,1)
  Ypm = Opm = YYpm = YDpm = matrix(0,1,2)
  
  ATX = adapt.thres.cov(X)
  eigenv <- eigen(ATX, symmetric = TRUE)
  e.vec <- as.matrix(eigenv$vectors)
  sqrS <- e.vec %*% diag(sqrt(eigenv$values^(-1)), ncol = d) %*% t(e.vec)
  
  XX = X %*% sqrS
  XD = X %*% diag(diag(sqrS))
  temp = getp(X,L=L,BB=BB)
  temp2 = getp(XX,L=L,BB=BB)
  temp3 = getp(XD,L=L,BB=BB)
  Ypm[1,] = temp$Yp
  YYpm[1,] = temp2$Yp  #with sigma_x^{-1/2} plug in
  YDpm[1,] = temp3$Yp   #with D_x^{-1/2} plug in
  Opm[1,] = temp$Op
  aY = getpow(Ypm)  # our test
  aO = getpow(Opm)  # eFR
  aYY = getpow(YYpm)  # our test with sigma_x^{-1/2} plug in
  aYD = getpow(YDpm)
  return(mean(Ypm))
}


library(psych)
library(tilting)
elliptical_test = function(X){
  m = length(X[,1])
  d = length(X[1,])
  n1 = floor(m/2)
  n2 = m-n1
  
  
  X1 = X[1:n1,]
  X2 = X[(n1+1):m,]
  
  
  # data1 
  b1 = colMeans(X1^4)
  b2 = 3 * (colMeans(X1^2)^2)
  T1 = mean(b1/b2) + 2/n1
  
  # data2
  sigHat = (1/n2)*t(X2)%*% X2
  omega_n = tr(sigHat)^2
  gamma_n = tr(sigHat%*%sigHat) - omega_n/n2
  l2norm = col.norm(t(X2))^2
  nu_n = sum((l2norm-mean(l2norm))^2)/(n2-1)
  T2 = (nu_n + omega_n)/(omega_n + 2*gamma_n)
  
  
  T_all = sqrt(d*n1) *(T1-T2)
  
  # all data
  sigHat = (1/m)*t(X)%*% X
  b1 = tr(sigHat)
  b2 = tr(sigHat%*%sigHat)
  b3 = tr(sigHat%*%sigHat%*%sigHat)
  b4 = tr(sigHat%*%sigHat%*%sigHat%*%sigHat)
  A = diag(1/sqrt(diag(sigHat)))%*%matrix(sigHat,d,d)%*%diag(1/sqrt(diag(sigHat)))
  
  var2 = d *(16*b4+8* (b2 - b1^2/m)^2)/ (b1^2+2*(b2 - b1^2/m))^2
  
  l2norm = col.norm(t(X))^2
  nu_n = sum((l2norm-mean(l2norm))^2)/(m-1)
  a1h = (nu_n +b1^2)/(b1^2 + 2*(b2 - b1^2/m))
  a2h = mean(col.norm(t(X))^6)/(b1^3+6*(b2 - b1^2/m)*b1+ 8 *b3)
  a3h = mean(col.norm(t(X))^8)/(b1^4+12*(b2 - b1^2/m)*b1^2+32*b1*b3 + 12 *(b2 - b1^2/m)^2 + 48 *b4)
  c_t = a3h + a1h - 2*a2h
  
  if (c_t > log(d) * d^(-3/4)){
    c_t = log(d) * d^(-3/4)
  } else if (c_t < -log(d) * d^(-3/4)){
    c_t = -log(d) * d^(-3/4)
  }
  
  
  var1 = 8/3*mean(A^4)*d * a3h + 8*c_t *mean(A^2)*d 
  sigma2 = var1 + var2
  
  quant_stat = pnorm(-abs(T_all/sqrt(sigma2))) *2
  return(c(T_all,sigma2,quant_stat))
}


