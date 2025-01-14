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
library(ExtDist)

file = getwd()
source(paste0(file,"/funcs.r"))
source(paste0(file,"/funcs2.R"))


Comparison_of_elliptical_test = function(rp,setting_distribution,k,setting_sigma,B,num,m = 400){
  d = m * rp #r = 0.5,1,1.5
  BB = 500  #repeat times for Step 2 in Algo 1
  
  
  
  Ypm_1 = matrix(0,B,2)
  T_all = sigma2 = quant_stat = rep(0,B)
  switch (setting_sigma,
          'Sig1' = {
            lambda_pop = diag( rep(1,d) )
            lambda_pop[1:5,1:5] = lambda_pop[1:5,1:5]*5
            Q = qr(matrix(rnorm(d*d, mean = 0, sd = 1),d,d))
            Q = qr.Q(Q)
            Sigma = Q%*%lambda_pop%*%t(Q)
          },
          'Sig2' = {
            Sigma = matrix(0,d,d)
            
            for (i in 1:d) {
              for (j in 1:d) {
                Sigma[i,j] = 0.1^{abs(i-j)}
              }
            }
          },
          'Sig3' = {
            k1= 1/4
            lambda_pop = diag( seq(1,d)^(-k1) )
            Q = qr(matrix(rnorm(d*d, mean = 0, sd = 1),d,d))
            Q = qr.Q(Q)
            Sigma = Q%*%lambda_pop%*%t(Q)
          },
          'Sig4' = {
            Sigma = diag( rep(1,d) )
          }
  )
  
  
  for (ii in 1:B){
    switch (setting_distribution,
            'dist1' = {
              s = 1
              Y = matrix(rLaplace(m*d, 0, s),m,d)/sqrt(2/s^2)
              Z = matrix(rnorm(m*d),m,d)
              X = Re(((1-k)^0.5 * Z + k^{0.5} * Y) %*% sqrtm(Sigma))
            },
            'dist2' = {
              a = 2
              b = 1.5
              Y = (matrix(rbeta(m*d, a, b),m,d)- (a/(a+b)))/sqrt(a*b/(a+b)^2/(a+b+1))
              Z = matrix(rnorm(m*d),m,d)
              X = Re(((1-k)^0.5 * Z + k^{0.5} * Y) %*% sqrtm(Sigma))
            }
    )
    
    
    X = Re(X)
    
    #normality test 
    
    ## L = 1
    L = 1  #estimate p-values L times
    # ATX = adapt.thres.cov(X)
    # eigenv <- eigen(ATX, symmetric = TRUE)
    # e.vec <- as.matrix(eigenv$vectors)
    # sqrS <- e.vec %*% diag(sqrt(eigenv$values^(-1)), ncol = d) %*% t(e.vec)
    # 
    # XX = X %*% sqrS
    # XD = X %*% diag(diag(sqrS))
    temp = getp(X,L=L,BB=BB)
    # temp2 = getp(XX,L=L,BB=BB)
    # temp3 = getp(XD,L=L,BB=BB)
    Ypm_1[ii,] = temp$Yp
    # YYpm_1[ii,] = temp2$Yp  #with sigma_x^{-1/2} plug in
    # YDpm_1[ii,] = temp3$Yp   #with D_x^{-1/2} plug in
    # Opm_1[ii,] = temp$Op
    
    ## L = 5
    # start.time <- Sys.time()
    # L=5
    # ATX = adapt.thres.cov(X)
    # eigenv <- eigen(ATX, symmetric = TRUE)
    # e.vec <- as.matrix(eigenv$vectors)
    # sqrS <- e.vec %*% diag(sqrt(eigenv$values^(-1)), ncol = d) %*% t(e.vec)
    # 
    # XX = X %*% sqrS
    # XD = X %*% diag(diag(sqrS))
    # temp = getp(X,L=L,BB=BB)
    # temp2 = getp(XX,L=L,BB=BB)
    # temp3 = getp(XD,L=L,BB=BB)
    # Ypm_5[ii,] = temp$Yp
    # YYpm_5[ii,] = temp2$Yp  #with sigma_x^{-1/2} plug in
    # YDpm_5[ii,] = temp3$Yp   #with D_x^{-1/2} plug in
    # Opm_5[ii,] = temp$Op

    
    
    # elliptical test
    n1 = m/2
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
    
    
    T_all[ii] = sqrt(d*n1) *(T1-T2)
    
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
    sigma2[ii] = var1 + var2
    
    quant_stat[ii] = pnorm(-abs(T_all[ii]/sqrt(sigma2[ii]))) *2
  }
  write.csv(Ypm_1,file = paste('norm',k,m,d,setting_distribution,setting_sigma,num,".csv", sep="_"), row.names = FALSE)
  write.csv(cbind(T_all,sigma2,quant_stat),paste("elliptical",k,m,d,setting_distribution,setting_sigma,num,".csv", sep="_"), row.names = FALSE)
}



