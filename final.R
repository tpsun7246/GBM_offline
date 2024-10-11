library(dplyr)
library(LSMRealOptions)
library(tidyverse)
library(dplyr)
library(quantmod)
library(beepr)

rm(list = ls())
#################simulation
SDE.data = function(t, n, mu, sigma, Y0){
  dt = 1/n
  
  t.seq = seq(0,t,by=dt)[-1]
  Wt.seq = sqrt(dt)*cumsum(rnorm(n*t,0,1))
  Y = Y0 + mu*t.seq + sigma*Wt.seq
  Y = c(Y0,Y)
  
  index = seq(1,n*t+1,by=n)
  O = Y[index[-length(index)]]
  C = Y[index[-1]]
  U = c()
  L = c()
  for (i in 1:t){
    dd = Y[index[i]:index[i+1]]
    U[i] = max(dd)
    L[i] = min(dd)
  }
  data = matrix(c(O, C, U, L), ncol = 4, byrow = FALSE)
  colnames(data) = c("O","C","U","L")
  
  return(data)
}










l=function(mu,sigma,data){
  f=c()
  
  
  o = data[,"O"]
  c = data[,"C"]
  u = data[,"U"]
  l = data[,"L"]
  k = seq(-10,10,by=1)
  sigma=exp(sigma)
  s=sigma^2
  
  
  
  for (i in 1:length(u)) {
    
    
    g1 = (4*k*(k+1))/(sqrt(2*pi)*sigma^3)*(1-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    g2 = (4*k^2)/(sqrt(2*pi)*sigma^3)*(1-((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    logf1 = -((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    logf2 = -((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    h1=exp(logf1)*g1
    h2=exp(logf2)*g2
    
    
    f[i]=log(sum(h1)-sum(h2))
    
  }
  f=sum(f)
  return(f)
} # likelihood



partial_ell_sigma=function(mu, sigma, data){
  o = data[,"O"]
  c = data[,"C"]
  u = data[,"U"]
  l = data[,"L"]
  k=seq(-10,10,1)
  sigma=exp(sigma)
  s=sigma^2
  li=c()
  
  for (i in 1:length(u)) {
    g1= (4*k*(k+1))/(sqrt(2*pi)*sigma^3)*(1-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    g2 = (4*k^2)/(sqrt(2*pi)*sigma^3)*(1-((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    logh1=-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    h1=exp(logh1)
    logh2 = -((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    h2=exp(logh2)
    z=sum(g1 * h1 - g2* h2)
    d_g1=(-3/(2*sigma^2))*g1+(4*k*(k+1))/sqrt(2*pi)*(c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    
    d_h1=-(h1*logh1)/sigma^2
    
    d_g2=(-3/(2*sigma^2))*g2+(4*k^2)/sqrt(2*pi)*(c[i]-o[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    
    d_h2=-(h2*logh2)/sigma^2
    
    li[i]=(sum(d_g1*h1+d_h1*g1)-sum(d_g2*h2+d_h2*g2))/z
    
  }
  
  li=sum(li)
  return(li)
}







ddsigma= function(mu,sigma,data){
  o = data[,"O"]
  c = data[,"C"]
  u = data[,"U"]
  l = data[,"L"]
  df=c()
  h=c()
  dh=c()
  k = seq(-10,10,by=1)
  sigma=exp(sigma)
  s=sigma^2
  
  
  
  for (i in 1:length(u)) {
    
    
    
    g1 = (4*k*(k+1))/(sqrt(2*pi)*sigma^3)*(1-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    g2 = (4*k^2)/(sqrt(2*pi)*sigma^3)*(1-((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    logf1 = -((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    logf2 = -((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    f1=exp(logf1)
    f2=exp(logf2)
    dg1 = (-3/(2*sigma^2))*g1+(4*k*(k+1))/sqrt(2*pi)*(c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    dg2 = (-3/(2*sigma^2))*g2+(4*k^2)/sqrt(2*pi)*(c[i]-o[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    df1 = -(exp(logf1)*logf1)/sigma^2
    df2 = -(exp(logf2)*logf2)/sigma^2
    
    a=sum(dg1*exp(logf1)+g1*df1)
    b=sum(dg2*exp(logf2)+g2*df2)
    cc=sum(g1*f1)-sum(g2*f2)
    da=sum(3*(s*(dg1*f1+g1*df1)+g1*f1)/(2*s^2)+(4*k*(k+1))/(sqrt(2*pi))*((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)*(-(7*s^(-9/2)*f1)/2+s^(-7/2)*df1)-(s*(df1*logf1*g1+df1*g1+f1*dg1*logf1)-logf1*g1*f1)/(s^2))
    db=sum(3*(s*(dg2*f2+g2*df2)+g2*f2)/(2*s^2)+(4*k^(2))/(sqrt(2*pi))*((c[i]-o[i]-2*k*(u[i]-l[i]))^2)*(-(7*s^(-9/2)*f2)/2+s^(-7/2)*df2)-(s*(df2*logf2*g2+df2*g2+f2*dg2*logf2)-logf2*g2*f2)/(s^2))
    h[i]=a-b
    dh[i]=da-db
    df[i]=(dh[i]*cc-h[i]^2)/(cc^2)
  }
  df=sum(df)
  #ifelse(df == -Inf,-10^(30),df)
  return(df)
  
  
}# second order partial



f=function(data,mu,sig,tol=10^(-6),max_iter=10000,j){
  
  i=0
  repeat {
    
    
    dsigma_sq1=partial_ell_sigma(mu = mu, sigma = sig,data = data)
    ddsigma_sq1=ddsigma(mu = mu, sigma = sig,data = data)
    
    dsig1=2*dsigma_sq1*exp(sig)
    ddsig1=4*ddsigma_sq1*exp(2*sig)+2*dsigma_sq1
    
    ds=dsig1*exp(sig)
    dds=ddsig1*exp(2*sig)+dsig1*exp(sig)
    
    value1 = ds
    
    value2 = dds
    
    if(is.infinite(value1)|is.infinite(value2)){
      
      break
    }
    
    
    sig_new <- sig -  sum(value1)/sum(value2)
    
    #ifelse(sig_new<0,runif(1,0.02,0.04),sig_new)
    
    i <- i + 1
    if(is.nan(sig_new)){
      break
    }
    
    if(exp(sig_new)<0){
      break
    }
    
    if (abs(sig_new - sig) < tol|| i >= max_iter ) {
      
      break
    }
    
    
    if(exp(sig)>1){
      break
    }
    
    
    sig <- sig_new
    print(c(j,exp(sig),i))
  }
  return(sig)
}



point.est = function(data,sig0,j){
  O = data[,"O"]
  C = data[,"C"]
  U = data[,"U"]
  L = data[,"L"]
  mu.hat = mean(C-O)
  sigma.hat =f(data=data,mu=mu.hat,sig = sig0,j=j) 
  if(is.nan(sigma.hat)){
    browser()
  }
  
  return(list(mu.hat = mu.hat, sigma.hat = sigma.hat))
}

change.point.detect.u = function(mu1,mu2,sigma1,sigma2,tau,n,Y0,n.sim,n.check,sig0){
  est.matrix = matrix(NA, ncol = 5, nrow = n.sim)
  colnames(est.matrix) = c("tau","mu1","mu2","sigma21","sigma22")
  true.value = c(tau,mu1,mu2,sigma1^2,sigma2^2)
  data_set=list()
  sig.ll=list()
  for (j in 1:n.sim){
    
    data1 = SDE.data(tau,2500000,mu1,sigma1,Y0=Y0)
    data2 = SDE.data(n-tau,2500000,mu2,sigma2,Y0=data1[tau,2])
    data = rbind(data1,data2)
    data_set[[j]]=data
    sig.ll1=matrix(0,nrow =2,ncol=length(sig0) )
    sig.ll2=matrix(0,nrow =2,ncol=length(sig0) )
    
    
    ll=c()
    for (i in (tau-5):(tau+5)){
      
      
      
      d1 = data[1:i,]
      d2 = data[(i+1):n,]
      a=unlist(lapply(sig0, function(x) point.est(sig0 = x, data = d1,j=j)))
      b=unlist(lapply(sig0, function(x) point.est(sig0 = x, data = d2,j=j)))
      est1.i = matrix(a,nrow = length(sig0),ncol = 2,byrow = T)
      
      est2.i = matrix(b,nrow = length(sig0),ncol = 2,byrow = T)
      
      l1 <- apply(est1.i,1, function(x) l(data = d1, mu = x[1], sig = x[2]))
      
      index1 <- which.max(l1)
      l1 <- max(l1,na.rm = TRUE)
      
      l2 = apply(est2.i, 1, function(x) l(data = d2, mu = x[1], sig = x[2])) 
      
      index2=which.max(l2)
      l2=max(l2,na.rm = TRUE)
      ll[i-(n.check-1)]=l1+l2
      
      
      
      
      
      
    }
    
    tau.hat = which.max(ll)+n.check-1
    
    
    a=unlist(lapply(sig0, function(x) point.est(sig0 = x, data = data[1:tau.hat,],j=j)))
    b=unlist(lapply(sig0, function(x) point.est(sig0 = x, data = data[(tau.hat+1):n,],j=j)))
    est1 = matrix(a,nrow = length(sig0),ncol = 2,byrow = T)
    
    est2 = matrix(b,nrow = length(sig0),ncol = 2,byrow = T)
    
    l1 <- apply(est1,1, function(x) l(data = data[1:tau.hat,], mu = x[1], sig = x[2]))
    
    sig.ll1[1,]=exp(est1[,2])
    sig.ll1[2,]=l1
    index1 <- which.max(l1)
    sig.ll2[1,]=exp(est2[,2])
    sig.ll2[2,]=l2
    l2 = apply(est2, 1, function(x) l(data = data[(tau.hat+1):n,], mu = x[1], sig = x[2])) 
    index2=which.max(l2)
    
    est1.mu=est1[index1,1]
    est1.sig=est1[index1,2]
    est2.mu=est2[index2,1]
    est2.sig=est2[index2,2]
    
    est.matrix[j,] = c(tau.hat, est1.mu, est2.mu,  exp(est1.sig)^2, exp(est2.sig)^2)
    sig.ll[[j]]=list(sig.ll1=sig.ll1,sig.ll2=sig.ll2)
    
  }
  
  RMSE = lapply(1:ncol(est.matrix), function(i) sqrt(mean((est.matrix[,i]-true.value[i])^2))) %>% as.numeric()
  RE=RMSE/true.value
  est = as.numeric(apply(est.matrix, 2, mean))
  
  value = data.frame(    
    true.value = true.value,
    estimate = est,
    RMSE = RMSE,
    RE=RE,
    #con_interval=con_interval,
    row.names = c("tau","mu1","mu2","sigma21","sigma22")
  )
  return(list( sig.ll=sig.ll,est.matrix = est.matrix,value = value,data_set=data_set))
}
################################################################################
## main code

#######################
## point estimate





#######################sigma
##  changepoint detection

## setting

n = 250
mu1 =  0.0008 ## mu before change point
mu2 = 0.0008  ## mu after change point
sigma1 = 0.013 ## sigma before change point
sigma2 = 0.028 ## sigma after change point



tau = 25 ## true change point
Y0 = log(100) ## initial stock price
n.sim = 1000 ## simulation times
n.check = tau-5 ## check the likelihood from #n.check to #(data.length-n.check)

sig0=log(c(0.01,0.02,0.03))
p1=Sys.time()
set.seed(112)
result = change.point.detect.u(mu1,mu2,sigma1,sigma2,tau,n,Y0,n.sim,n.check,sig0)
options(scipen = 999)
result
p2=Sys.time()
p2-p1







n = 250
mu1 =  0.0008 ## mu before change point
mu2 = 0.0008  ## mu after change point
sigma1 = 0.013 ## sigma before change point
sigma2 = 0.028 ## sigma after change point



tau = 50 ## true change point
Y0 = log(100) ## initial stock price
n.sim = 1000 ## simulation times
n.check = tau-5 ## check the likelihood from #n.check to #(data.length-n.check)

sig0=log(c(0.01,0.02,0.03))
p1=Sys.time()
set.seed(112)
result = change.point.detect.u(mu1,mu2,sigma1,sigma2,tau,n,Y0,n.sim,n.check,sig0)
options(scipen = 999)
result
p2=Sys.time()
p2-p1






#######################################################

n = 250
mu1 =  0.0008 ## mu before change point
mu2 = 0.0008  ## mu after change point
sigma1 = 0.013 ## sigma before change point
sigma2 = 0.028 ## sigma after change point



tau = 83 ## true change point
Y0 = log(100) ## initial stock price
n.sim = 1000 ## simulation times
n.check = tau-5 ## check the likelihood from #n.check to #(data.length-n.check)

sig0=log(c(0.01,0.02,0.03))
p1=Sys.time()
set.seed(112)
result = change.point.detect.u(mu1,mu2,sigma1,sigma2,tau,n,Y0,n.sim,n.check,sig0)
options(scipen = 999)
result
p2=Sys.time()
p2-p1






#######################################################


n = 250
mu1 =  0.0008 ## mu before change point
mu2 = 0.0008  ## mu after change point
sigma1 = 0.013 ## sigma before change point
sigma2 = 0.028 ## sigma after change point



tau = 125 ## true change point
Y0 = log(100) ## initial stock price
n.sim = 1000 ## simulation times
n.check = tau-5 ## check the likelihood from #n.check to #(data.length-n.check)

sig0=log(c(0.01,0.02,0.03))
p1=Sys.time()
set.seed(112)
result = change.point.detect.u(mu1,mu2,sigma1,sigma2,tau,n,Y0,n.sim,n.check,sig0)
options(scipen = 999)
result
p2=Sys.time()
p2-p1






#############################

#######################mu
##  changepoint detection

## setting

n = 250
mu1 =  0.0008 ## mu before change point
mu2 = 0.004  ## mu after change point
sigma1 = 0.013 ## sigma before change point
sigma2 = 0.013 ## sigma after change point



tau = 25 ## true change point
Y0 = log(100) ## initial stock price
n.sim = 1000 ## simulation times
n.check = tau-5 ## check the likelihood from #n.check to #(data.length-n.check)

sig0=log(c(0.01,0.02,0.03))
p1=Sys.time()
set.seed(112)
result = change.point.detect.u(mu1,mu2,sigma1,sigma2,tau,n,Y0,n.sim,n.check,sig0)
options(scipen = 999)
result
p2=Sys.time()
p2-p1







n = 250
mu1 =  0.0008 ## mu before change point
mu2 = 0.004  ## mu after change point
sigma1 = 0.013 ## sigma before change point
sigma2 = 0.013 ## sigma after change point



tau = 50 ## true change point
Y0 = log(100) ## initial stock price
n.sim = 1000 ## simulation times
n.check = tau-5 ## check the likelihood from #n.check to #(data.length-n.check)

sig0=log(c(0.01,0.02,0.03))
p1=Sys.time()
set.seed(112)
result = change.point.detect.u(mu1,mu2,sigma1,sigma2,tau,n,Y0,n.sim,n.check,sig0)
options(scipen = 999)
result
p2=Sys.time()
p2-p1






#######################################################

n = 250
mu1 =  0.0008 ## mu before change point
mu2 = 0.004  ## mu after change point
sigma1 = 0.013 ## sigma before change point
sigma2 = 0.013 ## sigma after change point



tau = 83 ## true change point
Y0 = log(100) ## initial stock price
n.sim = 1000 ## simulation times
n.check = tau-5 ## check the likelihood from #n.check to #(data.length-n.check)

sig0=log(c(0.01,0.02,0.03))
p1=Sys.time()
set.seed(112)
result = change.point.detect.u(mu1,mu2,sigma1,sigma2,tau,n,Y0,n.sim,n.check,sig0)
options(scipen = 999)
result
p2=Sys.time()
p2-p1






#######################################################mu


n = 250
mu1 =  0.0008 ## mu before change point
mu2 = 0.004  ## mu after change point
sigma1 = 0.013 ## sigma before change point
sigma2 = 0.013 ## sigma after change point



tau = 125 ## true change point
Y0 = log(100) ## initial stock price
n.sim = 1000 ## simulation times
n.check = tau-5 ## check the likelihood from #n.check to #(data.length-n.check)

sig0=log(c(0.01,0.02,0.03))
p1=Sys.time()
set.seed(112)
result = change.point.detect.u(mu1,mu2,sigma1,sigma2,tau,n,Y0,n.sim,n.check,sig0)
options(scipen = 999)
result
p2=Sys.time()
p2-p1





#################################epmircial


rm(list = ls())






#######################################################################################################
l=function(mu,sigma,data){
  f=c()
  o = data[,1]%>%as.numeric()
  u = data[,2]%>%as.numeric()
  l = data[,3]%>%as.numeric()
  c = data[,4]%>%as.numeric()
  k = seq(-10,10,by=1)
  sigma=exp(sigma)
  s=sigma^2
  
  
  
  for (i in 1:length(u)) {
    
    
    g1 = (4*k*(k+1))/(sqrt(2*pi)*sigma^3)*(1-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    g2 = (4*k^2)/(sqrt(2*pi)*sigma^3)*(1-((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    logf1 = -((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    logf2 = -((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    h1=exp(logf1)*g1
    h2=exp(logf2)*g2
    
    
    f[i]=log(sum(h1)-sum(h2))
    
  }
  f=sum(f)
  return(f)
} # likelihood



dsigma = function(mu, sigma, data){
  o = data[,1]%>%as.numeric()
  u = data[,2]%>%as.numeric()
  l = data[,3]%>%as.numeric()
  c = data[,4]%>%as.numeric()
  f=c()
  k = seq(-10,10,by=1)
  sigma=exp(sigma)
  s=sigma^2
  
  
  
  
  for (i in 1:length(u)) {
    
    
    g1 = (4*k*(k+1))/(sqrt(2*pi)*sigma^3)*(1-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    g2 = (4*k^2)/(sqrt(2*pi)*sigma^3)*(1-((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    logf1 = -((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    logf2 = -((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    
    dg1 = (-3/(2*sigma^2))*g1+(4*k*(k+1))/sqrt(2*pi)*(c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    dg2 = (-3/(2*sigma^2))*g2+(4*k^2)/sqrt(2*pi)*(c[i]-o[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    df1 = -(exp(logf1)*logf1)/sigma^2
    df2 = -(exp(logf2)*logf2)/sigma^2
    f[i] = (sum(dg1*exp(logf1)+g1*df1)-sum(dg2*exp(logf2)+g2*df2))/(sum(g1*exp(logf1))-sum(g2*exp(logf2)))
    
  }
  
  
  
  f=sum(f)
  #ifelse(df == -Inf,-10^(30),df)
  
  return(f)
}# first order partial



ddsigma= function(mu,sigma,data){
  o = data[,1]%>%as.numeric()
  u = data[,2]%>%as.numeric()
  l = data[,3]%>%as.numeric()
  c = data[,4]%>%as.numeric()
  df=c()
  h=c()
  dh=c()
  k = seq(-10,10,by=1)
  sigma=exp(sigma)
  s=sigma^2
  
  
  
  for (i in 1:length(u)) {
    
    
    
    g1 = (4*k*(k+1))/(sqrt(2*pi)*sigma^3)*(1-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    g2 = (4*k^2)/(sqrt(2*pi)*sigma^3)*(1-((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/sigma^2)
    logf1 = -((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    logf2 = -((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/(2*sigma^2)-mu^2/(2*sigma^2)+(mu*(c[i]-o[i]))/sigma^2
    f1=exp(logf1)
    f2=exp(logf2)
    dg1 = (-3/(2*sigma^2))*g1+(4*k*(k+1))/sqrt(2*pi)*(c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    dg2 = (-3/(2*sigma^2))*g2+(4*k^2)/sqrt(2*pi)*(c[i]-o[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    df1 = -(exp(logf1)*logf1)/sigma^2
    df2 = -(exp(logf2)*logf2)/sigma^2
    
    a=sum(dg1*exp(logf1)+g1*df1)
    b=sum(dg2*exp(logf2)+g2*df2)
    cc=sum(g1*f1)-sum(g2*f2)
    da=sum(3*(s*(dg1*f1+g1*df1)+g1*f1)/(2*s^2)+(4*k*(k+1))/(sqrt(2*pi))*((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)*(-(7*s^(-9/2)*f1)/2+s^(-7/2)*df1)-(s*(df1*logf1*g1+df1*g1+f1*dg1*logf1)-logf1*g1*f1)/(s^2))
    db=sum(3*(s*(dg2*f2+g2*df2)+g2*f2)/(2*s^2)+(4*k^(2))/(sqrt(2*pi))*((c[i]-o[i]-2*k*(u[i]-l[i]))^2)*(-(7*s^(-9/2)*f2)/2+s^(-7/2)*df2)-(s*(df2*logf2*g2+df2*g2+f2*dg2*logf2)-logf2*g2*f2)/(s^2))
    h[i]=a-b
    dh[i]=da-db
    df[i]=(dh[i]*cc-h[i]^2)/(cc^2)
  }
  df=sum(df)
  #ifelse(df == -Inf,-10^(30),df)
  return(df)
  
  
}# second order partial



f=function(data,mu,sig,tol=10^(-6),max_iter=10000){
  
  i=0
  repeat {
    
    dsigma_sq=dsigma(mu = mu, sigma = sig,data = data)
    ddsigma_sq=ddsigma(mu = mu, sigma = sig,data = data)
    
    dsig=2*dsigma_sq*exp(sig)
    ddsig=4*ddsigma_sq*exp(2*sig)+2*dsigma_sq*exp(sig)
    
    ds=2*dsigma_sq*exp(2*sig)
    dds=ddsig*exp(2*sig)+dsig*exp(sig)
    
    value1 = ds
    
    value2 = dds
    
    if(is.infinite(value1)|is.infinite(value2)){
      
      break
    }
    
    
    sig_new <- sig -  sum(value1)/sum(value2)
    
    #ifelse(sig_new<0,runif(1,0.02,0.04),sig_new)
    
    i <- i + 1
    if(is.nan(sig_new)){
      break
    }
    
    if(exp(sig_new)<0){
      break
    }
    
    if (abs(sig_new - sig) <= tol|| i >= max_iter ) {
      
      break
    }
    
    
    
    
    sig <- sig_new
    print(c(sig,i))
  }
  return(sig)
}



point.est = function(data,sig0){
  o = data[,1]%>%as.numeric()
  u = data[,2]%>%as.numeric()
  l = data[,3]%>%as.numeric()
  c = data[,4]%>%as.numeric()
  mu.hat = mean(c-o)
  sigma.hat =f(data=data,mu=mu.hat,sig = sig0) 
  if(is.nan(sigma.hat)){
    browser()
  }
  
  return(list(mu.hat = mu.hat, sigma.hat = sigma.hat))
}


#########################################################
partial_ell_mu1_mu1 <- function(mu1, sigma1, data) {
  o = data[,1]%>%as.numeric()
  u = data[,2]%>%as.numeric()
  l = data[,3]%>%as.numeric()
  c = data[,4]%>%as.numeric()
  k=seq(-10,10,1)
  li=c()
  for (i in 1:length(u)) {
    g1= (4*k*(k+1))/(sqrt(2*pi)*sigma1^3)*(1-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/sigma1^2)
    g2 = (4*k^2)/(sqrt(2*pi)*sigma1^3)*(1-((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/sigma1^2)
    logh1=-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/(2*sigma1^2)-mu1^2/(2*sigma1^2)+(mu1*(c[i]-o[i]))/sigma1^2
    h1=exp(logh1)
    logh2 = -((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/(2*sigma1^2)-mu1^2/(2*sigma1^2)+(mu1*(c[i]-o[i]))/sigma1^2
    h2=exp(logh2)
    li[i]=sum(g1*h1*(1/sigma1^2)-g2*h2*(1/sigma1^2))/sum(g1 * h1 - g2* h2)
  }
  
  return(sum(li))
}


partial_ell_mu2_mu2 <- function(mu2, sigma2, data) {
  o = data[,1]%>%as.numeric()
  u = data[,2]%>%as.numeric()
  l = data[,3]%>%as.numeric()
  c = data[,4]%>%as.numeric()
  k=seq(-10,10,1)
  li=c()
  for (i in 1:length(u)) {
    g1= (4*k*(k+1))/(sqrt(2*pi)*sigma2^3)*(1-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/sigma2^2)
    g2 = (4*k^2)/(sqrt(2*pi)*sigma2^3)*(1-((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/sigma2^2)
    logh1=-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/(2*sigma2^2)-mu2^2/(2*sigma2^2)+(mu2*(c[i]-o[i]))/sigma2^2
    h1=exp(logh1)
    logh2 = -((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/(2*sigma2^2)-mu2^2/(2*sigma2^2)+(mu2*(c[i]-o[i]))/sigma2^2
    h2=exp(logh2)
    li[i]=sum(g1*h1*(1/sigma2^2)-g2*h2*(1/sigma2^2))/sum(g1 * h1 - g2* h2)
  }
  
  return(sum(li))
}


partial_ell_sigma=function(mu1, sigma1, data){
  o = data[,1]%>%as.numeric()
  u = data[,2]%>%as.numeric()
  l = data[,3]%>%as.numeric()
  c = data[,4]%>%as.numeric()
  k=seq(-10,10,1)
  s=sigma1^2
  li=c()
  
  for (i in 1:length(u)) {
    g1= (4*k*(k+1))/(sqrt(2*pi)*sigma1^3)*(1-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/sigma1^2)
    g2 = (4*k^2)/(sqrt(2*pi)*sigma1^3)*(1-((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/sigma1^2)
    logh1=-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/(2*sigma1^2)-mu1^2/(2*sigma1^2)+(mu1*(c[i]-o[i]))/sigma1^2
    h1=exp(logh1)
    logh2 = -((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/(2*sigma1^2)-mu1^2/(2*sigma1^2)+(mu1*(c[i]-o[i]))/sigma1^2
    h2=exp(logh2)
    z=sum(g1 * h1 - g2* h2)
    d_g1=(-3/(2*sigma1^2))*g1+(4*k*(k+1))/sqrt(2*pi)*(c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    
    d_h1=-(h1*logh1)/sigma1^2
    
    d_g2=(-3/(2*sigma1^2))*g2+(4*k^2)/sqrt(2*pi)*(c[i]-o[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    
    d_h2=-(h2*logh2)/sigma1^2
    
    li[i]=(sum(d_g1*h1+d_h1*g1)-sum(d_g2*h2+d_h2*g2))/z
    
  }
  
  li=sum(li)
  return(li)
}

partial_ell_sigma_sigma <- function(mu1, sigma1, data) {
  o = data[,1]%>%as.numeric()
  u = data[,2]%>%as.numeric()
  l = data[,3]%>%as.numeric()
  c = data[,4]%>%as.numeric()
  k=seq(-10,10,1)
  s=sigma1^2
  li=c()
  for (i in 1:length(u)) {
    g1= (4*k*(k+1))/(sqrt(2*pi)*sigma1^3)*(1-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/sigma1^2)
    g2 = (4*k^2)/(sqrt(2*pi)*sigma1^3)*(1-((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/sigma1^2)
    logh1=-((c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2)/(2*sigma1^2)-mu1^2/(2*sigma1^2)+(mu1*(c[i]-o[i]))/sigma1^2
    h1=exp(logh1)
    logh2 = -((c[i]-o[i]-2*k*(u[i]-l[i]))^2)/(2*sigma1^2)-mu1^2/(2*sigma1^2)+(mu1*(c[i]-o[i]))/sigma1^2
    h2=exp(logh2)
    z=sum(g1 * h1 - g2* h2)
    d_g1=(-3/(2*sigma1^2))*g1+(4*k*(k+1))/sqrt(2*pi)*(c[i]+o[i]-2*u[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    dd_g1= -(3 / 2 ) * ((sigma1^2 * d_g1 - g1) / (sigma1^2)^2) -7 / 2 * 4 * k * (k + 1) / sqrt(2 * pi) * ((c[i] + o[i] - 2 * u[i] - 2 * k * (u[i] - l[i]))^2) * (sigma1^2)^(-9 / 2)
    d_h1=-(h1*logh1)/sigma1^2
    dd_h1=(s*(d_h1*logh1+d_h1)-h1*logh1)/(s)^2
    d_g2=(-3/(2*sigma1^2))*g2+(4*k^2)/sqrt(2*pi)*(c[i]-o[i]-2*k*(u[i]-l[i]))^2*s^(-7/2)
    dd_g2= -(3 / 2 ) * ((sigma1^2 * d_g2 - g2) / (sigma1^2)^2) -7 / 2 * 4 * k^2 / sqrt(2 * pi) * ((c[i] - o[i]  - 2 * k * (u[i] - l[i]))^2) * (sigma1^2)^(-9 / 2)
    d_h2=-(h2*logh2)/sigma1^2
    dd_h2=(s*(d_h2*logh2+d_h2)-h2*logh2)/(s)^2
    li[i]=sum(dd_g1*h1+2*d_g1*d_h1+dd_h1*g1)/z-sum(dd_g2*h2+2*d_g2*d_h2+dd_h2*g2)/z-((sum(d_g1*h1+d_h1*g1-d_g2*h2-d_h2*g2))^2)/(z)^2
    
  }
  
  li=sum(li)
  return(li)
}






SPX_2020 <- getSymbols("^GSPC",auto.assign = FALSE, from = "2021-12-31",to="2022-05-20")
n=length(SPX_2020$GSPC.Adjusted)
sig1 <- log(SPX_2020$GSPC.Adjusted[1:22])%>%sd()
sig2 <- log(SPX_2020$GSPC.Adjusted[(n-22):n])%>%sd()

rt = diff(log(SPX_2020$GSPC.Adjusted))
sc.rt=ts(rt*100,frequency = 1)
plot.ts(sc.rt,ylab="log return")
abline(v=201,col="red")
data = log(SPX_2020[,1:4])
plot(as.numeric(data[,2]))



change.point.detect.real = function(sig0,n.check){
  est.matrix = matrix(NA, ncol = 5, nrow = 1)
  colnames(est.matrix) = c("tau","mu1","mu2","sigma21","sigma22")
  
  
  data = log(SPX_2020[,1:4])
  
  n=dim(data)[1]
  
  ll=c()
  for (i in n.check:(n-3)){
    
    
    
    d1 = data[1:i,]
    d2 = data[(i+1):n,]
    a=unlist(lapply(sig0, function(x) point.est(sig0 = x, data = d1)))
    b=unlist(lapply(sig0, function(x) point.est(sig0 = x, data = d2)))
    est1.i = matrix(a,nrow = length(sig0),ncol = 2 ,byrow = T)
    
    est2.i = matrix(b,nrow = length(sig0),ncol = 2,byrow = T)
    
    l1 <- apply(est1.i,1, function(x) l(data = d1, mu = x[1], sig = x[2]))
    
    index1 <- which.max(l1)
    l1 <- max(l1,na.rm = TRUE)
    
    l2 = apply(est2.i, 1, function(x) l(data = d2, mu = x[1], sig = x[2])) 
    index2=which.max(l2)
    
    l2=max(l2,na.rm = TRUE)
    ll[i-(n.check-1)]=l1+l2
    
    
    
    
    
    
  }
  
  
  tau.hat = which.max(ll)+n.check-1
  
  
  a=unlist(lapply(sig0, function(x) point.est(sig0 = x, data = data[1:tau.hat,])))
  b=unlist(lapply(sig0, function(x) point.est(sig0 = x, data = data[(tau.hat+1):n,])))
  est1 = matrix(a,nrow = length(sig0),ncol = 2,byrow = T)
  
  est2 = matrix(b,nrow = length(sig0),ncol = 2,byrow = T)
  
  l1 <- apply(est1,1, function(x) l(data = data[1:tau.hat,], mu = x[1], sig = x[2]))
  
  index1 <- which.max(l1)
  
  l2 = apply(est2, 1, function(x) l(data = data[(tau.hat+1):n,], mu = x[1], sig = x[2])) 
  index2=which.max(l2)
  
  est1.mu=est1[index1,1]
  est1.sig=est1[index1,2]
  est2.mu=est2[index2,1]
  est2.sig=est2[index2,2]
  
  est.matrix = c(tau.hat, est1.mu, est2.mu,  exp(est1.sig)^2, exp(est2.sig)^2)
  
  beep(2)
  
  
  return(list(est.matrix = est.matrix,ll=ll))
}


n.check=3
sig0=log(c(sig1,sig2))
#sig0=log(runif(5))
p1=Sys.time()
result = change.point.detect.real(sig0,n.check)
p2=Sys.time()
p2-p1
options(scipen = 999)
result






