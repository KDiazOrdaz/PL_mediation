#functions for the minimum of two chisquared rvs
qmin_chisq <- function(p,lower.tail = TRUE){
  if (!lower.tail){ p = 1-p}
  x = qchisq(1-sqrt(1-p),df=1)
  return(x)
}

pmin_chisq <-function(q,lower.tail = TRUE){
  x = 1 - (1 - pchisq(q,df=1))^2
  if (!lower.tail){ x = 1-x}
  return(x)
}

dmin_chisq <-function(x){
  2*(1-pchisq(x,df=1))*dchisq(x,df=1)
}

rmin_chisq <- function(n){
  apply(matrix(rchisq(2*n,df=1),ncol=2),1,min)
}

#and for the max of two chisquared rvs
qmax_chisq <- function(p,lower.tail = TRUE){
  if (!lower.tail){ p = 1-p}
  x = qchisq(sqrt(p),df=1)
  return(x)
}

pmax_chisq <-function(q,lower.tail = TRUE){
  x = pchisq(q,df=1)^2
  if (!lower.tail){ x = 1-x}
  return(x)
}

dmax_chisq <-function(x){
  2*pchisq(x,df=1)*dchisq(x,df=1)
}

rmax_chisq <- function(n){
  apply(matrix(rchisq(2*n,df=1),ncol=2),1,max)
}

## For the product normal distribution
dprod_norm <- function(x){
  besselK(abs(x),0)/pi
}

pprod_norm <- function(q){
  require(RandomFieldsUtils) #required for struvel function
  Int_pdf <- function(z) {
    ifelse(z==0,0,(z/2)*(besselK(z,0)*struveL(z,-1) +besselK(z,1)*struveL(z,0)))
    }
  #Int_pdf is the indefinite integral of dprod_norm (without the abs value bit)
  #Source: http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/21/01/01/0004/
  #tends to zero at zero
  #tends to 1/2 at infinity
  (1/2) + sign(q)*Int_pdf(abs(q))
}

qprod_norm <- function(p){
  #invert pprodnorm with newton raphson
  x = qnorm(p)
  for (i in 1:10){#enough precision for 8 decimal places
    x = x - (pprod_norm(x) - p)/dprod_norm(x)
  }
  (x)
}

rprod_norm <- function(n){
  rnorm(n)*rnorm(n)
}

## some useful plotting functions
PP.plot <- function(set,pfunc,...){#set of data, quantile function, arguments of qfunc
  z = sort(set)
  n = length(z)
  x = pfunc(z ,...)
  Fn = ecdf(z)
  y = Fn(z)
  
  plot(x,y,main='PP-plot',
       xlab='Theoretical',ylab='Observed')
  abline(0,1,col='red')
}

QQ.plot <- function(set,qfunc,...){#set of data, quantile function, arguments of qfunc
  y = sort(set)
  n = length(y)
  x = qfunc((1:n)/(n+1) ,...)
  plot(x,y,main='QQ-plot',
       xlab='Theoretical',ylab='Observed')
  abline(0,1,col='red')
}

QQ.plot_g <- function(set,qfunc,...){#set of data, quantile function, arguments of qfunc
  y = sort(set)
  n = length(y)
  x = qfunc((1:n)/(n+1) ,...)
  df = data.frame(Theoretical = x, Observed = y)
  require(ggplot2)
  
  p = ggplot(data=df, mapping=aes(x=Theoretical,y=Observed)) + 
    geom_point() + geom_abline(slope=1,intercept = 0,col='red')+
    theme_classic()
  
  return(p)
}


## MC estimates of KL divergences

KL_divergence_chisq <- function(dfunc,N.mc,...){
  x = rchisq(N.mc,1)
  mean(log(dchisq(x,1)/dfunc(x,...)))
}

KL_divergence_norm <- function(dfunc,N.mc,...){
  x = rnorm(N.mc)
  mean(log(dnorm(x)/dfunc(x,...)))
}





