####### Functions for uncoinstrained G-estimation and Bootstrap variance estimation ####
Vec_J_bin <- function(df,theta){
  ## Read in parameters and construct residuals ##
  beta = theta[1:3]
  gamma = theta[-(1:3)]
  beta1 = beta[1]
  beta2 = beta[2]
  delta  = beta[3]
  
  p = length(gamma)/6
  gamma_x1 = gamma[1:p]
  gamma_x2 = gamma[(p+1)  :(2*p)]
  gamma_m1 = gamma[(2*p+1):(3*p)]
  gamma_m2 = gamma[(3*p+1):(4*p)]
  gamma_y1 = gamma[(4*p+1):(5*p)]
  gamma_y2 = gamma[(5*p+1):(6*p)]
  
  X = df$X
  M = df$M
  Z = df$Z
  Y = df$Y
  N = length(M)
  Z.matrix = cbind(1,Z)
  
  odds_x1  = Z.matrix%*%gamma_x1
  odds_x2  = Z.matrix%*%gamma_x2
  f1       = Z.matrix%*%gamma_m1
  f2       = Z.matrix%*%gamma_m2
  g1       = Z.matrix%*%gamma_y1 
  g2       = Z.matrix%*%gamma_y2 
  
  e1 = exp(odds_x1)
  e2 = exp(odds_x2)
  
  #Assume X is binomail with nt trials
  nt=1
  
  x.res1   = X-nt/(1+exp(-1*odds_x1))
  x.res2   = X-nt/(1+exp(-1*odds_x2))
  m.res1   = M - beta1*X - f1 
  m.res2   = M - beta1*X - f2 
  y.res1   = Y - beta2*M - delta*X - g1
  y.res2   = Y - beta2*M - delta*X - g2
  
  hdot1     = nt*e1/(1+e1)^2
  hdot2     = nt*e2/(1+e2)^2
  
  #### Estimating Equations and first derivatives ####
  
  U1 = sum(x.res1*m.res1)/N
  U2 = sum(m.res2*y.res1)/N
  U3 = sum(x.res2*y.res2)/N
  U = c(U1,U2,U3)
  
  V1 = t(Z.matrix) %*% (m.res1*hdot1)/N
  V2 = t(Z.matrix) %*% x.res1/N
  V3 = t(Z.matrix) %*% y.res1/N
  V4 = t(Z.matrix) %*% m.res2/N
  V5 = t(Z.matrix) %*% (y.res2*hdot2)/N
  V6 = t(Z.matrix) %*% x.res2/N
  
  sig11 = sum(x.res1^2*m.res1^2)/N
  sig22 = sum(m.res2^2*y.res1^2)/N
  sig33 = sum(x.res2^2*y.res2^2)/N
  
  sig12 = sum(x.res1*m.res1*m.res2*y.res1)/N
  sig13 = sum(x.res1*m.res1*x.res2*y.res2)/N
  sig23 = sum(m.res2*y.res1*x.res2*y.res2)/N
  
  Sig = matrix(c(sig11,sig12,sig13,
                 sig12,sig22,sig23,
                 sig13,sig23,sig33),ncol=3,byrow = T)#/N
  
  dU1dbeta = -c(sum(X*x.res1),0,0)/N
  dU2dbeta = -c(sum(X*y.res1),sum(M*m.res2),sum(X*m.res2))/N
  dU3dbeta = -c(0,sum(M*x.res2),sum(X*x.res2))/N
  
  p0 = rep.int(0,p)
  dU1dgamma = -c(V1,p0,V2,p0,p0,p0)
  dU2dgamma = -c(p0,p0,p0,V3,V4,p0)
  dU3dgamma = -c(p0,V5,p0,p0,p0,V6)
  
  dUdbeta = matrix(c(dU1dbeta,dU2dbeta,dU3dbeta),nrow=3,byrow = TRUE)
  dUdgamma = matrix(c(dU1dgamma,dU2dgamma,dU3dgamma),nrow=3,byrow = TRUE)
  
  dUdtheta = cbind(dUdbeta,dUdgamma)
  
  #### nuisance parameters ####
  dV1dbeta = - t(Z.matrix) %*%cbind((X*hdot1),0,0)
  dV2dbeta = -                cbind(p0       ,0,0)
  dV3dbeta = - t(Z.matrix) %*%cbind(0        ,M,X)
  dV4dbeta = - t(Z.matrix) %*%cbind(X        ,0,0)
  dV5dbeta = - t(Z.matrix) %*%cbind(0        ,(M*hdot2),(X*hdot2))
  dV6dbeta = -                cbind(p0       ,0,0)
  
  dVdbeta = matrix(rbind(dV1dbeta,dV2dbeta,dV3dbeta,dV4dbeta,dV5dbeta,dV6dbeta),ncol=3)/N
  
  pp0 = matrix(0,nrow=p,ncol=p)
  
  Z.sq = t(Z.matrix) %*%Z.matrix
  Z.hdot1 = t(Z.matrix) %*%(Z.matrix*c(hdot1))
  Z.hdot2 = t(Z.matrix) %*%(Z.matrix*c(hdot1))
  
  Z.hdotdot1 = -t(Z.matrix) %*% (Z.matrix*c(-((hdot1*(e1-1)/(e1+1))*m.res1)))
  Z.hdotdot2 = -t(Z.matrix) %*% (Z.matrix*c(-((hdot2*(e2-1)/(e2+1))*y.res2)))
  
  dV1dgamma = -cbind(Z.hdotdot1,pp0,Z.hdot1,pp0,pp0,pp0)
  dV2dgamma = -cbind(Z.hdot1,pp0,pp0,pp0,pp0,pp0)
  dV3dgamma = -cbind(pp0,pp0,pp0,pp0,Z.sq,pp0)
  dV4dgamma = -cbind(pp0,pp0,pp0,Z.sq,pp0,pp0)
  dV5dgamma = -cbind(pp0,Z.hdotdot2,pp0,pp0,pp0,Z.hdot2)
  dV6dgamma = -cbind(pp0,Z.hdot2,pp0,pp0,pp0,pp0)
  
  dVdgamma = matrix(rbind(dV1dgamma,dV2dgamma,dV3dgamma,dV4dgamma,dV5dgamma,dV6dgamma),ncol=6*p)/N
  
  
  #### Return output ####
  vec = c(U,V1,V2,V3,V4,V5,V6)
  J   = rbind(cbind(dUdbeta,dUdgamma),cbind(dVdbeta,dVdgamma))
  
  inv_dUdbeta = solve(dUdbeta)
  PhiPhiT = inv_dUdbeta%*%Sig%*%t(inv_dUdbeta)
  
  return(list(vec=vec,J=J,
              var=diag(PhiPhiT)/N))
}


fit_beta <-   function(df){
  theta =  MLE_thetas_bin(df)$theta
  Max.it = 1000
  prec = 1e-12
  converged = FALSE
  i = 1
  while(i<Max.it& !converged){
    step = Vec_J_bin(df,theta)
    theta = theta - solve(step$J,step$vec)
    if(sum((abs(step$vec)>=prec))==0){converged = TRUE}
    i = i+1
  }
  var.hat = Vec_J_bin(df,theta)$var
  nide.hat = theta[1]*theta[2]
  nide.var = theta[1]^2*var.hat[2] + theta[2]^2*var.hat[1]
  return(c(theta[1:3],nide.hat,var.hat,nide.var))
}


boot_beta <- function(df,N.boot=1000){
  N = dim(df)[1]
  boot_reps = t(replicate(N.boot,{fit_beta(df[sample.int(N,replace = TRUE),])}))
  return(apply(boot_reps[,1:4],2,var))
}  

